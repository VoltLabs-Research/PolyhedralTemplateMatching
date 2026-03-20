#include <volt/analysis/structure_analysis.h>
#include <volt/polyhedral_template_matching.h>

#include <ptm_constants.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

#include <cmath>
#include <limits>
#include <vector>

namespace Volt {

namespace {

bool setupPTM(AnalysisContext& context, Volt::PTM& ptm, size_t particleCount){
    ptm.setCalculateDefGradient(false);
    ptm.setRmsdCutoff(std::numeric_limits<double>::infinity());
    return ptm.prepare(context.positions->constDataPoint3(), particleCount, context.simCell);
}

void storeOrientationData(AnalysisContext& context, const PTM::Kernel& kernel, size_t atomIndex){
    auto quaternion = kernel.orientation();
    double* orientation = context.ptmOrientation->dataDouble() + 4 * atomIndex;

    orientation[0] = quaternion.x();
    orientation[1] = quaternion.y();
    orientation[2] = quaternion.z();
    orientation[3] = quaternion.w();
}

void storeDeformationGradient(AnalysisContext& context, const PTM::Kernel& kernel, size_t atomIndex) {
    if(!context.ptmDeformationGradient){
        return;
    }

    const auto& F = kernel.deformationGradient();
    double* destination = context.ptmDeformationGradient->dataDouble() + 9 * atomIndex;
    const double* source = F.elements();

    for(int k = 0; k < 9; ++k){
        destination[k] = source[k];
    }
}

}

void StructureAnalysis::computeMaximumNeighborDistanceFromPTM(){
    const size_t N = _context.atomCount();
    if(N == 0){
        _maximumNeighborDistance = 0.0;
        return;
    }

    const auto* positions = _context.positions->constDataPoint3();
    const auto& inverseMatrix = _context.simCell.inverseMatrix();
    const auto& directMatrix = _context.simCell.matrix();
    const int* counts = _context.neighborCounts->constDataInt();
    const int* offsets = _context.neighborOffsets->constDataInt();
    const int* indices = _context.neighborIndices->constDataInt();

    double maxDistance = tbb::parallel_reduce(
        tbb::blocked_range<size_t>(0, N),
        0.0,
        [&](const tbb::blocked_range<size_t>& range, double maxSoFar) -> double {
            for(size_t i = range.begin(); i < range.end(); ++i){
                const int neighborCount = counts[i];
                double localMaxDistance = 0.0;
                const int start = offsets[i];
                for(int j = 0; j < neighborCount; ++j){
                    int neighbor = indices[start + j];

                    Vector3 delta = positions[neighbor] - positions[i];
                    double f[3];
                    for(int d = 0; d < 3; ++d){
                        f[d] = inverseMatrix.prodrow(delta, d);
                        f[d] -= std::round(f[d]);
                    }

                    Vector3 minimumImage;
                    minimumImage = directMatrix.column(0) * f[0];
                    minimumImage += directMatrix.column(1) * f[1];
                    minimumImage += directMatrix.column(2) * f[2];
                    double distance = minimumImage.length();
                    if(distance > localMaxDistance){
                        localMaxDistance = distance;
                    }
                }
                if(localMaxDistance > maxSoFar){
                    maxSoFar = localMaxDistance;
                }
            }
            return maxSoFar;
        },
        [](double a, double b) -> double { return std::max(a, b); }
    );

    spdlog::debug("Maximum neighbor distance (from PTM): {}", maxDistance);
    _maximumNeighborDistance = maxDistance;
}

void StructureAnalysis::determineLocalStructuresWithPTM() {
    const size_t N = _context.atomCount();
    if(!N){
        return;
    }

    Volt::PTM ptm;
    if(!setupPTM(_context, ptm, N)){
        throw std::runtime_error("Error trying to initialize PTM.");
    }

    _context.ptmOrientation = std::make_shared<ParticleProperty>(N, DataType::Double, 4, 0.0, true);
    _context.ptmRmsd = std::make_shared<ParticleProperty>(N, DataType::Double, 1, 0.0, true);
    _context.correspondencesCode = std::make_shared<ParticleProperty>(N, DataType::Int64, 1, 0, true);
    if(ptm.calculateDefGradient()){
        _context.ptmDeformationGradient = std::make_shared<ParticleProperty>(N, DataType::Double, 9, 0.0, true);
    }else{
        _context.ptmDeformationGradient.reset();
    }

    std::fill(_context.neighborCounts->dataInt(),
              _context.neighborCounts->dataInt() + _context.neighborCounts->size(), 0);
    std::fill(_context.structureTypes->dataInt(),
              _context.structureTypes->dataInt() + _context.structureTypes->size(), LATTICE_OTHER);

    std::vector<uint64_t> cached(N, 0ull);
    std::vector<int> localCounts(N, 0);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, N), [&](const auto &range){
        PTM::Kernel kernel(ptm);
        for(size_t i = range.begin(); i < range.end(); ++i){
            kernel.cacheNeighbors(i, &cached[i]);
            StructureType type = kernel.identifyStructure(i, cached);

            const double rmsd = kernel.rmsd();
            _context.ptmRmsd->setDouble(i, rmsd);

            auto* correspondences = reinterpret_cast<uint64_t*>(_context.correspondencesCode->data());
            correspondences[i] = kernel.correspondencesCode();

            if(type == StructureType::OTHER || rmsd > _rmsd){
                continue;
            }

            _context.structureTypes->setInt(i, type);
            const int neighborCount = kernel.numTemplateNeighbors();
            localCounts[i] = neighborCount;
            _context.neighborCounts->setInt(i, neighborCount);
            storeOrientationData(_context, kernel, i);
            storeDeformationGradient(_context, kernel, i);
            _context.templateIndex->setInt(i, kernel.bestTemplateIndex());
        }
    });

    auto* offsets = _context.neighborOffsets->dataInt();
    offsets[0] = 0;
    for(size_t i = 0; i < N; ++i){
        offsets[i + 1] = offsets[i] + localCounts[i];
    }

    const size_t totalNeighbors = static_cast<size_t>(offsets[N]);
    _context.neighborIndices = std::make_shared<ParticleProperty>(
        totalNeighbors, DataType::Int, 1, 0, false);
    auto* indices = _context.neighborIndices->dataInt();

    tbb::parallel_for(tbb::blocked_range<size_t>(0, N), [&](const auto& range){
        PTM::Kernel kernel(ptm);
        for(size_t i = range.begin(); i < range.end(); ++i){
            const int count = localCounts[i];
            if(count == 0){
                continue;
            }
            kernel.cacheNeighbors(i, &cached[i]);
            kernel.identifyStructure(i, cached);
            const int start = offsets[i];
            for(int j = 0; j < count; ++j){
                indices[start + j] = kernel.getTemplateNeighbor(j).index;
            }
        }
    });

    for(size_t i = 0; i < N; ++i){
        if(_context.neighborCounts->getInt(i) == 0){
            _context.neighborCounts->setInt(i, localCounts[i]);
        }
    }
}

}
