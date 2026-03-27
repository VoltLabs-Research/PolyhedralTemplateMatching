#include <volt/analysis/structure_analysis.h>
#include <volt/polyhedral_template_matching.h>

#include <ptm_constants.h>
#include <ptm_initialize_data.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

#include <cmath>
#include <array>
#include <memory>
#include <limits>
#include <vector>

namespace Volt {

constexpr std::array<int, 6> kSimpleCubicTemplateToCanonicalNeighborSlot = {5, 4, 3, 2, 1, 0};
constexpr std::array<int, 6> kSimpleCubicCanonicalToTemplateSlot = {5, 4, 3, 2, 1, 0};

bool setupPTM(StructureContext& context, Volt::PTM& ptm, size_t particleCount){
    ptm.setCalculateDefGradient(false);
    ptm.setRmsdCutoff(std::numeric_limits<double>::infinity());
    return ptm.prepare(context.positions->constDataPoint3(), particleCount, context.simCell);
}

namespace{

struct PtmSymmetryPermutation{
    std::vector<int> permutation;
    Matrix3 transformation = Matrix3::Identity();
    std::vector<int> inverseProduct;
};

struct PtmCrystalData{
    int coordinationNumber = 0;
    std::vector<Vector3> latticeVectors;
    std::vector<std::array<int, 2>> commonNeighbors;
    std::vector<PtmSymmetryPermutation> symmetries;
};

int normalizedStructureType(int structureType){
    switch(static_cast<StructureType>(structureType)){
        case StructureType::CUBIC_DIAMOND_FIRST_NEIGH:
        case StructureType::CUBIC_DIAMOND_SECOND_NEIGH:
            return static_cast<int>(StructureType::CUBIC_DIAMOND);
        case StructureType::HEX_DIAMOND_FIRST_NEIGH:
        case StructureType::HEX_DIAMOND_SECOND_NEIGH:
            return static_cast<int>(StructureType::HEX_DIAMOND);
        default:
            return structureType;
    }
}

const ptm::refdata_t* ptmReferenceForStructureType(int structureType){
    const StructureType normalized = static_cast<StructureType>(normalizedStructureType(structureType));
    if(normalized == StructureType::OTHER){
        return nullptr;
    }
    const int ptmType = PTM::toPtmStructureType(normalized);
    if(ptmType == PTM_MATCH_NONE){
        return nullptr;
    }
    return ptm::refdata[ptmType];
}

int canonicalTemplateIndex(const ptm::refdata_t& ref){
    if(ref.num_conventional_mappings > 0 && ref.template_indices != nullptr){
        return ref.template_indices[0];
    }
    return 0;
}

const int8_t (*symmetryMappings(const ptm::refdata_t& ref))[PTM_MAX_POINTS]{
    if(ref.num_conventional_mappings > 0 && ref.mapping_conventional != nullptr){
        return ref.mapping_conventional;
    }
    return ref.mapping;
}

int symmetryMappingCount(const ptm::refdata_t& ref){
    if(ref.num_conventional_mappings > 0 && ref.mapping_conventional != nullptr){
        return ref.num_conventional_mappings;
    }
    return ref.num_mappings;
}

std::vector<Vector3> buildCanonicalLatticeVectors(int structureType, const ptm::refdata_t& ref){
    const int numNeighbors = ref.num_nbrs;
    std::vector<Vector3> latticeVectors(static_cast<std::size_t>(numNeighbors), Vector3::Zero());
    const double (*points)[3] = ref.points[canonicalTemplateIndex(ref)];

    if(static_cast<StructureType>(normalizedStructureType(structureType)) == StructureType::SC){
        for(int templateSlot = 0; templateSlot < numNeighbors; ++templateSlot){
            const int canonicalSlot = kSimpleCubicTemplateToCanonicalNeighborSlot[templateSlot];
            const double* point = points[templateSlot + 1];
            latticeVectors[canonicalSlot] = Vector3(point[0], point[1], point[2]);
        }
        return latticeVectors;
    }

    for(int neighborSlot = 0; neighborSlot < numNeighbors; ++neighborSlot){
        const double* point = points[neighborSlot + 1];
        latticeVectors[static_cast<std::size_t>(neighborSlot)] = Vector3(point[0], point[1], point[2]);
    }
    return latticeVectors;
}

std::vector<int> buildSymmetryPermutation(
    int structureType,
    const ptm::refdata_t& ref,
    int mappingIndex
){
    const int numNeighbors = ref.num_nbrs;
    const int8_t* mapping = symmetryMappings(ref)[mappingIndex];
    std::vector<int> permutation(static_cast<std::size_t>(numNeighbors), 0);

    if(static_cast<StructureType>(normalizedStructureType(structureType)) == StructureType::SC){
        for(int canonicalSlot = 0; canonicalSlot < numNeighbors; ++canonicalSlot){
            const int templateSlot = kSimpleCubicCanonicalToTemplateSlot[canonicalSlot];
            const int mappedTemplateSlot = mapping[templateSlot + 1] - 1;
            permutation[static_cast<std::size_t>(canonicalSlot)] =
                kSimpleCubicTemplateToCanonicalNeighborSlot[static_cast<std::size_t>(mappedTemplateSlot)];
        }
        return permutation;
    }

    for(int neighborSlot = 0; neighborSlot < numNeighbors; ++neighborSlot){
        permutation[static_cast<std::size_t>(neighborSlot)] = mapping[neighborSlot + 1] - 1;
    }
    return permutation;
}

std::array<int, 3> findNonCoplanarIndices(const std::vector<Vector3>& latticeVectors){
    std::array<int, 3> indices{-1, -1, -1};
    Matrix3 basis = Matrix3::Zero();
    int found = 0;

    for(int vectorIndex = 0; vectorIndex < static_cast<int>(latticeVectors.size()) && found < 3; ++vectorIndex){
        basis.column(found) = latticeVectors[static_cast<std::size_t>(vectorIndex)];

        if(found == 1){
            if(basis.column(0).cross(basis.column(1)).squaredLength() <= EPSILON){
                continue;
            }
        }else if(found == 2){
            if(std::abs(basis.determinant()) <= EPSILON){
                continue;
            }
        }

        indices[static_cast<std::size_t>(found++)] = vectorIndex;
    }

    if(found != 3){
        throw std::runtime_error("Unable to determine a non-coplanar PTM basis.");
    }

    return indices;
}

std::vector<std::array<int, 2>> buildCommonNeighbors(
    int structureType,
    const std::vector<Vector3>& latticeVectors
){
    const int numNeighbors = static_cast<int>(latticeVectors.size());
    std::vector<std::array<int, 2>> commonNeighbors(
        static_cast<std::size_t>(numNeighbors),
        std::array<int, 2>{-1, -1}
    );

    for(int neighborIndex = 0; neighborIndex < numNeighbors; ++neighborIndex){
        Matrix3 basis = Matrix3::Zero();
        basis.column(0) = latticeVectors[static_cast<std::size_t>(neighborIndex)];

        if(static_cast<StructureType>(normalizedStructureType(structureType)) == StructureType::SC){
            for(int i1 = 0; i1 < numNeighbors; ++i1){
                if(i1 == neighborIndex){
                    continue;
                }
                basis.column(1) = latticeVectors[static_cast<std::size_t>(i1)];
                for(int i2 = i1 + 1; i2 < numNeighbors; ++i2){
                    if(i2 == neighborIndex){
                        continue;
                    }
                    basis.column(2) = latticeVectors[static_cast<std::size_t>(i2)];
                    if(std::abs(basis.determinant()) > EPSILON){
                        commonNeighbors[static_cast<std::size_t>(neighborIndex)] = {i1, i2};
                        goto next_neighbor;
                    }
                }
            }
        }

        {
            double minBondDistanceSq = std::numeric_limits<double>::max();
            for(int otherIndex = 0; otherIndex < numNeighbors; ++otherIndex){
                if(otherIndex == neighborIndex){
                    continue;
                }
                const double distSq = (
                    latticeVectors[static_cast<std::size_t>(neighborIndex)] -
                    latticeVectors[static_cast<std::size_t>(otherIndex)]
                ).squaredLength();
                if(distSq > EPSILON && distSq < minBondDistanceSq){
                    minBondDistanceSq = distSq;
                }
            }

            if(!std::isfinite(minBondDistanceSq)){
                goto next_neighbor;
            }

            for(int i1 = 0; i1 < numNeighbors; ++i1){
                if(i1 == neighborIndex){
                    continue;
                }
                const double d1 = (
                    latticeVectors[static_cast<std::size_t>(neighborIndex)] -
                    latticeVectors[static_cast<std::size_t>(i1)]
                ).squaredLength();
                if(std::abs(d1 - minBondDistanceSq) > 1e-6){
                    continue;
                }
                basis.column(1) = latticeVectors[static_cast<std::size_t>(i1)];

                for(int i2 = i1 + 1; i2 < numNeighbors; ++i2){
                    if(i2 == neighborIndex){
                        continue;
                    }
                    const double d2 = (
                        latticeVectors[static_cast<std::size_t>(neighborIndex)] -
                        latticeVectors[static_cast<std::size_t>(i2)]
                    ).squaredLength();
                    if(std::abs(d2 - minBondDistanceSq) > 1e-6){
                        continue;
                    }
                    basis.column(2) = latticeVectors[static_cast<std::size_t>(i2)];
                    if(std::abs(basis.determinant()) > EPSILON){
                        commonNeighbors[static_cast<std::size_t>(neighborIndex)] = {i1, i2};
                        goto next_neighbor;
                    }
                }
            }
        }

    next_neighbor:
        continue;
    }

    return commonNeighbors;
}

PtmCrystalData buildCrystalData(int structureType){
    PtmCrystalData data;
    const ptm::refdata_t* ref = ptmReferenceForStructureType(structureType);
    if(!ref){
        return data;
    }

    data.coordinationNumber = ref->num_nbrs;
    data.latticeVectors = buildCanonicalLatticeVectors(structureType, *ref);
    data.commonNeighbors = buildCommonNeighbors(structureType, data.latticeVectors);

    const auto basisIndices = findNonCoplanarIndices(data.latticeVectors);
    Matrix3 basis = Matrix3::Zero();
    basis.column(0) = data.latticeVectors[static_cast<std::size_t>(basisIndices[0])];
    basis.column(1) = data.latticeVectors[static_cast<std::size_t>(basisIndices[1])];
    basis.column(2) = data.latticeVectors[static_cast<std::size_t>(basisIndices[2])];
    Matrix3 basisInverse = basis.inverse();

    const int count = symmetryMappingCount(*ref);
    data.symmetries.reserve(static_cast<std::size_t>(count));
    for(int mappingIndex = 0; mappingIndex < count; ++mappingIndex){
        PtmSymmetryPermutation symmetry;
        symmetry.permutation = buildSymmetryPermutation(structureType, *ref, mappingIndex);
        symmetry.transformation = Matrix3::Zero();
        symmetry.transformation.column(0) = data.latticeVectors[
            static_cast<std::size_t>(symmetry.permutation[static_cast<std::size_t>(basisIndices[0])])
        ];
        symmetry.transformation.column(1) = data.latticeVectors[
            static_cast<std::size_t>(symmetry.permutation[static_cast<std::size_t>(basisIndices[1])])
        ];
        symmetry.transformation.column(2) = data.latticeVectors[
            static_cast<std::size_t>(symmetry.permutation[static_cast<std::size_t>(basisIndices[2])])
        ];
        symmetry.transformation = symmetry.transformation * basisInverse;

        bool duplicate = false;
        for(const auto& existing : data.symmetries){
            if(existing.transformation.equals(symmetry.transformation)){
                duplicate = true;
                break;
            }
        }

        if(!duplicate){
            data.symmetries.push_back(std::move(symmetry));
        }
    }

    for(std::size_t s1 = 0; s1 < data.symmetries.size(); ++s1){
        data.symmetries[s1].inverseProduct.reserve(data.symmetries.size());
        for(std::size_t s2 = 0; s2 < data.symmetries.size(); ++s2){
            const Matrix3 inverseProduct =
                data.symmetries[s2].transformation.inverse() *
                data.symmetries[s1].transformation;
            int matchIndex = 0;
            for(std::size_t candidate = 0; candidate < data.symmetries.size(); ++candidate){
                if(data.symmetries[candidate].transformation.equals(inverseProduct)){
                    matchIndex = static_cast<int>(candidate);
                    break;
                }
            }
            data.symmetries[s1].inverseProduct.push_back(matchIndex);
        }
    }

    return data;
}

class PtmCrystalInfoProvider final : public StructureAnalysisCrystalInfo{
public:
    PtmCrystalInfoProvider(){
        initialize(static_cast<int>(StructureType::SC));
        initialize(static_cast<int>(StructureType::FCC));
        initialize(static_cast<int>(StructureType::HCP));
        initialize(static_cast<int>(StructureType::BCC));
        initialize(static_cast<int>(StructureType::ICO));
        initialize(static_cast<int>(StructureType::CUBIC_DIAMOND));
        initialize(static_cast<int>(StructureType::HEX_DIAMOND));
        initialize(static_cast<int>(StructureType::GRAPHENE));
    }

    int findClosestSymmetryPermutation(int structureType, const Matrix3& rotation) const override{
        const auto& symmetries = dataFor(structureType).symmetries;
        int bestIndex = 0;
        double bestDeviation = std::numeric_limits<double>::max();

        for(std::size_t symmetryIndex = 0; symmetryIndex < symmetries.size(); ++symmetryIndex){
            const Matrix3& symmetry = symmetries[symmetryIndex].transformation;
            double deviation = 0.0;
            for(int row = 0; row < 3; ++row){
                for(int column = 0; column < 3; ++column){
                    const double diff = rotation(row, column) - symmetry(row, column);
                    deviation += diff * diff;
                }
            }
            if(deviation < bestDeviation){
                bestDeviation = deviation;
                bestIndex = static_cast<int>(symmetryIndex);
            }
        }

        return bestIndex;
    }

    int coordinationNumber(int structureType) const override{
        return dataFor(structureType).coordinationNumber;
    }

    int commonNeighborIndex(int structureType, int neighborIndex, int commonNeighborSlot) const override{
        return dataFor(structureType).commonNeighbors[static_cast<std::size_t>(neighborIndex)]
            [static_cast<std::size_t>(commonNeighborSlot)];
    }

    int symmetryPermutationCount(int structureType) const override{
        return static_cast<int>(dataFor(structureType).symmetries.size());
    }

    int symmetryPermutationEntry(int structureType, int symmetryIndex, int neighborIndex) const override{
        return dataFor(structureType).symmetries[static_cast<std::size_t>(symmetryIndex)]
            .permutation[static_cast<std::size_t>(neighborIndex)];
    }

    const Matrix3& symmetryTransformation(int structureType, int symmetryIndex) const override{
        return dataFor(structureType).symmetries[static_cast<std::size_t>(symmetryIndex)].transformation;
    }

    int symmetryInverseProduct(int structureType, int symmetryIndex, int transformationIndex) const override{
        return dataFor(structureType).symmetries[static_cast<std::size_t>(symmetryIndex)]
            .inverseProduct[static_cast<std::size_t>(transformationIndex)];
    }

    const Vector3& latticeVector(int structureType, int latticeVectorIndex) const override{
        return dataFor(structureType).latticeVectors[static_cast<std::size_t>(latticeVectorIndex)];
    }

private:
    void initialize(int structureType){
        try{
            _data[static_cast<std::size_t>(normalizedStructureType(structureType))] =
                buildCrystalData(structureType);
        }catch(const std::exception&){
            _data[static_cast<std::size_t>(normalizedStructureType(structureType))] = PtmCrystalData{};
        }
    }

    const PtmCrystalData& dataFor(int structureType) const{
        return _data[static_cast<std::size_t>(normalizedStructureType(structureType))];
    }

    std::array<PtmCrystalData, static_cast<std::size_t>(StructureType::NUM_STRUCTURE_TYPES)> _data{};
};

std::shared_ptr<const StructureAnalysisCrystalInfo> ptmCrystalInfoProvider(){
    static const auto provider = std::make_shared<PtmCrystalInfoProvider>();
    return provider;
}

}

void StructureAnalysis::computeMaximumNeighborDistanceFromPTM(){
    const size_t N = _context.atomCount();
    if(N == 0){
        _context.maximumNeighborDistance = 0.0;
        return;
    }

    const auto* positions = _context.positions->constDataPoint3();
    const auto& inverseMatrix = _context.simCell.inverseMatrix();
    const auto& directMatrix = _context.simCell.matrix();
    const int* counts = _context.neighborCounts->constDataInt();
    const int* offsets = _context.neighborOffsets->constDataInt();
    const int* indices = _context.neighborIndices->constDataInt();

    double maxDistance = tbb::parallel_reduce(tbb::blocked_range<size_t>(0, N),  0.0,
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

    spdlog::info("Maximum neighbor distance (from PTM): {}", maxDistance);
    _context.maximumNeighborDistance = maxDistance;
}

PTMComputationData StructureAnalysis::determineLocalStructuresWithPTM(double rmsdCutoff) {
    const size_t N = _context.atomCount();
    if(!N){
        return {
            std::make_shared<ParticleProperty>(0, DataType::Double, 1, 0.0, true),
            std::make_shared<ParticleProperty>(0, DataType::Double, 4, 0.0, true),
            std::make_shared<ParticleProperty>(0, DataType::Int64, 1, 0, true)
        };
    }
    setCrystalInfoProvider(ptmCrystalInfoProvider());

    Volt::PTM ptm;
    if(!setupPTM(_context, ptm, N)){
        throw std::runtime_error("Error trying to initialize PTM.");
    }

    auto orientations = std::make_shared<ParticleProperty>(N, DataType::Double, 4, 0.0, true);
    auto ptmRmsd = std::make_shared<ParticleProperty>(N, DataType::Double, 1, 0.0, true);
    auto correspondences = std::make_shared<ParticleProperty>(N, DataType::Int64, 1, 0, true);

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
            ptmRmsd->setDouble(i, rmsd);

            auto* encodedCorrespondences = reinterpret_cast<uint64_t*>(correspondences->data());
            encodedCorrespondences[i] = kernel.correspondencesCode();

            if(type == StructureType::OTHER || rmsd > rmsdCutoff){
                continue;
            }

            _context.structureTypes->setInt(i, type);

            const int neighborCount = kernel.numTemplateNeighbors();
            localCounts[i] = neighborCount;
            _context.neighborCounts->setInt(i, neighborCount);
            auto quaternion = kernel.orientation();
            double* orientation = orientations->dataDouble() + 4 * i;

            orientation[0] = quaternion.x();
            orientation[1] = quaternion.y();
            orientation[2] = quaternion.z();
            orientation[3] = quaternion.w();
        }
    });

    auto* offsets = _context.neighborOffsets->dataInt();
    offsets[0] = 0;
    for(size_t i = 0; i < N; ++i){
        offsets[i + 1] = offsets[i] + localCounts[i];
    }

    const size_t totalNeighbors = static_cast<size_t>(offsets[N]);
    _context.neighborIndices = std::make_shared<ParticleProperty>(totalNeighbors, DataType::Int, 1, 0, false);
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
            
            if(_context.structureTypes->getInt(i) == StructureType::SC && count == static_cast<int>(kSimpleCubicTemplateToCanonicalNeighborSlot.size())) {
                for(int templateSlot = 0; templateSlot < count; ++templateSlot){
                    const int canonicalSlot = kSimpleCubicTemplateToCanonicalNeighborSlot[templateSlot];
                    indices[start + canonicalSlot] = kernel.getTemplateNeighbor(templateSlot).index;
                }
            } else {
                for(int j = 0; j < count; ++j){
                    indices[start + j] = kernel.getTemplateNeighbor(j).index;
                }
            }
        }
    });

    for(size_t i = 0; i < N; ++i){
        if(_context.neighborCounts->getInt(i) == 0){
            _context.neighborCounts->setInt(i, localCounts[i]);
        }
    }

    return {ptmRmsd, orientations, correspondences};
}

}
