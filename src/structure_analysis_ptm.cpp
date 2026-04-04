#include <volt/coordination_structures.h>
#include <volt/analysis/nearest_neighbor_finder.h>

#include "ptm_structure_analysis_detail.h"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

#include <cmath>
#include <array>
#include <memory>
#include <limits>
#include <mutex>
#include <vector>

namespace Volt::PtmStructureAnalysisDetail {

bool setupPTM(StructureContext& context, Volt::PTM& ptm, size_t particleCount){
    ptm.setCalculateDefGradient(false);
    ptm.setRmsdCutoff(std::numeric_limits<double>::infinity());
    ptm.setInputCrystalStructure(context.inputCrystalType);
    return ptm.prepare(context.positions->constDataPoint3(), particleCount, context.simCell);
}

std::uint64_t fullSymmetryMask(int symmetryCount){
    if(symmetryCount <= 0){
        return 0;
    }
    if(symmetryCount >= 63){
        return std::numeric_limits<std::uint64_t>::max();
    }
    return (std::uint64_t{1} << symmetryCount) - 1;
}

void ensureCoordinationStructuresInitialized(){
    static std::once_flag initFlag;
    std::call_once(initFlag, []() {
        CoordinationStructures::initializeStructures();
    });
}

PtmCrystalData toPtmCrystalData(const AdaptedCrystalTopology& adapted){
    PtmCrystalData data;
    data.coordinationNumber = adapted.coordinationNumber;
    data.latticeVectors = adapted.latticeVectors;
    data.commonNeighbors.assign(
        adapted.commonNeighbors.begin(),
        adapted.commonNeighbors.begin() + adapted.coordinationNumber
    );
    data.symmetries.reserve(adapted.symmetries.size());
    for(const auto& symmetry : adapted.symmetries){
        PtmSymmetryPermutation ptmSymmetry;
        ptmSymmetry.transformation = symmetry.transformation;
        ptmSymmetry.inverseProduct = symmetry.inverseProduct;
        ptmSymmetry.permutation.assign(
            symmetry.permutation.begin(),
            symmetry.permutation.begin() + adapted.coordinationNumber
        );
        data.symmetries.push_back(std::move(ptmSymmetry));
    }
    return data;
}

PtmCrystalData buildSharedAdaptedPtmCrystalData(int structureType){
    PtmCrystalData data;
    const SharedCrystalTopology* shared = sharedCrystalTopology(structureType);
    const ptm::refdata_t* ref = ptmReferenceForStructureType(structureType);
    if(!shared || !ref){
        return data;
    }

    AdaptedCrystalTopology adapted;
    if(!adaptSharedCrystalTopology(*shared, buildCanonicalLatticeVectors(structureType, *ref), adapted)){
        return data;
    }

    return toPtmCrystalData(adapted);
}

PtmCrystalData buildCnaOrderedDiamondCrystalData(int structureType){
    ensureCoordinationStructuresInitialized();

    PtmCrystalData data;
    const CoordinationStructure& coord = CoordinationStructures::getCoordStruct(structureType);
    const LatticeStructure& lattice = CoordinationStructures::getLatticeStruct(structureType);
    data.coordinationNumber = coord.numNeighbors;
    data.latticeVectors = lattice.latticeVectors;
    data.commonNeighbors.resize(static_cast<std::size_t>(coord.numNeighbors), std::array<int, 2>{-1, -1});
    for(int neighborIndex = 0; neighborIndex < coord.numNeighbors; ++neighborIndex){
        data.commonNeighbors[static_cast<std::size_t>(neighborIndex)] = {
            coord.commonNeighbors[neighborIndex][0],
            coord.commonNeighbors[neighborIndex][1]
        };
    }
    data.symmetries.reserve(lattice.permutations.size());
    for(const auto& symmetry : lattice.permutations){
        PtmSymmetryPermutation ptmSymmetry;
        ptmSymmetry.transformation = symmetry.transformation;
        ptmSymmetry.inverseProduct = symmetry.inverseProduct;
        ptmSymmetry.permutation.assign(
            symmetry.permutation.begin(),
            symmetry.permutation.begin() + coord.numNeighbors
        );
        data.symmetries.push_back(std::move(ptmSymmetry));
    }

    return data;
}

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

    if(static_cast<StructureType>(normalizedStructureType(structureType)) == StructureType::BCC){
        double maxAbsComponent = 0.0;
        for(const Vector3& vector : latticeVectors){
            maxAbsComponent = std::max(maxAbsComponent, std::abs(vector.x()));
            maxAbsComponent = std::max(maxAbsComponent, std::abs(vector.y()));
            maxAbsComponent = std::max(maxAbsComponent, std::abs(vector.z()));
        }

        if(maxAbsComponent > EPSILON){
            const double scale = maxAbsComponent;
            for(Vector3& vector : latticeVectors){
                vector /= scale;
            }
        }
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

PtmCrystalInfoProvider::PtmCrystalInfoProvider(){
    initialize(static_cast<int>(StructureType::SC));
    initialize(static_cast<int>(StructureType::FCC));
    initialize(static_cast<int>(StructureType::HCP));
    initialize(static_cast<int>(StructureType::BCC));
    initialize(static_cast<int>(StructureType::ICO));
    initialize(static_cast<int>(StructureType::CUBIC_DIAMOND));
    initialize(static_cast<int>(StructureType::HEX_DIAMOND));
    initialize(static_cast<int>(StructureType::GRAPHENE));
}

int PtmCrystalInfoProvider::findClosestSymmetryPermutation(int structureType, const Matrix3& rotation) const{
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

int PtmCrystalInfoProvider::coordinationNumber(int structureType) const{
    return dataFor(structureType).coordinationNumber;
}

int PtmCrystalInfoProvider::commonNeighborIndex(int structureType, int neighborIndex, int commonNeighborSlot) const{
    return dataFor(structureType).commonNeighbors[static_cast<std::size_t>(neighborIndex)]
        [static_cast<std::size_t>(commonNeighborSlot)];
}

int PtmCrystalInfoProvider::symmetryPermutationCount(int structureType) const{
    return static_cast<int>(dataFor(structureType).symmetries.size());
}

int PtmCrystalInfoProvider::symmetryPermutationEntry(int structureType, int symmetryIndex, int neighborIndex) const{
    return dataFor(structureType).symmetries[static_cast<std::size_t>(symmetryIndex)]
        .permutation[static_cast<std::size_t>(neighborIndex)];
}

const Matrix3& PtmCrystalInfoProvider::symmetryTransformation(int structureType, int symmetryIndex) const{
    return dataFor(structureType).symmetries[static_cast<std::size_t>(symmetryIndex)].transformation;
}

int PtmCrystalInfoProvider::symmetryInverseProduct(int structureType, int symmetryIndex, int transformationIndex) const{
    return dataFor(structureType).symmetries[static_cast<std::size_t>(symmetryIndex)]
        .inverseProduct[static_cast<std::size_t>(transformationIndex)];
}

const Vector3& PtmCrystalInfoProvider::latticeVector(int structureType, int latticeVectorIndex) const{
    return dataFor(structureType).latticeVectors[static_cast<std::size_t>(latticeVectorIndex)];
}

void PtmCrystalInfoProvider::initialize(int structureType){
    const StructureType normalized = static_cast<StructureType>(normalizedStructureType(structureType));
    if(normalized == StructureType::CUBIC_DIAMOND || normalized == StructureType::HEX_DIAMOND){
        _data[static_cast<std::size_t>(normalizedStructureType(structureType))] =
            buildCnaOrderedDiamondCrystalData(structureType);
        return;
    }

    try{
        _data[static_cast<std::size_t>(normalizedStructureType(structureType))] =
            buildCrystalData(structureType);
    }catch(const std::exception&){
        _data[static_cast<std::size_t>(normalizedStructureType(structureType))] = PtmCrystalData{};
    }
}

const PtmCrystalData& PtmCrystalInfoProvider::dataFor(int structureType) const{
    return _data[static_cast<std::size_t>(normalizedStructureType(structureType))];
}

std::shared_ptr<const StructureAnalysisCrystalInfo> ptmCrystalInfoProvider(){
    static const auto provider = std::make_shared<PtmCrystalInfoProvider>();
    return provider;
}

} // namespace Volt::PtmStructureAnalysisDetail

namespace Volt {

void computeMaximumNeighborDistanceFromPTM(StructureAnalysis& analysis){
    StructureContext& context = analysis.context();
    const size_t N = context.atomCount();
    if(N == 0){
        context.maximumNeighborDistance = 0.0;
        return;
    }

    const auto* positions = context.positions->constDataPoint3();
    const auto& inverseMatrix = context.simCell.inverseMatrix();
    const auto& directMatrix = context.simCell.matrix();
    const int* counts = context.neighborCounts->constDataInt();
    const int* offsets = context.neighborOffsets->constDataInt();
    const int* indices = context.neighborIndices->constDataInt();

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

    context.maximumNeighborDistance = maxDistance;
}

void determineLocalStructuresWithPTM(
    StructureAnalysis& analysis,
    double rmsdCutoff
) {
    StructureContext& context = analysis.context();
    const size_t N = context.atomCount();
    if(!N){
        return;
    }
    analysis.setCrystalInfoProvider(PtmStructureAnalysisDetail::ptmCrystalInfoProvider());

    Volt::PTM ptm;
    if(!PtmStructureAnalysisDetail::setupPTM(context, ptm, N)){
        throw std::runtime_error("Error trying to initialize PTM.");
    }

    std::fill(context.neighborCounts->dataInt(),
              context.neighborCounts->dataInt() + context.neighborCounts->size(), 0);
    std::fill(context.structureTypes->dataInt(),
              context.structureTypes->dataInt() + context.structureTypes->size(), LATTICE_OTHER);
    std::fill(
        context.atomAllowedSymmetryMasks->dataInt64(),
        context.atomAllowedSymmetryMasks->dataInt64() + context.atomAllowedSymmetryMasks->size(),
        0
    );

    std::vector<uint64_t> cached(N, 0ull);
    std::vector<int> localCounts(N, 0);
    std::vector<std::array<int, MAX_NEIGHBORS>> canonicalDiamondNeighbors;
    std::vector<unsigned char> canonicalDiamondShellValid;

    tbb::parallel_for(tbb::blocked_range<size_t>(0, N), [&](const auto &range){
        PTM::Kernel kernel(ptm);

        for(size_t i = range.begin(); i < range.end(); ++i){
            kernel.cacheNeighbors(i, &cached[i]);
            StructureType type = kernel.identifyStructure(i, cached);
            if(type == StructureType::OTHER || kernel.rmsd() > rmsdCutoff){
                continue;
            }

            context.structureTypes->setInt(i, type);
            context.atomAllowedSymmetryMasks->setInt64(
                i,
                static_cast<std::int64_t>(
                    PtmStructureAnalysisDetail::fullSymmetryMask(
                        PtmStructureAnalysisDetail::ptmCrystalInfoProvider()->symmetryPermutationCount(type)
                    )
                )
            );

            const int neighborCount = kernel.numTemplateNeighbors();
            localCounts[i] = neighborCount;
            context.neighborCounts->setInt(i, neighborCount);
        }
    });

    if(
        context.inputCrystalType == LATTICE_CUBIC_DIAMOND ||
        context.inputCrystalType == LATTICE_HEX_DIAMOND
    ){
        PtmStructureAnalysisDetail::ensureCoordinationStructuresInitialized();
        canonicalDiamondNeighbors.resize(N);
        canonicalDiamondShellValid.assign(N, 0);
        for(auto& row : canonicalDiamondNeighbors){
            row.fill(-1);
        }

        auto cnaOrderedStructureTypes = std::make_shared<ParticleProperty>(N, DataType::Int, 1, 0, true);
        CoordinationStructures diamondCoordinationStructures(
            cnaOrderedStructureTypes.get(),
            context.inputCrystalType,
            true,
            context.simCell
        );
        NearestNeighborFinder diamondOrderingFinder(MAX_NEIGHBORS);
        if(!diamondOrderingFinder.prepare(context.positions, context.simCell, context.particleSelection)){
            throw std::runtime_error("Error trying to prepare PTM diamond ordering finder.");
        }

        tbb::parallel_for(tbb::blocked_range<size_t>(0, N), [&](const auto& range){
            for(size_t i = range.begin(); i < range.end(); ++i){
                const int structureType = context.structureTypes->getInt(i);
                if(
                    structureType != StructureType::CUBIC_DIAMOND &&
                    structureType != StructureType::HEX_DIAMOND
                ){
                    continue;
                }
                if(localCounts[i] != 16){
                    continue;
                }

                int cnaOrderedCount = 0;
                if(diamondCoordinationStructures.determineLocalStructure(
                    diamondOrderingFinder,
                    static_cast<int>(i),
                    &cnaOrderedCount,
                    canonicalDiamondNeighbors[i].data()
                ) > 0.0){
                    const int cnaOrderedType = cnaOrderedStructureTypes->getInt(i);
                    if(
                        cnaOrderedCount == 16 &&
                        (cnaOrderedType == StructureType::CUBIC_DIAMOND || cnaOrderedType == StructureType::HEX_DIAMOND)
                    ){
                        if(cnaOrderedType != structureType){
                            context.structureTypes->setInt(i, cnaOrderedType);
                        }
                        canonicalDiamondShellValid[i] = 1;
                        continue;
                    }
                }
            }
        });

        for(size_t i = 0; i < N; ++i){
            const int structureType = context.structureTypes->getInt(i);
            if(
                structureType != StructureType::CUBIC_DIAMOND &&
                structureType != StructureType::HEX_DIAMOND
            ){
                continue;
            }

            if(canonicalDiamondShellValid[i]){
                continue;
            }

            context.structureTypes->setInt(i, LATTICE_OTHER);
            context.neighborCounts->setInt(i, 0);
            localCounts[i] = 0;
        }
    }

    auto* offsets = context.neighborOffsets->dataInt();
    offsets[0] = 0;
    for(size_t i = 0; i < N; ++i){
        offsets[i + 1] = offsets[i] + localCounts[i];
    }

    const size_t totalNeighbors = static_cast<size_t>(offsets[N]);
    context.neighborIndices = std::make_shared<ParticleProperty>(totalNeighbors, DataType::Int, 1, 0, false);
    auto* indices = context.neighborIndices->dataInt();

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

            const int structureType = context.structureTypes->getInt(i);
            if(
                (structureType == StructureType::CUBIC_DIAMOND || structureType == StructureType::HEX_DIAMOND) &&
                count == 16 &&
                canonicalDiamondShellValid[i]
            ){
                for(int neighborSlot = 0; neighborSlot < count; ++neighborSlot){
                    indices[start + neighborSlot] = canonicalDiamondNeighbors[i][static_cast<std::size_t>(neighborSlot)];
                }
                continue;
            }
            
            if(structureType == StructureType::SC && count == static_cast<int>(PtmStructureAnalysisDetail::kSimpleCubicTemplateToCanonicalNeighborSlot.size())) {
                for(int templateSlot = 0; templateSlot < count; ++templateSlot){
                    const int canonicalSlot = PtmStructureAnalysisDetail::kSimpleCubicTemplateToCanonicalNeighborSlot[templateSlot];
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
        if(context.neighborCounts->getInt(i) == 0){
            context.neighborCounts->setInt(i, localCounts[i]);
        }
    }
}

}
