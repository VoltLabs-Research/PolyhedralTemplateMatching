#include <volt/analysis/internal/ptm_crystal_info_provider.h>

#include <volt/analysis/crystal_symmetry_utils.h>
#include <volt/analysis/crystal_topology_library.h>
#include <volt/analysis/ptm.h>
#include <volt/structures/crystal_topology_registry.h>

#include <ptm_constants.h>
#include <ptm_initialize_data.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Volt::PtmStructureAnalysisDetail {
namespace {

const Vector3& zeroVector(){
    static const Vector3 vector = Vector3::Zero();
    return vector;
}

int normalizedStructureType(int structureType){
    if(const auto* topology = crystalTopologyByStructureType(structureType)){
        return topology->structureType;
    }
    return structureType;
}

bool shouldPreserveNativePtmReference(int structureType){
    // FCC is consumed downstream with PTM's historical neighbor/symmetry gauge.
    // Re-adapting it through the shared topology registry changes the reconstructed
    // dump and cluster transitions enough to alter DXA output.
    switch(static_cast<StructureType>(normalizedStructureType(structureType))){
        case StructureType::FCC:
            return true;
        default:
            return false;
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

std::vector<Vector3> templateVectorsForReference(const ptm::refdata_t& ref){
    const int templateIndex = canonicalTemplateIndex(ref);
    const double (*points)[3] = ref.points[templateIndex];

    std::vector<Vector3> vectors(static_cast<std::size_t>(ref.num_nbrs), Vector3::Zero());
    for(int templateSlot = 0; templateSlot < ref.num_nbrs; ++templateSlot){
        const double* point = points[templateSlot + 1];
        vectors[static_cast<std::size_t>(templateSlot)] = Vector3(point[0], point[1], point[2]);
    }
    return vectors;
}

Vector3 normalizedDirection(const Vector3& vector){
    const double length = vector.length();
    return length > EPSILON ? (vector / length) : Vector3::Zero();
}

std::vector<int> identityMapping(int count){
    std::vector<int> mapping(static_cast<std::size_t>(count), 0);
    std::iota(mapping.begin(), mapping.end(), 0);
    return mapping;
}

std::vector<int> buildTemplateToCanonicalMapping(
    const std::vector<Vector3>& templateVectors,
    const std::vector<Vector3>& canonicalVectors
){
    if(templateVectors.size() != canonicalVectors.size()){
        throw std::runtime_error("PTM template/canonical vector count mismatch.");
    }

    struct CandidatePair{
        int templateSlot = -1;
        int canonicalSlot = -1;
        double error = std::numeric_limits<double>::max();
    };

    std::vector<CandidatePair> candidates;
    candidates.reserve(templateVectors.size() * canonicalVectors.size());
    for(int templateSlot = 0; templateSlot < static_cast<int>(templateVectors.size()); ++templateSlot){
        const Vector3 templateDirection = normalizedDirection(templateVectors[static_cast<std::size_t>(templateSlot)]);
        for(int canonicalSlot = 0; canonicalSlot < static_cast<int>(canonicalVectors.size()); ++canonicalSlot){
            const Vector3 canonicalDirection = normalizedDirection(canonicalVectors[static_cast<std::size_t>(canonicalSlot)]);
            candidates.push_back({
                templateSlot,
                canonicalSlot,
                (templateDirection - canonicalDirection).squaredLength()
            });
        }
    }

    std::sort(candidates.begin(), candidates.end(), [](const CandidatePair& left, const CandidatePair& right){
        if(std::abs(left.error - right.error) > 1e-12){
            return left.error < right.error;
        }
        if(left.templateSlot != right.templateSlot){
            return left.templateSlot < right.templateSlot;
        }
        return left.canonicalSlot < right.canonicalSlot;
    });

    std::vector<int> mapping(templateVectors.size(), -1);
    std::vector<unsigned char> templateAssigned(templateVectors.size(), 0);
    std::vector<unsigned char> canonicalAssigned(canonicalVectors.size(), 0);
    for(const CandidatePair& candidate : candidates){
        if(candidate.error > 1e-4){
            break;
        }
        if(templateAssigned[static_cast<std::size_t>(candidate.templateSlot)] ||
           canonicalAssigned[static_cast<std::size_t>(candidate.canonicalSlot)]){
            continue;
        }
        mapping[static_cast<std::size_t>(candidate.templateSlot)] = candidate.canonicalSlot;
        templateAssigned[static_cast<std::size_t>(candidate.templateSlot)] = 1;
        canonicalAssigned[static_cast<std::size_t>(candidate.canonicalSlot)] = 1;
    }

    if(std::find(mapping.begin(), mapping.end(), -1) != mapping.end()){
        throw std::runtime_error("Unable to align PTM template order with canonical topology.");
    }

    return mapping;
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

std::vector<std::array<int, 2>> buildCommonNeighbors(const std::vector<Vector3>& latticeVectors){
    const int numNeighbors = static_cast<int>(latticeVectors.size());
    std::vector<std::array<int, 2>> commonNeighbors(
        static_cast<std::size_t>(numNeighbors),
        std::array<int, 2>{-1, -1}
    );

    for(int neighborIndex = 0; neighborIndex < numNeighbors; ++neighborIndex){
        Matrix3 basis = Matrix3::Zero();
        basis.column(0) = latticeVectors[static_cast<std::size_t>(neighborIndex)];

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
            continue;
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

    next_neighbor:
        continue;
    }

    return commonNeighbors;
}

std::vector<PtmSymmetryPermutation> buildSymmetriesFromReference(
    const ptm::refdata_t& ref,
    const std::vector<Vector3>& canonicalVectors,
    const std::vector<int>& templateToCanonical
){
    std::vector<PtmSymmetryPermutation> symmetries;
    if(canonicalVectors.empty()){
        return symmetries;
    }

    const auto basisIndices = findNonCoplanarIndices(canonicalVectors);
    Matrix3 basis = Matrix3::Zero();
    basis.column(0) = canonicalVectors[static_cast<std::size_t>(basisIndices[0])];
    basis.column(1) = canonicalVectors[static_cast<std::size_t>(basisIndices[1])];
    basis.column(2) = canonicalVectors[static_cast<std::size_t>(basisIndices[2])];
    Matrix3 basisInverse;
    if(!basis.inverse(basisInverse)){
        throw std::runtime_error("Unable to invert PTM canonical basis.");
    }

    const int count = symmetryMappingCount(ref);
    symmetries.reserve(static_cast<std::size_t>(count));
    for(int mappingIndex = 0; mappingIndex < count; ++mappingIndex){
        const int8_t* mapping = symmetryMappings(ref)[mappingIndex];

        PtmSymmetryPermutation symmetry;
        symmetry.permutation.assign(canonicalVectors.size(), -1);
        for(int templateSlot = 0; templateSlot < static_cast<int>(templateToCanonical.size()); ++templateSlot){
            const int fromCanonical = templateToCanonical[static_cast<std::size_t>(templateSlot)];
            const int mappedTemplateSlot = mapping[templateSlot + 1] - 1;
            if(fromCanonical < 0 || mappedTemplateSlot < 0 ||
               mappedTemplateSlot >= static_cast<int>(templateToCanonical.size())){
                continue;
            }
            const int toCanonical = templateToCanonical[static_cast<std::size_t>(mappedTemplateSlot)];
            if(toCanonical < 0){
                continue;
            }
            symmetry.permutation[static_cast<std::size_t>(fromCanonical)] = toCanonical;
        }

        if(std::find(symmetry.permutation.begin(), symmetry.permutation.end(), -1) != symmetry.permutation.end()){
            continue;
        }

        symmetry.transformation = Matrix3::Zero();
        symmetry.transformation.column(0) = canonicalVectors[
            static_cast<std::size_t>(symmetry.permutation[static_cast<std::size_t>(basisIndices[0])])
        ];
        symmetry.transformation.column(1) = canonicalVectors[
            static_cast<std::size_t>(symmetry.permutation[static_cast<std::size_t>(basisIndices[1])])
        ];
        symmetry.transformation.column(2) = canonicalVectors[
            static_cast<std::size_t>(symmetry.permutation[static_cast<std::size_t>(basisIndices[2])])
        ];
        symmetry.transformation = symmetry.transformation * basisInverse;

        bool duplicate = false;
        for(const auto& existing : symmetries){
            if(existing.transformation.equals(symmetry.transformation)){
                duplicate = true;
                break;
            }
        }
        if(!duplicate){
            symmetries.push_back(std::move(symmetry));
        }
    }

    AnalysisSymmetryUtils::calculateSymmetryProducts(symmetries);
    return symmetries;
}

PtmCrystalData buildSharedCrystalData(int structureType){
    PtmCrystalData data;
    const SharedCrystalTopology* topology = sharedCrystalTopology(structureType);
    if(!topology){
        return data;
    }

    if(const ptm::refdata_t* ref = ptmReferenceForStructureType(structureType)){
        const std::vector<Vector3> templateVectors = templateVectorsForReference(*ref);
        if(static_cast<int>(templateVectors.size()) == topology->coordinationNumber){
            AdaptedCrystalTopology adapted;
            if(adaptSharedCrystalTopology(*topology, templateVectors, adapted)){
                data.coordinationNumber = adapted.coordinationNumber;
                data.latticeVectors = adapted.latticeVectors;
                data.commonNeighbors.resize(
                    static_cast<std::size_t>(adapted.coordinationNumber),
                    std::array<int, 2>{-1, -1}
                );
                for(int neighborIndex = 0; neighborIndex < adapted.coordinationNumber; ++neighborIndex){
                    data.commonNeighbors[static_cast<std::size_t>(neighborIndex)] =
                        adapted.commonNeighbors[static_cast<std::size_t>(neighborIndex)];
                }
                data.symmetries.reserve(adapted.symmetries.size());
                for(const auto& symmetry : adapted.symmetries){
                    PtmSymmetryPermutation converted;
                    converted.transformation = symmetry.transformation;
                    converted.inverseProduct = symmetry.inverseProduct;
                    converted.permutation.assign(
                        symmetry.permutation.begin(),
                        symmetry.permutation.begin() + adapted.coordinationNumber
                    );
                    data.symmetries.push_back(std::move(converted));
                }
                data.templateToCanonicalNeighborSlot = identityMapping(adapted.coordinationNumber);
                return data;
            }
        }
    }

    data.coordinationNumber = topology->coordinationNumber;
    data.latticeVectors.assign(
        topology->latticeVectors.begin(),
        topology->latticeVectors.begin() + topology->coordinationNumber
    );
    data.commonNeighbors.resize(
        static_cast<std::size_t>(topology->coordinationNumber),
        std::array<int, 2>{-1, -1}
    );
    for(int neighborIndex = 0; neighborIndex < topology->coordinationNumber; ++neighborIndex){
        data.commonNeighbors[static_cast<std::size_t>(neighborIndex)] =
            topology->commonNeighbors[static_cast<std::size_t>(neighborIndex)];
    }
    data.symmetries.reserve(topology->symmetries.size());
    for(const auto& symmetry : topology->symmetries){
        PtmSymmetryPermutation converted;
        converted.transformation = symmetry.transformation;
        converted.inverseProduct = symmetry.inverseProduct;
        converted.permutation.assign(
            symmetry.permutation.begin(),
            symmetry.permutation.begin() + topology->coordinationNumber
        );
        data.symmetries.push_back(std::move(converted));
    }

    if(const ptm::refdata_t* ref = ptmReferenceForStructureType(structureType)){
        data.templateToCanonicalNeighborSlot = buildTemplateToCanonicalMapping(
            templateVectorsForReference(*ref),
            data.latticeVectors
        );
    }else{
        data.templateToCanonicalNeighborSlot = identityMapping(topology->coordinationNumber);
    }

    return data;
}

PtmCrystalData buildReferenceCrystalData(int structureType){
    PtmCrystalData data;
    const ptm::refdata_t* ref = ptmReferenceForStructureType(structureType);
    if(!ref){
        return data;
    }

    data.coordinationNumber = ref->num_nbrs;
    data.latticeVectors = templateVectorsForReference(*ref);
    data.commonNeighbors = buildCommonNeighbors(data.latticeVectors);
    data.templateToCanonicalNeighborSlot = identityMapping(data.coordinationNumber);
    data.symmetries = buildSymmetriesFromReference(*ref, data.latticeVectors, data.templateToCanonicalNeighborSlot);
    return data;
}

const PtmCrystalData& emptyCrystalData(){
    static const PtmCrystalData empty;
    return empty;
}

std::shared_ptr<const PtmCrystalInfoProvider> ptmCrystalInfoProviderImpl(){
    static const auto provider = std::make_shared<PtmCrystalInfoProvider>();
    return provider;
}

} // namespace

PtmCrystalInfoProvider::PtmCrystalInfoProvider(){
    initialize(static_cast<int>(StructureType::ICO));
    initialize(static_cast<int>(StructureType::GRAPHENE));
}

int PtmCrystalInfoProvider::findClosestSymmetryPermutation(int structureType, const Matrix3& rotation) const{
    return AnalysisSymmetryUtils::findClosestSymmetryPermutation(dataFor(structureType).symmetries, rotation);
}

int PtmCrystalInfoProvider::coordinationNumber(int structureType) const{
    return dataFor(structureType).coordinationNumber;
}

int PtmCrystalInfoProvider::commonNeighborIndex(int structureType, int neighborIndex, int commonNeighborSlot) const{
    const auto& commonNeighbors = dataFor(structureType).commonNeighbors;
    if(neighborIndex < 0 || neighborIndex >= static_cast<int>(commonNeighbors.size()) ||
       commonNeighborSlot < 0 || commonNeighborSlot > 1){
        return -1;
    }
    return commonNeighbors[static_cast<std::size_t>(neighborIndex)][static_cast<std::size_t>(commonNeighborSlot)];
}

int PtmCrystalInfoProvider::symmetryPermutationCount(int structureType) const{
    return static_cast<int>(dataFor(structureType).symmetries.size());
}

int PtmCrystalInfoProvider::symmetryPermutationEntry(int structureType, int symmetryIndex, int neighborIndex) const{
    const auto& data = dataFor(structureType);
    if(symmetryIndex < 0 || symmetryIndex >= static_cast<int>(data.symmetries.size())){
        return neighborIndex;
    }
    const auto& permutation = data.symmetries[static_cast<std::size_t>(symmetryIndex)].permutation;
    if(neighborIndex < 0 || neighborIndex >= static_cast<int>(permutation.size())){
        return neighborIndex;
    }
    return permutation[static_cast<std::size_t>(neighborIndex)];
}

const Matrix3& PtmCrystalInfoProvider::symmetryTransformation(int structureType, int symmetryIndex) const{
    const auto& data = dataFor(structureType);
    if(symmetryIndex < 0 || symmetryIndex >= static_cast<int>(data.symmetries.size())){
        static const Matrix3 identity = Matrix3::Identity();
        return identity;
    }
    return data.symmetries[static_cast<std::size_t>(symmetryIndex)].transformation;
}

int PtmCrystalInfoProvider::symmetryInverseProduct(int structureType, int symmetryIndex, int transformationIndex) const{
    const auto& data = dataFor(structureType);
    if(symmetryIndex < 0 || symmetryIndex >= static_cast<int>(data.symmetries.size())){
        return 0;
    }
    const auto& inverseProduct = data.symmetries[static_cast<std::size_t>(symmetryIndex)].inverseProduct;
    if(transformationIndex < 0 || transformationIndex >= static_cast<int>(inverseProduct.size())){
        return 0;
    }
    return inverseProduct[static_cast<std::size_t>(transformationIndex)];
}

const Vector3& PtmCrystalInfoProvider::latticeVector(int structureType, int latticeVectorIndex) const{
    const auto& data = dataFor(structureType);
    if(latticeVectorIndex < 0 || latticeVectorIndex >= static_cast<int>(data.latticeVectors.size())){
        return zeroVector();
    }
    return data.latticeVectors[static_cast<std::size_t>(latticeVectorIndex)];
}

int PtmCrystalInfoProvider::templateToCanonicalNeighborSlot(int structureType, int templateSlot) const{
    const auto& mapping = dataFor(structureType).templateToCanonicalNeighborSlot;
    if(templateSlot < 0 || templateSlot >= static_cast<int>(mapping.size())){
        return templateSlot;
    }
    return mapping[static_cast<std::size_t>(templateSlot)];
}

void PtmCrystalInfoProvider::initialize(int structureType) const{
    const int normalizedType = normalizedStructureType(structureType);
    if(_data.find(normalizedType) != _data.end()){
        return;
    }

    if(shouldPreserveNativePtmReference(normalizedType)){
        try{
            _data[normalizedType] = buildReferenceCrystalData(normalizedType);
            return;
        }catch(const std::exception&){
            _data[normalizedType] = PtmCrystalData{};
            return;
        }
    }

    if(sharedCrystalTopology(normalizedType)){
        _data[normalizedType] = buildSharedCrystalData(normalizedType);
        return;
    }

    try{
        _data[normalizedType] = buildReferenceCrystalData(normalizedType);
    }catch(const std::exception&){
        _data[normalizedType] = PtmCrystalData{};
    }
}

const PtmCrystalData& PtmCrystalInfoProvider::dataFor(int structureType) const{
    const int normalizedType = normalizedStructureType(structureType);
    if(_data.find(normalizedType) == _data.end()){
        initialize(normalizedType);
    }
    const auto it = _data.find(normalizedType);
    return it != _data.end() ? it->second : emptyCrystalData();
}

std::shared_ptr<const StructureAnalysisCrystalInfo> ptmCrystalInfoProvider(){
    return ptmCrystalInfoProviderImpl();
}

int ptmTemplateToCanonicalNeighborSlot(int structureType, int templateSlot){
    return ptmCrystalInfoProviderImpl()->templateToCanonicalNeighborSlot(structureType, templateSlot);
}

} 
