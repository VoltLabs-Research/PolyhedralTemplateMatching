#pragma once

#include <volt/analysis/shared_crystal_topology.h>
#include <volt/polyhedral_template_matching.h>
#include <volt/ptm_structure_analysis.h>

#include <ptm_constants.h>
#include <ptm_initialize_data.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace Volt::PtmStructureAnalysisDetail{

inline constexpr std::array<int, 6> kSimpleCubicTemplateToCanonicalNeighborSlot = {5, 4, 3, 2, 1, 0};
inline constexpr std::array<int, 6> kSimpleCubicCanonicalToTemplateSlot = {5, 4, 3, 2, 1, 0};

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

bool setupPTM(StructureContext& context, Volt::PTM& ptm, std::size_t particleCount);
std::uint64_t fullSymmetryMask(int symmetryCount);
void ensureCoordinationStructuresInitialized();
PtmCrystalData toPtmCrystalData(const AdaptedCrystalTopology& adapted);
PtmCrystalData buildSharedAdaptedPtmCrystalData(int structureType);
PtmCrystalData buildCnaOrderedDiamondCrystalData(int structureType);
int normalizedStructureType(int structureType);
const ptm::refdata_t* ptmReferenceForStructureType(int structureType);
int canonicalTemplateIndex(const ptm::refdata_t& ref);
const int8_t (*symmetryMappings(const ptm::refdata_t& ref))[PTM_MAX_POINTS];
int symmetryMappingCount(const ptm::refdata_t& ref);
std::vector<Vector3> buildCanonicalLatticeVectors(int structureType, const ptm::refdata_t& ref);
std::vector<int> buildSymmetryPermutation(int structureType, const ptm::refdata_t& ref, int mappingIndex);
std::array<int, 3> findNonCoplanarIndices(const std::vector<Vector3>& latticeVectors);
std::vector<std::array<int, 2>> buildCommonNeighbors(
    int structureType,
    const std::vector<Vector3>& latticeVectors
);
PtmCrystalData buildCrystalData(int structureType);

class PtmCrystalInfoProvider final : public StructureAnalysisCrystalInfo{
public:
    PtmCrystalInfoProvider();

    int findClosestSymmetryPermutation(int structureType, const Matrix3& rotation) const override;
    int coordinationNumber(int structureType) const override;
    int commonNeighborIndex(int structureType, int neighborIndex, int commonNeighborSlot) const override;
    int symmetryPermutationCount(int structureType) const override;
    int symmetryPermutationEntry(int structureType, int symmetryIndex, int neighborIndex) const override;
    const Matrix3& symmetryTransformation(int structureType, int symmetryIndex) const override;
    int symmetryInverseProduct(int structureType, int symmetryIndex, int transformationIndex) const override;
    const Vector3& latticeVector(int structureType, int latticeVectorIndex) const override;

private:
    void initialize(int structureType);
    const PtmCrystalData& dataFor(int structureType) const;

    std::array<PtmCrystalData, static_cast<std::size_t>(StructureType::NUM_STRUCTURE_TYPES)> _data{};
};

std::shared_ptr<const StructureAnalysisCrystalInfo> ptmCrystalInfoProvider();

}
