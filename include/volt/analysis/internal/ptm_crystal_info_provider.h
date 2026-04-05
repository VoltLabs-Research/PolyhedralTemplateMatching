#pragma once

#include <volt/analysis/structure_analysis.h>

#include <array>
#include <cstddef>
#include <memory>
#include <vector>

namespace Volt::PtmStructureAnalysisDetail{

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
