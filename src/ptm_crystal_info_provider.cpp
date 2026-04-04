#include "ptm_crystal_info_provider.h"
#include "ptm_structure_analysis_detail.h"

#include <volt/coordination_structures.h>
#include <volt/polyhedral_template_matching.h>

#include <ptm_constants.h>
#include <ptm_initialize_data.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace Volt::PtmStructureAnalysisDetail {

inline constexpr std::array<int, 6> kSimpleCubicCanonicalToTemplateSlot = {5, 4, 3, 2, 1, 0};

PtmCrystalData buildCanonicalDiamondCrystalData(int structureType){
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
            latticeVectors[static_cast<std::size_t>(canonicalSlot)] = Vector3(point[0], point[1], point[2]);
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

std::vector<int> buildSymmetryPermutation(int structureType, const ptm::refdata_t& ref, int mappingIndex){
    const int numNeighbors = ref.num_nbrs;
    const int8_t* mapping = symmetryMappings(ref)[mappingIndex];
    std::vector<int> permutation(static_cast<std::size_t>(numNeighbors), 0);

    if(static_cast<StructureType>(normalizedStructureType(structureType)) == StructureType::SC){
        for(int canonicalSlot = 0; canonicalSlot < numNeighbors; ++canonicalSlot){
            const int templateSlot = kSimpleCubicCanonicalToTemplateSlot[static_cast<std::size_t>(canonicalSlot)];
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
        const Matrix3& symmetry = symmetries[static_cast<std::size_t>(symmetryIndex)].transformation;
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
    const int normalizedType = normalizedStructureType(structureType);
    const auto normalizedIndex = static_cast<std::size_t>(normalizedType);
    const StructureType normalized = static_cast<StructureType>(normalizedType);
    if(normalized == StructureType::CUBIC_DIAMOND || normalized == StructureType::HEX_DIAMOND){
        _data[normalizedIndex] = buildCanonicalDiamondCrystalData(structureType);
        return;
    }

    try{
        _data[normalizedIndex] = buildCrystalData(structureType);
    }catch(const std::exception&){
        _data[normalizedIndex] = PtmCrystalData{};
    }
}

const PtmCrystalData& PtmCrystalInfoProvider::dataFor(int structureType) const{
    return _data[static_cast<std::size_t>(normalizedStructureType(structureType))];
}

std::shared_ptr<const StructureAnalysisCrystalInfo> ptmCrystalInfoProvider(){
    static const auto provider = std::make_shared<PtmCrystalInfoProvider>();
    return provider;
}

} 
