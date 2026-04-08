#include <volt/analysis/ptm_cluster_rule_provider.h>

#include <volt/analysis/structure_analysis.h>
#include <volt/structures/cluster.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace Volt{

namespace{

constexpr double kScOrientationThresholdRadians = 8.0 * PI / 180.0;
constexpr double kNeighborVectorTolerance = 1e-4;

bool symmetryAllowed(std::uint64_t mask, int symmetryIndex){
    return mask == 0 || (
        symmetryIndex >= 0 &&
        symmetryIndex < 63 &&
        (mask & (std::uint64_t{1} << symmetryIndex)) != 0
    );
}

double matrixTrace(const Matrix3& matrix){
    return matrix(0, 0) + matrix(1, 1) + matrix(2, 2);
}

bool orthonormalizeMatrix(const Matrix3& input, Matrix3& output){
    Vector3 c0 = input.column(0);
    Vector3 c1 = input.column(1);
    Vector3 c2 = input.column(2);

    const double l0 = c0.length();
    if(l0 <= EPSILON){
        return false;
    }
    c0 /= l0;

    c1 -= c0 * c0.dot(c1);
    const double l1 = c1.length();
    if(l1 <= EPSILON){
        return false;
    }
    c1 /= l1;

    c2 -= c0 * c0.dot(c2);
    c2 -= c1 * c1.dot(c2);
    const double l2 = c2.length();
    if(l2 <= EPSILON){
        return false;
    }
    c2 /= l2;

    output = Matrix3(c0, c1, c2);
    if(output.determinant() < 0.0){
        output.column(2) = -output.column(2);
    }
    return true;
}

}

PtmClusterRuleProvider::PtmClusterRuleProvider(
    std::shared_ptr<const std::vector<PtmLocalAtomState>> atomStates
) :
    _atomStates(std::move(atomStates)){}

void PtmClusterRuleProvider::initializeClusterSeed(
    const StructureAnalysis&,
    const AnalysisContext&,
    Cluster& cluster,
    int,
    int
) const{
    cluster.symmetryTransformation = 0;
}

const PtmLocalAtomState* PtmClusterRuleProvider::stateFor(int atomIndex) const{
    if(!_atomStates || atomIndex < 0 || atomIndex >= static_cast<int>(_atomStates->size())){
        return nullptr;
    }
    const auto& state = (*_atomStates)[static_cast<std::size_t>(atomIndex)];
    return state.valid ? &state : nullptr;
}

Matrix3 PtmClusterRuleProvider::quaternionToMatrix(const Quaternion& orientation){
    const Quaternion normalized = orientation.normalized();
    return Matrix3(
        normalized * Vector3(1.0, 0.0, 0.0),
        normalized * Vector3(0.0, 1.0, 0.0),
        normalized * Vector3(0.0, 0.0, 1.0)
    );
}

ClusterRuleDecision PtmClusterRuleProvider::tryAssignNeighbor(
    const StructureAnalysis& analysis,
    const AnalysisContext& context,
    const Cluster&,
    int currentAtomIndex,
    int neighborAtomIndex,
    int neighborIndex,
    int structureType,
    int& outNeighborSymmetry
) const{
    if(structureType != StructureType::SC){
        return ClusterRuleDecision::Unhandled;
    }

    const auto* currentState = stateFor(currentAtomIndex);
    const auto* neighborState = stateFor(neighborAtomIndex);
    if(!currentState || !neighborState){
        return ClusterRuleDecision::Unhandled;
    }
    if(context.structureTypes->getInt(neighborAtomIndex) != structureType){
        return ClusterRuleDecision::Rejected;
    }

    const int reverseSlot = analysis.findNeighbor(neighborAtomIndex, currentAtomIndex);
    if(reverseSlot < 0){
        return ClusterRuleDecision::Rejected;
    }

    const int currentSymmetry = std::max(0, context.atomSymmetryPermutations->getInt(currentAtomIndex));
    const int currentBondSlot = analysis.symmetryPermutationEntry(structureType, currentSymmetry, neighborIndex);
    if(currentBondSlot < 0){
        return ClusterRuleDecision::Rejected;
    }

    const Matrix3 rotationMatrix = quaternionToMatrix(
        currentState->orientation.inverse() * neighborState->orientation
    );
    const double minTrace = 1.0 + 2.0 * std::cos(kScOrientationThresholdRadians);
    const std::uint64_t allowedMask = static_cast<std::uint64_t>(
        context.atomAllowedSymmetryMasks->getInt64(static_cast<std::size_t>(neighborAtomIndex))
    );
    const int symmetryCount = analysis.symmetryPermutationCount(structureType);
    const Vector3 currentBond = analysis.latticeVector(structureType, currentBondSlot);

    double bestTrace = -std::numeric_limits<double>::infinity();
    int bestNeighborSymmetry = -1;

    for(int neighborSymmetry = 0; neighborSymmetry < symmetryCount; ++neighborSymmetry){
        if(!symmetryAllowed(allowedMask, neighborSymmetry)){
            continue;
        }

        const int neighborBondSlot = analysis.symmetryPermutationEntry(structureType, neighborSymmetry, reverseSlot);
        if(neighborBondSlot < 0){
            continue;
        }
        const Vector3 neighborBond = analysis.latticeVector(structureType, neighborBondSlot);
        if(!(currentBond + neighborBond).isZero(kNeighborVectorTolerance)){
            continue;
        }

        const int relativeSymmetry = analysis.symmetryInverseProduct(
            structureType,
            neighborSymmetry,
            currentSymmetry
        );
        const Matrix3 product = Matrix3(
            rotationMatrix * analysis.symmetryTransformation(structureType, relativeSymmetry).transposed()
        );
        const double trace = matrixTrace(product);
        if(trace > minTrace && trace > bestTrace){
            bestTrace = trace;
            bestNeighborSymmetry = neighborSymmetry;
        }
    }

    if(bestNeighborSymmetry < 0){
        return ClusterRuleDecision::Rejected;
    }

    outNeighborSymmetry = bestNeighborSymmetry;
    return ClusterRuleDecision::Accepted;
}

ClusterRuleDecision PtmClusterRuleProvider::tryCalculateTransition(
    const StructureAnalysis& analysis,
    const AnalysisContext& context,
    int atomIndex,
    int neighborAtomIndex,
    int,
    Matrix3& outTransition
) const{
    const int structureType = context.structureTypes->getInt(atomIndex);
    if(structureType != StructureType::SC || context.structureTypes->getInt(neighborAtomIndex) != structureType){
        return ClusterRuleDecision::Unhandled;
    }

    const auto* atomState = stateFor(atomIndex);
    const auto* neighborState = stateFor(neighborAtomIndex);
    if(!atomState || !neighborState){
        return ClusterRuleDecision::Unhandled;
    }

    const int atomSymmetry = context.atomSymmetryPermutations->getInt(atomIndex);
    const int neighborSymmetry = context.atomSymmetryPermutations->getInt(neighborAtomIndex);
    if(atomSymmetry < 0 || neighborSymmetry < 0){
        return ClusterRuleDecision::Rejected;
    }

    const Matrix3 atomOrientation = Matrix3(
        quaternionToMatrix(atomState->orientation) *
        analysis.symmetryTransformation(structureType, atomSymmetry).transposed()
    );
    const Matrix3 neighborOrientation = Matrix3(
        quaternionToMatrix(neighborState->orientation) *
        analysis.symmetryTransformation(structureType, neighborSymmetry).transposed()
    );

    Matrix3 rawTransition = Matrix3(neighborOrientation.transposed() * atomOrientation);
    if(!orthonormalizeMatrix(rawTransition, outTransition)){
        return ClusterRuleDecision::Rejected;
    }
    return ClusterRuleDecision::Accepted;
}

}
