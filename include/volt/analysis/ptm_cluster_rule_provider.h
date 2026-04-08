#pragma once

#include <volt/analysis/cluster_rule_provider.h>
#include <volt/analysis/ptm_local_atom_state.h>

#include <memory>
#include <vector>

namespace Volt{

class PtmClusterRuleProvider final : public ClusterRuleProvider{
public:
    explicit PtmClusterRuleProvider(
        std::shared_ptr<const std::vector<PtmLocalAtomState>> atomStates
    );

    void initializeClusterSeed(
        const StructureAnalysis& analysis,
        const AnalysisContext& context,
        Cluster& cluster,
        int seedAtomIndex,
        int structureType
    ) const override;

    ClusterRuleDecision tryAssignNeighbor(
        const StructureAnalysis& analysis,
        const AnalysisContext& context,
        const Cluster& cluster,
        int currentAtomIndex,
        int neighborAtomIndex,
        int neighborIndex,
        int structureType,
        int& outNeighborSymmetry
    ) const override;

    ClusterRuleDecision tryCalculateTransition(
        const StructureAnalysis& analysis,
        const AnalysisContext& context,
        int atomIndex,
        int neighborAtomIndex,
        int neighborIndex,
        Matrix3& outTransition
    ) const override;

private:
    const PtmLocalAtomState* stateFor(int atomIndex) const;
    static Matrix3 quaternionToMatrix(const Quaternion& orientation);

    std::shared_ptr<const std::vector<PtmLocalAtomState>> _atomStates;
};

}
