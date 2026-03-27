#pragma once

#include <volt/analysis/analysis_context.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/cluster_builder.h>

#include <deque>

namespace Volt{

class PTMClusterBuilder final : public ClusterBuilder{
public:
    PTMClusterBuilder(
        StructureAnalysis& sa,
        AnalysisContext& context,
        std::shared_ptr<const ParticleProperty> ptmRmsd,
        std::shared_ptr<const ParticleProperty> orientations
    );

    bool areOrientationsCompatible(int atom1, int atom2, int structureType);

    void build(bool dissolveSmallClusters = false);

    void initOrientation(Cluster* cluster, size_t seedAtomIndex);
    void grow(
        Cluster* cluster, 
        std::deque<int>& atomsToVisit, 
        int structureType
    );

    Quaternion getAtomOrientation(int atom) const;

private:
    std::shared_ptr<const ParticleProperty> _ptmRmsd;
    std::shared_ptr<const ParticleProperty> _orientations;
    void reorientAtomsToAlign();
    void connectClusters();
    void formSuperClusters();
    void processDefectCluster(Cluster* defectCluster);
    bool calculateMisorientation(
        int atomIndex,
        int neighbor,
        int neighborIndex,
        Matrix3& outTransition
    ) const;
};

}
