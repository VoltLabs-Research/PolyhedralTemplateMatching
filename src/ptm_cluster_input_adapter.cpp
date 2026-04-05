#include <volt/ptm_cluster_input_adapter.h>

#include "ptm_structure_analysis_detail.h"

#include <cstdint>

namespace Volt{

void PTMClusterInputAdapter::prepare(StructureAnalysis& analysis, AnalysisContext& context){
    // PTM should provide cluster inputs, not replace the generic cluster-growth logic.
    analysis.setClusterRuleProvider(nullptr);

    if(!context.atomAllowedSymmetryMasks){
        context.atomAllowedSymmetryMasks = std::make_shared<ParticleProperty>(
            context.atomCount(),
            DataType::Int64,
            1,
            0,
            true
        );
    }

    std::fill(
        context.atomSymmetryPermutations->dataInt(),
        context.atomSymmetryPermutations->dataInt() + context.atomSymmetryPermutations->size(),
        -1
    );
    for(size_t atomIndex = 0; atomIndex < context.atomCount(); ++atomIndex){
        const int structureType = context.structureTypes->getInt(atomIndex);
        if(structureType == LATTICE_OTHER || analysis.numberOfNeighbors(static_cast<int>(atomIndex)) == 0){
            continue;
        }

        if(context.atomAllowedSymmetryMasks->getInt64(atomIndex) != 0){
            continue;
        }

        context.atomAllowedSymmetryMasks->setInt64(
            atomIndex,
            static_cast<std::int64_t>(
                PtmStructureAnalysisDetail::fullSymmetryMask(analysis.symmetryPermutationCount(structureType))
            )
        );
    }
}

}
