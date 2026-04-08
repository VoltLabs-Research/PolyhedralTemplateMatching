#pragma once

#include <volt/analysis/ptm_local_atom_state.h>
#include <volt/analysis/structure_analysis.h>

#include <memory>
#include <vector>

namespace Volt{

void computeMaximumNeighborDistanceFromPTM(StructureAnalysis& analysis);
void determineLocalStructuresWithPTM(
    StructureAnalysis& analysis,
    double rmsdCutoff,
    std::shared_ptr<std::vector<PtmLocalAtomState>> atomStates = nullptr
);

}
