#pragma once

#include <volt/analysis/structure_analysis.h>

namespace Volt{

void computeMaximumNeighborDistanceFromPTM(StructureAnalysis& analysis);
void determineLocalStructuresWithPTM(StructureAnalysis& analysis, double rmsdCutoff);

}
