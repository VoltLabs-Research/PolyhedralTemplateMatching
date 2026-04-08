#pragma once

#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/ptm.h>
#include <volt/analysis/internal/ptm_crystal_info_provider.h>

#include <cstddef>

namespace Volt::PtmStructureAnalysisDetail{

bool setupPTM(StructureContext& context, Volt::PTM& ptm, std::size_t particleCount);

}
