#pragma once

#include <volt/core/volt.h>
#include <volt/analysis/analysis_context.h>
#include <volt/analysis/polyhedral_template_matching.h>
#include <nlohmann/json.hpp>
#include <map>

namespace Volt{
 
using json = nlohmann::json;

class PolyhedralTemplateMatchingEngine{
public:
    PolyhedralTemplateMatchingEngine(AnalysisContext& context, double rmsdCutoff);

    void perform();

    json buildMainListing() const;
    std::string getStructureTypeName(int structureType) const;

private:
    void storeOrientationData(const PTM::Kernel& kernel, size_t atomIndex);

    bool setupPTM(PTM& ptm) const;

    AnalysisContext& _context;
    double  _rmsd;
    mutable std::map<int, int> _structureStatistics;
    mutable bool _statisticsValid = false;
};

}