#pragma once

#include <volt/core/volt.h>
#include <nlohmann/json.hpp>
#include <volt/core/lammps_parser.h>
#include <string>

namespace Volt{

using json = nlohmann::json;

class PolyhedralTemplateMatchingService{
public:
    PolyhedralTemplateMatchingService();

    void setRMSD(double rmsd);
    
    json compute(
        const LammpsParser::Frame& frame,
        const std::string& outputBase = ""
    );

private:
    double _rmsd;
};
    
}
