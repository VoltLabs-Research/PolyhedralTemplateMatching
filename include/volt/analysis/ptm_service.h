#pragma once

#include <volt/core/volt.h>
#include <nlohmann/json.hpp>
#include <volt/core/lammps_parser.h>
#include <volt/structures/crystal_structure_types.h>
#include <string>

namespace Volt{

using json = nlohmann::json;

class PolyhedralTemplateMatchingService{
public:
    PolyhedralTemplateMatchingService();

    void setInputCrystalStructure(LatticeStructureType structureType);
    void setRMSD(double rmsd);
    void setDissolveSmallClusters(bool dissolveSmallClusters);
    
    json compute(
        const LammpsParser::Frame& frame,
        const std::string& outputBase = ""
    );

private:
    LatticeStructureType _inputCrystalStructure;
    double _rmsd;
    bool _dissolveSmallClusters;
};
    
}
