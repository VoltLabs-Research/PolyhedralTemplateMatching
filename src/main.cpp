#include <volt/cli/common.h>
#include <volt/analysis/ptm_service.h>
#include <volt/structures/crystal_structure_types.h>

using namespace Volt;
using namespace Volt::CLI;

void showUsage(const std::string& name){
    printUsageHeader(name, "Volt - Polyhedral Template Matching");
    std::cerr
        << "  --crystalStructure <type>     Crystal structure. (SC|FCC|HCP|BCC|CUBIC_DIAMOND|HEX_DIAMOND) [default: FCC]\n"
        << "  --rmsd <float>                RMSD threshold for PTM. [default: 0.1]\n"
        << "  --no-write                    Run analysis without writing dump or cluster tables.\n"
        << "  --dissolveSmallClusters       Mark small clusters as OTHER after building clusters.\n";
    printHelpOption();
}

int main(int argc, char* argv[]){
    std::string filename, outputBase;
    auto opts = parseArgs(argc, argv, filename, outputBase);
    if(const int startupStatus = handleHelpOrMissingInput(argc, argv, opts, filename, showUsage);
       startupStatus >= 0){
        return startupStatus;
    }

    initLogging("volt-polyhedral-template-matching");

    LammpsParser::Frame frame;
    if(!parseFrame(filename, frame)) return 1;

    const bool noWrite = hasOption(opts, "--no-write");
    if(noWrite){
        outputBase.clear();
        spdlog::info("Output writing disabled for this run.");
    }else{
        outputBase = deriveOutputBase(filename, outputBase);
        spdlog::info("Output base: {}", outputBase);
    }

    PolyhedralTemplateMatchingService analyzer;
    LatticeStructureType crystalStructure = LATTICE_FCC;
    const std::string crystalStructureOption = getString(opts, "--crystalStructure", "FCC");
    if(!parseLatticeStructureType(crystalStructureOption, crystalStructure)){
        spdlog::warn("Unknown crystal structure '{}', defaulting to FCC.", crystalStructureOption);
        crystalStructure = LATTICE_FCC;
    }
    analyzer.setInputCrystalStructure(crystalStructure);
    analyzer.setRMSD(getDouble(opts, "--rmsd", 0.1));
    analyzer.setDissolveSmallClusters(hasOption(opts, "--dissolveSmallClusters"));

    spdlog::info("Starting PTM analysis...");
    json result = analyzer.compute(frame, outputBase);
    if(result.value("is_failed", false)){
        spdlog::error("Analysis failed: {}", result.value("error", "Unknown error"));
        return 1;
    }

    spdlog::info("PTM analysis completed.");
    return 0;
}
