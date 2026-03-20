#include <volt/cli/common.h>
#include <volt/polyhedral_template_matching_service.h>

using namespace Volt;
using namespace Volt::CLI;

void showUsage(const std::string& name){
    printUsageHeader(name, "Volt - Polyhedral Template Matching");
    std::cerr
        << "  --rmsd <float>                RMSD threshold for PTM. [default: 0.1]\n"
        << "  --threads <int>               Max worker threads (TBB/OMP). [default: auto]\n";
    printHelpOption();
}

int main(int argc, char* argv[]){
    // TODO: hmmmmmmm... perhaps another approach for
    // defined required arguments, this hack is ugly
    if(argc < 2){
        showUsage(argv[0]);
        return 1;
    }

    std::string filename, outputBase;
    auto opts = parseArgs(argc, argv, filename, outputBase);

    if(hasOption(opts, "--help") || filename.empty()){
        showUsage(argv[0]);
        return filename.empty() ? 1 : 0;
    }

    auto parallel = initParallelism(opts, false);
    initLogging("volt-polyhedral-template-matching", parallel.threads);

    LammpsParser::Frame frame;
    if(!parseFrame(filename, frame)) return 1;

    outputBase = deriveOutputBase(filename, outputBase);
    spdlog::info("Output base: {}", outputBase);

    PolyhedralTemplateMatchingService analyzer;
    analyzer.setRMSD(getDouble(opts, "--rmsd", 0.1));

    spdlog::info("Starting PTM analysis...");
    json result = analyzer.compute(frame, outputBase);
    if(result.value("is_failed", false)){
        spdlog::error("Analysis failed: {}", result.value("error", "Unknown error"));
        return 1;
    }

    spdlog::info("PTM analysis completed.");
    return 0;
}
