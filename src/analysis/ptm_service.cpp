#include <volt/analysis/reconstructed_analysis_pipeline.h>
#include <volt/analysis/structure_analysis_context.h>
#include <volt/analysis/cluster_graph_builder.h>
#include <volt/analysis/cluster_graph_io.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/core/analysis_result.h>
#include <volt/analysis/ptm_service.h>
#include <volt/analysis/ptm_structure_analysis.h>
#include <volt/analysis/ptm_cluster_input_adapter.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/particle_property.h>
#include <volt/structures/crystal_structure_types.h>

namespace Volt{

PolyhedralTemplateMatchingService::PolyhedralTemplateMatchingService()
    : _inputCrystalStructure(LATTICE_FCC)
    , _rmsd(0.1)
    , _dissolveSmallClusters(false){}

void PolyhedralTemplateMatchingService::setInputCrystalStructure(LatticeStructureType structureType){
    _inputCrystalStructure = structureType;
}

void PolyhedralTemplateMatchingService::setRMSD(double rmsd){
    _rmsd = rmsd;
}

void PolyhedralTemplateMatchingService::setDissolveSmallClusters(bool dissolveSmallClusters){
    _dissolveSmallClusters = dissolveSmallClusters;
}

json PolyhedralTemplateMatchingService::compute(
    const LammpsParser::Frame& frame,
    const std::string& outputBase
){
    std::string frameError;
    auto session = AnalysisPipelineUtils::prepareAnalysisSession(
        frame,
        _inputCrystalStructure,
        &frameError
    );
    if(!session){
        return AnalysisResult::failure(frameError);
    }
    AnalysisContext& context = session->context;

    try{
        StructureAnalysis analysis(context);
        determineLocalStructuresWithPTM(analysis, _rmsd);
        computeMaximumNeighborDistanceFromPTM(analysis);
        PTMClusterInputAdapter clusterInputAdapter;
        clusterInputAdapter.prepare(analysis, context);
        ClusterBuilder clusterBuilder(analysis, context);
        clusterBuilder.build(_dissolveSmallClusters);

        json result = AnalysisResult::success();

        if(!AnalysisPipelineUtils::appendClusterOutputs(
            frame,
            outputBase,
            context,
            analysis,
            result,
            &frameError
        )){
            return AnalysisResult::failure(frameError);
        }

        return result;
    }catch(const std::exception& error){
        return AnalysisResult::failure(std::string("PTM analysis failed: ") + error.what());
    }
}

}
