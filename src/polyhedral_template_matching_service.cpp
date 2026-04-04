#include <volt/analysis/analysis_context.h>
#include <volt/analysis/cluster_builder.h>
#include <volt/analysis/cluster_graph_export.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/core/analysis_result.h>
#include <volt/polyhedral_template_matching_service.h>
#include <volt/ptm_structure_analysis.h>
#include <volt/ptm_cluster_input_adapter.h>
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
    FrameAdapter::PreparedAnalysisInput prepared;
    std::string frameError;
    if(!FrameAdapter::prepareAnalysisInput(frame, prepared, &frameError)){
        return AnalysisResult::failure(frameError);
    }

    auto positions = std::move(prepared.positions);

    auto structureTypes = std::make_unique<ParticleProperty>(frame.natoms, DataType::Int, 1, 0, true);
    std::vector<Matrix3> preferredOrientations{Matrix3::Identity()};

    AnalysisContext context(
        positions.get(),
        frame.simulationCell,
        _inputCrystalStructure,
        nullptr,
        structureTypes.get(),
        std::move(preferredOrientations)
    );

    try{
        StructureAnalysis analysis(context);
        determineLocalStructuresWithPTM(analysis, _rmsd);
        computeMaximumNeighborDistanceFromPTM(analysis);
        PTMClusterInputAdapter clusterInputAdapter;
        clusterInputAdapter.prepare(analysis, context);
        ClusterBuilder clusterBuilder(analysis, context);
        clusterBuilder.build(_dissolveSmallClusters);

        json result = AnalysisResult::success();

        if(!outputBase.empty()){
            const std::string dumpPath = outputBase + "_annotated.dump";
            if(!context.writeDumpWithContext(frame, dumpPath, &analysis)){
                return AnalysisResult::failure("Failed to write " + dumpPath);
            }

            ClusterGraphExportPaths clusterGraphPaths;
            if(!exportClusterGraph(analysis.clusterGraph(), outputBase, &clusterGraphPaths)){
                return AnalysisResult::failure("Failed to export cluster graph tables");
            }

            result["annotated_dump"] = dumpPath;
            result["clusters_table"] = clusterGraphPaths.clustersTablePath;
            result["cluster_transitions_table"] = clusterGraphPaths.clusterTransitionsTablePath;
        }

        return result;
    }catch(const std::exception& error){
        return AnalysisResult::failure(std::string("PTM analysis failed: ") + error.what());
    }
}

}
