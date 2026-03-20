#include <volt/polyhedral_template_matching_service.h>

#include <volt/analysis/analysis_context.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/core/analysis_result.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/particle_property.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/utilities/json_utils.h>
#include <spdlog/spdlog.h>
#include <algorithm>

using namespace Volt;
using namespace Volt::Particles;

namespace{

json buildPtmPerAtomProperties(
    const AnalysisContext& context,
    const LammpsParser::Frame& frame,
    const std::vector<int>& structureTypes
){
    const auto ptmOrientation = context.ptmOrientation;
    const auto correspondences = context.correspondencesCode;

    json perAtom = json::array();

    for(size_t index = 0; index < static_cast<size_t>(frame.natoms); ++index){
        json atom;
        atom["id"] = index < frame.ids.size()
            ? frame.ids[index]
            : static_cast<int>(index);
        atom["structure_type"] = structureTypes[index];

        if(index < frame.positions.size()){
            const auto& pos = frame.positions[index];
            atom["pos"] = {pos.x(), pos.y(), pos.z()};
        }else{
            atom["pos"] = {0.0, 0.0, 0.0};
        }

        if(correspondences){
            atom["correspondence"] = static_cast<std::uint64_t>(correspondences->getInt64(index));
        }

        if(ptmOrientation){
            atom["orientation"] = {
                ptmOrientation->getDoubleComponent(index, 0),
                ptmOrientation->getDoubleComponent(index, 1),
                ptmOrientation->getDoubleComponent(index, 2),
                ptmOrientation->getDoubleComponent(index, 3)
            };
        }

        perAtom.push_back(atom);
    }

    return perAtom;
}

json buildAtomsExport(
    const LammpsParser::Frame& frame,
    const std::vector<int>& structureTypes,
    const StructureAnalysis& analysis
){
    constexpr int K = static_cast<int>(StructureType::NUM_STRUCTURE_TYPES);

    std::vector<std::string> names(K);
    for(int st = 0; st < K; ++st){
        names[st] = analysis.getStructureTypeName(st);
    }

    std::vector<std::vector<size_t>> structureAtomIndices(K);
    for(size_t i = 0; i < static_cast<size_t>(frame.natoms); ++i){
        const int raw = structureTypes[i];
        const int st = (0 <= raw && raw < K) ? raw : static_cast<int>(StructureType::OTHER);
        structureAtomIndices[static_cast<size_t>(st)].push_back(i);
    }

    std::vector<int> structureOrder;
    structureOrder.reserve(K);
    for(int st = 0; st < K; ++st){
        if(!structureAtomIndices[static_cast<size_t>(st)].empty()){
            structureOrder.push_back(st);
        }
    }

    std::sort(structureOrder.begin(), structureOrder.end(), [&](int a, int b){
        return names[a] < names[b];
    });

    json atomsByStructure = json::object();
    for(int st : structureOrder){
        json atomsArray = json::array();
        for(size_t atomIndex : structureAtomIndices[static_cast<size_t>(st)]){
            const bool hasPosition = atomIndex < frame.positions.size();
            const auto atomId = atomIndex < frame.ids.size()
                ? frame.ids[atomIndex]
                : static_cast<int>(atomIndex);

            if(hasPosition){
                const auto& pos = frame.positions[atomIndex];
                atomsArray.push_back({
                    {"id", atomId},
                    {"pos", {pos.x(), pos.y(), pos.z()}}
                });
            }else{
                atomsArray.push_back({
                    {"id", atomId},
                    {"pos", {0.0, 0.0, 0.0}}
                });
            }
        }
        atomsByStructure[names[st]] = atomsArray;
    }

    json exportWrapper;
    exportWrapper["export"] = json::object();
    exportWrapper["export"]["AtomisticExporter"] = atomsByStructure;
    return exportWrapper;
}

}

PolyhedralTemplateMatchingService::PolyhedralTemplateMatchingService(): _rmsd(0.1){}

void PolyhedralTemplateMatchingService::setRMSD(double rmsd){
    _rmsd = rmsd;
}

json PolyhedralTemplateMatchingService::compute(
    const LammpsParser::Frame& frame,
    const std::string& outputBase
){
    if(frame.natoms <= 0){
        return AnalysisResult::failure("Invalid number of atoms");
    }

    if(!FrameAdapter::validateSimulationCell(frame.simulationCell)){
        return AnalysisResult::failure("Invalid simulation cell");
    }

    auto positions = FrameAdapter::createPositionPropertyShared(frame);
    if(!positions){
        return AnalysisResult::failure("Failed to create position property");
    }

    auto structureTypes = std::make_unique<ParticleProperty>(frame.natoms, DataType::Int, 1, 0, true);
    std::vector<Matrix3> preferredOrientations{Matrix3::Identity()};

    AnalysisContext context(
        positions.get(),
        frame.simulationCell,
        LATTICE_FCC,
        nullptr,
        structureTypes.get(),
        std::move(preferredOrientations)
    );

    try{
        StructureAnalysis analysis(
            context,
            false,
            StructureAnalysis::Mode::PTM,
            static_cast<float>(_rmsd)
        );
        analysis.determineLocalStructuresWithPTM();
        analysis.computeMaximumNeighborDistanceFromPTM();

        std::vector<int> atomStructureTypes(
            static_cast<size_t>(frame.natoms),
            static_cast<int>(StructureType::OTHER)
        );

        for(int atomIndex = 0; atomIndex < frame.natoms; ++atomIndex){
            atomStructureTypes[static_cast<size_t>(atomIndex)] = context.structureTypes->getInt(atomIndex);
        }

        json result;
        result["main_listing"] = analysis.buildMainListing();
        result["per-atom-properties"] = buildPtmPerAtomProperties(context, frame, atomStructureTypes);

        if(!outputBase.empty()){
            const std::string msgpackPath = outputBase + "_polyhedral_template_matching.msgpack";
            if(!JsonUtils::writeJsonMsgpackToFile(result, msgpackPath, false)){
                return AnalysisResult::failure("Failed to write " + msgpackPath);
            }

            const std::string atomsPath = outputBase + "_atoms.msgpack";
            if(!JsonUtils::writeJsonMsgpackToFile(buildAtomsExport(frame, atomStructureTypes, analysis), atomsPath, false)){
                return AnalysisResult::failure("Failed to write " + atomsPath);
            }
        }

        result["is_failed"] = false;
        return result;
    }catch(const std::exception& error){
        return AnalysisResult::failure(std::string("PTM analysis failed: ") + error.what());
    }
}
