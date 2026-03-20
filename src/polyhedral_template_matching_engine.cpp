#include <volt/polyhedral_template_matching_engine.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

namespace Volt{

PolyhedralTemplateMatchingEngine::PolyhedralTemplateMatchingEngine(AnalysisContext& context, double rmsd)
    : _context(context)
    , _rmsd(rmsd){}

bool PolyhedralTemplateMatchingEngine::setupPTM(PTM& ptm) const{
    // TODO: this should be an argument => "--calculte-def-gradient"
    ptm.setCalculateDefGradient(false);
    ptm.setRmsdCutoff(std::numeric_limits<double>::infinity());
    return ptm.prepare(_context.positions->constDataPoint3(), _context.atomCount(), _context.simCell);
}

void PolyhedralTemplateMatchingEngine::storeOrientationData(const PTM::Kernel& kernel, size_t atomIndex){
    auto quaternion = kernel.orientation();
    double* orientation = _context.orientation->dataDouble() + 4 * atomIndex;

    orientation[0] = quaternion.x();
    orientation[1] = quaternion.y();
    orientation[2] = quaternion.z();
    orientation[3] = quaternion.w();
}

// Performs structural identification
void PolyhedralTemplateMatchingEngine::perform(){
    const size_t N = _context.atomCount();
    if(N == 0){
        invalidateStatistics();
        return;
    }

    PTM ptm;
    if(!setupPTM(ptm)){
        throw std::runtime_error("Error trying to initialize PTM.");
    }

    _context.orientation = std::make_shared<ParticleProperty>(N, ParticleProperty::OrientationProperty, 4, true);
    _context.correspondences = std::make_shared<ParticleProperty>(N, DataType::Int64, 1, 0, true);

    std::fill(
        _context.structureTypes->dataInt(),
        _context.structureTypes->dataInt() + _context.structureTypes->size(),
        static_cast<int>(StructureType::OTHER)
    );

    std::vector<uint64_t> cached(N, 0ull);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, N), [&](const auto& r){
        PTM::Kernel kernel(ptm);
        for(size_t i = r.begin(); i < r.end(); ++i){
            kernel.cacheNeighbors(i, &cached[i]);
        }
    });

    tbb::parallel_for(tbb::blocked_range<size_t>(0, N), [&](const auto& r){
        PTM::Kernel kernel(ptm);
        auto* rawCorrespondence = reinterpret_cast<uint64_t*>(_context.correspondences->data());
        for(size_t i = r.begin(); i < r.end(); ++i){
            StructureType type = kernel.identifyStructure(i, cached);
            const double rmsd = kernel.rmsd();
            //->setDouble(i, rmsd);
            rawCorrespondence[i] = kernel.correspondencesCode();

            if(type == StructureType::OTHER || rmsd > _rmsdCutoff){
                continue;
            }

            _context.structureTypes->setInt(i, type);
            storeOrientationData(kernel, i);
        }
    });
}

std::string PolyhedralTemplateMatchingEngine::getStructureTypeName(int structureType) const{
    return std::string(structureTypeName(structureType));
}

json PolyhedralTemplateMatchingEngine::buildMainListing() const{
    const size_t N = _context.atomCount();
    for(size_t i = 0; i < N; ++i){
        int structureType = _context.structureTypes->getInt(i);
        _structureStatistics[structureType]++;
    }

    const int N = static_cast<int>(_context.atomCount());
    const double invN = (N > 0) ? (100.0 / static_cast<double>(N)) : 0.0;

    int totalIdentified = 0;
    int unidentified = 0;
    auto itOther = _structureStatistics.find(static_cast<int>(StructureType::OTHER));
    if(itOther != _structureStatistics.end()){
        unidentified = itOther->second;
    }

    json mainListing = json::object();
    mainListing["total_atoms"] = N;
    mainListing["analysis_method"] = "PTM";

    for(const auto& [structureType, count] : _structureStatistics){
        std::string name = getStructureTypeName(structureType);
        std::transform(name.begin(), name.end(), name.begin(), ::tolower);
        mainListing[name + "_count"] = count;
        mainListing[name + "_percentage"] = static_cast<double>(count) * invN;
        if(structureType != static_cast<int>(StructureType::OTHER) &&
           structureType != static_cast<int>(CoordinationStructureType::COORD_OTHER)){
            totalIdentified += count;
        }
    }

    mainListing["total_identified"] = totalIdentified;
    mainListing["total_unidentified"] = unidentified;
    mainListing["identification_rate"] = static_cast<double>(totalIdentified) * invN;
    mainListing["unique_structure_types"] = static_cast<int>(_structureStatistics.size());

    return mainListing;
}

}