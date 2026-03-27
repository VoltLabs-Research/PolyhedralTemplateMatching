#include <volt/ptm_cluster_builder.h>

#include <volt/polyhedral_template_matching.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <cmath>
#include <deque>
#include <vector>

namespace Volt{

namespace{

constexpr double kStrictOrientationThresholdDegrees = 3.0;
constexpr double kRelaxedOrientationThresholdDegrees = 8.0;

}

PTMClusterBuilder::PTMClusterBuilder(
    StructureAnalysis& sa,
    AnalysisContext& context,
    std::shared_ptr<const ParticleProperty> ptmRmsd,
    std::shared_ptr<const ParticleProperty> orientations
) : ClusterBuilder(sa, context)
  , _ptmRmsd(std::move(ptmRmsd))
  , _orientations(std::move(orientations)){}

bool PTMClusterBuilder::areOrientationsCompatible(int atom1, int atom2, int structureType){
    const Quaternion q1 = getAtomOrientation(atom1);
    const Quaternion q2 = getAtomOrientation(atom2);
    const Quaternion quaternionDifference = q1.inverse() * q2;

    const double rmsd1 = _ptmRmsd ? _ptmRmsd->getDouble(atom1) : 0.0;
    const double rmsd2 = _ptmRmsd ? _ptmRmsd->getDouble(atom2) : 0.0;
    const double averageRmsd = 0.5 * (rmsd1 + rmsd2);

    double thresholdDegrees = kRelaxedOrientationThresholdDegrees;
    if(structureType != static_cast<int>(StructureType::SC) && averageRmsd < 0.1){
        thresholdDegrees = kStrictOrientationThresholdDegrees;
    }

    const Matrix3 rotationMatrix = quaternionToMatrix(quaternionDifference);
    const double thresholdRadians = thresholdDegrees * M_PI / 180.0;
    const double minimumTrace = 1.0 + 2.0 * std::cos(thresholdRadians);

    for(int symmetryIndex = 0; symmetryIndex < _sa.symmetryPermutationCount(structureType); ++symmetryIndex){
        const Matrix3 product =
            rotationMatrix * _sa.symmetryTransformation(structureType, symmetryIndex).transposed();
        if((product(0, 0) + product(1, 1) + product(2, 2)) > minimumTrace){
            return true;
        }
    }

    return false;
}

Quaternion PTMClusterBuilder::getAtomOrientation(int atom) const{
    const double* orientation = _orientations->constDataDouble() + atom * 4;
    Quaternion quaternion(orientation[0], orientation[1], orientation[2], orientation[3]);
    quaternion.normalize();
    return quaternion;
}

void PTMClusterBuilder::grow(
    Cluster* cluster,
    std::deque<int>& atomsToVisit,
    int structureType
){
    while(!atomsToVisit.empty()){
        const int currentAtom = atomsToVisit.front();
        atomsToVisit.pop_front();

        const int neighborCount = _sa.numberOfNeighbors(currentAtom);
        for(int neighborIndex = 0; neighborIndex < neighborCount; ++neighborIndex){
            const int neighbor = _sa.getNeighbor(currentAtom, neighborIndex);
            if(neighbor < 0 || neighbor == currentAtom){
                continue;
            }
            if(_context.atomClusters->getInt(neighbor) != 0){
                continue;
            }
            if(_context.structureTypes->getInt(neighbor) != structureType){
                continue;
            }
            if(!areOrientationsCompatible(currentAtom, neighbor, structureType)){
                continue;
            }

            _context.atomClusters->setInt(neighbor, cluster->id);
            cluster->atomCount++;

            const Matrix3 clusterOrientation = cluster->orientation;
            const Matrix3 neighborOrientation = quaternionToMatrix(getAtomOrientation(neighbor));
            const Matrix3 localRotation = clusterOrientation.inverse() * neighborOrientation;
            const int symmetryIndex = _sa.findClosestSymmetryPermutation(structureType, localRotation);
            _context.atomSymmetryPermutations->setInt(neighbor, symmetryIndex);
            atomsToVisit.push_back(neighbor);
        }
    }
}

void PTMClusterBuilder::initOrientation(Cluster* cluster, size_t seedAtomIndex){
    cluster->orientation = quaternionToMatrix(getAtomOrientation(static_cast<int>(seedAtomIndex)));
}

void PTMClusterBuilder::reorientAtomsToAlign(){
    tbb::parallel_for(tbb::blocked_range<size_t>(0, _context.atomCount()), [this](const tbb::blocked_range<size_t>& r){
        for(size_t atomIndex = r.begin(); atomIndex != r.end(); ++atomIndex){
            int clusterId = _context.atomClusters->getInt(atomIndex);
            if(clusterId == 0){
                continue;
            }

            Cluster* cluster = _sa.clusterGraph().findCluster(clusterId);
            assert(cluster);
            if(cluster->symmetryTransformation == 0){
                continue;
            }

            int oldSymmetry = _context.atomSymmetryPermutations->getInt(atomIndex);
            int newSymmetry = _sa.symmetryInverseProduct(
                cluster->structure,
                oldSymmetry,
                cluster->symmetryTransformation
            );

            _context.atomSymmetryPermutations->setInt(atomIndex, newSymmetry);
        }
    });
}

bool PTMClusterBuilder::calculateMisorientation(
    int atomIndex,
    int neighbor,
    int neighborIndex,
    Matrix3& outTransition
) const{
    const int structureType = _context.structureTypes->getInt(atomIndex);
    const int neighborStructureType = _context.structureTypes->getInt(neighbor);
    if(structureType == static_cast<int>(StructureType::OTHER) ||
       neighborStructureType == static_cast<int>(StructureType::OTHER)){
        return false;
    }

    Matrix3 tm1, tm2;
    for(int i = 0; i < 3; ++i){
        int atomToMatch = atomIndex;
        if(i != 2){
            const int commonNeighborIndex = _sa.commonNeighborIndex(structureType, neighborIndex, i);
            if(commonNeighborIndex < 0){
                return false;
            }

            atomToMatch = _sa.getNeighbor(atomIndex, commonNeighborIndex);
            tm1.column(i) =
                _sa.latticeVector(
                    structureType,
                    _sa.symmetryPermutationEntry(
                        structureType,
                        _context.atomSymmetryPermutations->getInt(atomIndex),
                        commonNeighborIndex
                    )
                ) -
                _sa.latticeVector(
                    structureType,
                    _sa.symmetryPermutationEntry(
                        structureType,
                        _context.atomSymmetryPermutations->getInt(atomIndex),
                        neighborIndex
                    )
                );
        }else{
            tm1.column(i) = -_sa.latticeVector(
                structureType,
                _sa.symmetryPermutationEntry(
                    structureType,
                    _context.atomSymmetryPermutations->getInt(atomIndex),
                    neighborIndex
                )
            );
        }

        const int reverseNeighborIndex = _sa.findNeighbor(neighbor, atomToMatch);
        if(reverseNeighborIndex == -1){
            return false;
        }

        tm2.column(i) = _sa.latticeVector(
            neighborStructureType,
            _sa.symmetryPermutationEntry(
                neighborStructureType,
                _context.atomSymmetryPermutations->getInt(neighbor),
                reverseNeighborIndex
            )
        );
    }

    if(std::abs(tm1.determinant()) < EPSILON){
        return false;
    }

    Matrix3 tm1Inverse;
    if(!tm1.inverse(tm1Inverse)){
        return false;
    }

    outTransition = tm2 * tm1Inverse;
    return true;
}

void PTMClusterBuilder::connectClusters(){
    std::vector<std::vector<int>> extras(_context.atomCount());

    for(size_t atomIndex = 0; atomIndex < _context.atomCount(); ++atomIndex){
        int clusterId = _context.atomClusters->getInt(atomIndex);
        if(clusterId == 0){
            continue;
        }

        Cluster* cluster1 = _sa.clusterGraph().findCluster(clusterId);
        const int neighborCount = _sa.numberOfNeighbors(atomIndex);

        for(int neighborIndex = 0; neighborIndex < neighborCount; ++neighborIndex){
            const int neighbor = _sa.getNeighbor(static_cast<int>(atomIndex), neighborIndex);
            if(neighbor < 0 || neighbor == static_cast<int>(atomIndex)){
                continue;
            }

            const int neighborClusterId = _context.atomClusters->getInt(neighbor);
            if(neighborClusterId == 0){
                extras[neighbor].push_back(static_cast<int>(atomIndex));
                continue;
            }

            if(neighborClusterId == cluster1->id){
                continue;
            }

            Cluster* cluster2 = _sa.clusterGraph().findCluster(neighborClusterId);
            if(ClusterTransition* existing = cluster1->findTransition(cluster2)){
                existing->area++;
                existing->reverse->area++;
                continue;
            }

            Matrix3 transition;
            if(!calculateMisorientation(static_cast<int>(atomIndex), neighbor, neighborIndex, transition)){
                continue;
            }
            if(!transition.isOrthogonalMatrix()){
                continue;
            }

            ClusterTransition* transitionLink = _sa.clusterGraph().createClusterTransition(cluster1, cluster2, transition);
            transitionLink->area++;
            transitionLink->reverse->area++;
        }
    }

    _sa.appendNeighbors(extras);
    spdlog::info("Number of cluster transitions: {}", _sa.clusterGraph().clusterTransitions().size());
}

void PTMClusterBuilder::formSuperClusters(){
    const size_t oldTransitionCount = _sa.clusterGraph().clusterTransitions().size();

    for(Cluster* cluster : _sa.clusterGraph().clusters()){
        if(!cluster || cluster->id == 0){
            continue;
        }
        cluster->rank = 0;
    }

    for(Cluster* cluster : _sa.clusterGraph().clusters()){
        if(!cluster || cluster->id == 0){
            continue;
        }
        if(cluster->structure != _context.inputCrystalType){
            processDefectCluster(cluster);
        }
    }

    const size_t newTransitionCount = _sa.clusterGraph().clusterTransitions().size();

    for(size_t i = oldTransitionCount; i < newTransitionCount; ++i){
        ClusterTransition* transition = _sa.clusterGraph().clusterTransitions()[i];

        Cluster* parent1 = getParentGrain(transition->cluster1);
        Cluster* parent2 = getParentGrain(transition->cluster2);

        if(parent1 == parent2){
            continue;
        }

        ClusterTransition* parentTransition = transition;

        if(parent2 != transition->cluster2){
            parentTransition = _sa.clusterGraph().concatenateClusterTransitions(
                parentTransition,
                transition->cluster2->parentTransition
            );
        }

        if(parent1 != transition->cluster1){
            parentTransition = _sa.clusterGraph().concatenateClusterTransitions(
                transition->cluster1->parentTransition->reverse,
                parentTransition
            );
        }

        if(parent1->rank > parent2->rank){
            parent2->parentTransition = parentTransition->reverse;
            continue;
        }

        parent1->parentTransition = parentTransition;

        if(parent1->rank == parent2->rank){
            parent2->rank++;
        }
    }

    for(Cluster* cluster : _sa.clusterGraph().clusters()){
        getParentGrain(cluster);
    }
}

void PTMClusterBuilder::processDefectCluster(Cluster* defectCluster){
    for(ClusterTransition* transition = defectCluster->transitions; transition; transition = transition->next){
        if(transition->cluster2->structure != _context.inputCrystalType || transition->distance != 1){
            continue;
        }
        for(ClusterTransition* sibling = transition->next; sibling; sibling = sibling->next){
            if(sibling->cluster2->structure != _context.inputCrystalType || sibling->distance != 1){
                continue;
            }
            if(sibling->cluster2 == transition->cluster2){
                continue;
            }

            const Matrix3 misorientation = sibling->tm * transition->reverse->tm;
            for(int symIndex = 0; symIndex < _sa.symmetryPermutationCount(sibling->cluster2->structure); ++symIndex){
                if(_sa.symmetryTransformation(sibling->cluster2->structure, symIndex).equals(misorientation, 1e-6)){
                    _sa.clusterGraph().createClusterTransition(
                        transition->cluster2,
                        sibling->cluster2,
                        misorientation,
                        2
                    );
                    break;
                }
            }
        }
    }
}

void PTMClusterBuilder::build(bool dissolveSmallClusters){
    for(std::size_t seedAtomIndex = 0; seedAtomIndex < _context.atomCount(); ++seedAtomIndex){
        if(alreadyProcessedAtom(static_cast<int>(seedAtomIndex))){
            continue;
        }

        const int structureType = _context.structureTypes->getInt(seedAtomIndex);
        Cluster* cluster = startNew(static_cast<int>(seedAtomIndex), structureType);
        initOrientation(cluster, seedAtomIndex);
        cluster->symmetryTransformation = 0;
        _context.atomSymmetryPermutations->setInt(seedAtomIndex, 0);

        std::deque<int> atomsToVisit{static_cast<int>(seedAtomIndex)};
        grow(cluster, atomsToVisit, structureType);
    }

    reorientAtomsToAlign();
    connectClusters();
    formSuperClusters();
    if(dissolveSmallClusters){
        ClusterBuilder::dissolveSmallClusters();
    }
}

}
