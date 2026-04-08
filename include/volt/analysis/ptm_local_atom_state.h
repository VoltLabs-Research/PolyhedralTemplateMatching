#pragma once

#include <volt/math/quaternion.h>

namespace Volt{

struct PtmLocalAtomState{
    Quaternion orientation;
    double rmsd = 0.0;
    bool valid = false;

    PtmLocalAtomState() : orientation(Quaternion::Identity{}){}
};

}
