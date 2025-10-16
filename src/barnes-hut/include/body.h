#ifndef BODY_H
#define BODY_H

#include <cmath>
#include "util.h"

/**
 * a simple class that represents a body (star) in the simulation. contains information on 
 * position, velocity, acceleration, and mass.
*/
class Body {

    public:
        /** position in cartesian space, m/s */
        vec pos;
        /** velocity in cartesian space, m/s */
        vec vel;
        /** acceleration in cartesian space, m/s */
        vec acc;
        /** mass, kg */
        scalar mass;

        Body( );
        Body( scalar x, scalar y, scalar z, scalar vx, scalar vy, scalar vz, scalar m );

};

#endif