#ifndef BODY_H
#define BODY_H

#include <cmath>
#include "util.h"

class Body {

    public:
        vec pos;
        vec vel;
        vec acc;
        scalar mass;

        Body( );
        Body( scalar x, scalar y, scalar z, scalar vx, scalar vy, scalar vz, scalar m);

};

#endif