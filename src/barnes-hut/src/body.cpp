#include "body.h"

Body::Body( ) {
    
    this->pos = { };
    this->vel = { };
    this->acc = { };
    this->mass = 0.0;

}


Body::Body( scalar x, scalar y, scalar z, scalar vx, scalar vy, scalar vz, scalar m) {
    
    this->pos = {x, y, z};
    this->vel = {vx, vy, vz};
    this->acc = { };
    this->mass = m;

}