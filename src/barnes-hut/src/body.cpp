#include "body.h"
#include "util.h"

Body::Body( ) {
    
    this->pos = {0, 0, 0};
    this->vel = {0, 0, 0};
    this->acc = {0, 0, 0};
    this->mass = 0.0;

}


Body::Body( scalar x, scalar y, scalar z, scalar vx, scalar vy, scalar vz, scalar m ) {
    
    this->pos = {x, y, z};
    this->vel = {vx, vy, vz};
    this->acc = {0, 0, 0};
    this->mass = m;

}