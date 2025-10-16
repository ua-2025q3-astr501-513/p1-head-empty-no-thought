#include "body.h"
#include "util.h"

/**
 * base constructor, initializes all fields to zero.
*/
Body::Body( ) {
    
    this->pos = {0, 0, 0};
    this->vel = {0, 0, 0};
    this->acc = {0, 0, 0};
    this->mass = 0.0;

}

/**
 * constructor for a new Body.
 * 
 * @param x initial position, x [m]
 * @param y initial position, y [m]
 * @param z initial position, z [m]
 * @param vx initial velocity, x [m]
 * @param vy initial velocity, y [m]
 * @param vz initial velocity, z [m]
 * @param m mass [kg]
*/
Body::Body( scalar x, scalar y, scalar z, scalar vx, scalar vy, scalar vz, scalar m ) {
    
    this->pos = {x, y, z};
    this->vel = {vx, vy, vz};
    this->acc = {0, 0, 0};
    this->mass = m;

}