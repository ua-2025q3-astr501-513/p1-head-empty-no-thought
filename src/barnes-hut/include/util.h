#ifndef UTIL_H
#define UTIL_H

#define G 6.67e-8   // gravitational constant
typedef float scalar; // defining a 'scalar' data type

#include <cmath>

/**
 * helper struct that creates a vector with up to three coordinates and 
 * includes basic mathematical operations (to make our lives easier).
*/
struct vec {
    scalar x, y, z;

    vec() : x(0), y(0), z(0) {}
    vec( scalar X, scalar Y, scalar Z ): x(X), y(Y), z(Z) {}

    vec operator+(const vec &v )  const { return {x+v.x, y+v.y, z+v.z}; }
    vec operator-(const vec &v )  const { return {x-v.x, y-v.y, z-v.z}; }
    vec operator*( scalar f )  const { return {x*f, y*f, z*f}; }
    vec operator*( const vec &v )  const { return {x*v.x, y*v.y, z*v.z}; }
    vec& operator+=(const vec &v) { x+=v.x; y+=v.y; z+=v.z; return *this; }
    scalar norm() const { return std::sqrt(x*x + y*y + z*z); }
};

/**
 * returns the distance between two vectors.
*/
scalar distance( const vec a, const vec b) {
    vec diff =  (a*a) - (b*b);
    return sqrt( diff.x + diff.y + diff.z);
}

#endif