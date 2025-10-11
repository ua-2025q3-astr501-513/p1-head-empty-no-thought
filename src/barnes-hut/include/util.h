#ifndef UTIL_H
#define UTIL_H

#define G 6.67e-11      // m^3 * kg^-1 * s^-1 gravitational constant
#define MSUN 2e30       // kg, solar mass
#define AU 1.495E11     // astronomical unit, meters
#define PC 3.085E16     // parsec, meters

typedef float scalar; // defining a 'scalar' data type

#include <cmath>

/**
 * helper struct that creates a vector with up to three coordinates and 
 * includes basic mathematical operations (to make our lives easier).
*/
struct vec {
    scalar x, y, z;

    // constructors
    vec() : x(0), y(0), z(0) {}
    vec( scalar X, scalar Y, scalar Z ): x(X), y(Y), z(Z) {}

    // math helpers!
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
inline scalar distance( const vec a, const vec b) {
    vec diff =  (a*a) - (b*b);
    return sqrt( diff.x + diff.y + diff.z);
}

// NOTE. quadrants are labeled starting at the top left quadrant (0) 
// and numbered clockwise, then moving 'down' through the levels of the cube
inline int get_quadrant( scalar dx, vec corner, vec pos ) {
    int q =  0;
    scalar ndx = dx / 2;

    // checking x and y first, easy to handle z later
    if ( pos.x <= corner.x + ndx) {
        if (pos.y <= corner.y + ndx) 
            q = 0;
        else 
            q = 1;
    } else {
        if (pos.y <= corner.y + ndx) 
            q = 3;
        else 
            q = 2;
    }

    // now we do z in a simpler way
    if (pos.z > corner.z + ndx) 
        q += 4;
    
    return q;
}

inline vec get_new_corner( int q, vec corner, scalar dx) {
    scalar ndx = dx / 2;
    vec ncorner = corner;

    if (q % 4 == 1) {
        ncorner.y += ndx;
    } else if (q % 4 == 2) {
        ncorner.y += ndx;
        ncorner.x += ndx;
    } else if (q % 4 == 3) {
        ncorner.x += ndx;
    }

    if (q > 3) {
        ncorner.z += ndx;
    }

    return ncorner;
}

#endif