#ifndef UTIL_H
#define UTIL_H

#define DATPATH "../data"
#define G 6.67e-11      // gravitational constant [ m^3/kg/s] 
#define MSUN 2e30       // solar mass [kg]
#define AU 1.495E11     // astronomical unit [m]
#define PC 3.085E16     // parsec [m]
#define YR 3.15e7       // year [s]

typedef float scalar; // general data type, should be float or higher precision.

#include <cmath>

/**
 * helper struct that creates a vector in cartesian coordinates and 
 * includes basic mathematical operations (to make our lives easier).
 * 
 * @param x position, x [m]
 * @param y position, y [m]
 * @param z position, z [m]
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
 * 
 * @param a vector, cartesian coordinates
 * @param b vector, cartesian coordinates
 * 
 * @returns the Euclidean distance between a and b.
*/
inline scalar distance( const vec a, const vec b) {
    vec diff =  (a*a) - (b*b);
    return sqrt( diff.x + diff.y + diff.z);
}

/** 
 * calculates which quadrant/octant of a node's physical domain a 
 * child node should be placed in.
 * 
 * NOTE. quadrants are labeled starting at the top left quadrant (0) 
 * and numbered clockwise, then moving 'down' through the levels of the cube.
 * 
 * @param dx side length of the parent node
 * @param corner corner coordinates of the parent node
 * @param pos vector, position of the body that is being placed
 * 
 * @returns an integer from 0-7 indicating the correct octant.
*/
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

/** 
 * calculates the new corner coordinates of a child node.
 * 
 * NOTE. quadrants are labeled starting at the top left quadrant (0) 
 * and numbered clockwise, then moving 'down' through the levels of the cube.
 * 
 * @param q quadrant of the new node
 * @param corner corner coordinates of the parent node
 * @param dx side length of the parent node
 * 
 * @returns an vector with the correct corner coordinates for the child.
*/
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