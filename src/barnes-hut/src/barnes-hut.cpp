#include "tree.h"
#include "body.h"
#include "node.h"
#include "util.h"

int main( int argc, char *argv[] ) {

    int n = 10; // test suite of particles.

    scalar size = 5 * AU; // total simulation size
    vec c = { -size/2, -size/2, size/2};
    Octree *bhtree = new Octree( c.x, c.y, c.z, size ); // initializing our tree

    // initial positions (centered roughly around origin, range -2.5 to 2.5)
    scalar x[10]  = {  1.42, -1.35,  2.10, -1.94,  0.46, -2.23,  1.78, -0.34,  2.42, -1.32 };
    scalar y[10]  = { -1.15,  0.74, -2.02,  1.83, -0.58,  2.37, -0.91,  1.21, -1.76,  0.27 };
    scalar z[10]  = {  0.91, -1.91,  1.27, -0.44,  2.19, -1.68,  0.33,  1.56, -2.25,  0.82 };

    // small initial velocities (sum ≈ 0, range -0.05 to 0.05)
    scalar vx[10] = {  0.018, -0.037,  0.022,  0.005, -0.041,  0.017, -0.026,  0.033, -0.011,  0.020 };
    scalar vy[10] = { -0.008,  0.014, -0.043,  0.019,  0.007, -0.032,  0.025, -0.017,  0.038, -0.003 };
    scalar vz[10] = {  0.045, -0.012,  0.003, -0.027,  0.031, -0.009, -0.023,  0.041, -0.035, -0.014 };

    scalar m[10]  = { 0.001, 0.1, 1.2, 3.1, 0.5, 0.7, 0.008, 0.094, 1.6, 0.2}; 
    bhtree->build_tree(n, x, y, z, vx, vy, vz, m);

    bhtree->print_bodies( 0.0 );

    return 0;
        
}