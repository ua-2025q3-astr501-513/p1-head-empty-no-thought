#ifndef TREE_H
#define TREE_H

#include <math.h>
#include <vector>

#include "body.h"
#include "node.h"
#include "util.h"

class Octree {

    public:

        Node* root; /** pointer, root node*/
        scalar tsize; /** total simulation domain [m] */
        vec corner; /** coordinates of upper left corner, simulation domain [m] */
        int n; /** total number of particles in the simulation */

        scalar kenergy; /** total kinetic energy [J] */
        scalar penergy; /** total potential energy [J] */

        Body** nbody; /** list of pointers to all bodies in the simuation */

        Octree(); // default constructor
        Octree( scalar cx, scalar cy, scalar cz, scalar dx); // used to construct the root node (full simulation area)

        ~Octree( ); // destructor
    
        void build_tree(int n, scalar *xi, scalar *yi, scalar *zi, scalar *vxi, scalar *vyi, scalar *vzi, scalar *mass);
        void compute_forces( scalar theta, scalar dt);
        void print_bodies( int step );
        void save_step( int step, scalar time, scalar theta, const char *run );

    private: // to help us rebuild the tree during force calculations
        void rebuild_tree( );
};

#endif