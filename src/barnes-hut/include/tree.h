#ifndef TREE_H
#define TREE_H

#include <math.h>
#include <vector>

#include "body.h"
#include "node.h"
#include "util.h"

class Octree {

    public:
        Node* root;
        scalar tsize; // simulation region, full scale (no units currently)
        vec corner;
        int n; // total number of particles in the simulation

        scalar kenergy;
        scalar penergy;

        Body** nbody; // list of all bodies in the simuation

        Octree(); // default constructor
        Octree( scalar cx, scalar cy, scalar cz, scalar dx); // used to construct the root node (full simulation area)

        ~Octree( ); // destructor (think about what needs to be cleaned up here...)
    
        void build_tree(int n, scalar *xi, scalar *yi, scalar *zi, scalar *vxi, scalar *vyi, scalar *vzi, scalar *mass);
        void compute_forces( scalar theta, scalar dt);
        void print_bodies( int step );
        void save_step( int step, scalar time, scalar theta, const char *run, const char *prefix );

    private: // to help us rebuild the tree during force calculations
        void delete_nodes( Node *node ); // recursive
        void rebuild_tree( );
};

#endif