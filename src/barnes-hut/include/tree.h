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
        int n; // total number of particles in the simulation
        Body** nbody; // list of all bodies in the simuation

        Octree(); // default constructor
        Octree( scalar cx, scalar cy, scalar cz, scalar dx); // used to construct the root node (full simulation area)

        ~Octree( ); // destructor (think about what needs to be cleaned up here...)
    
        void build_tree(int n, scalar *xi, scalar *yi, scalar *zi, scalar *vxi, scalar *vyi, scalar *vzi, scalar *mass);
        void print_bodies( scalar step );

};

#endif