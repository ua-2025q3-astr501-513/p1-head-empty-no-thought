#include "tree.h"

Octree::Octree( ) {
    this->root = nullptr;
    this->tsize = 0;

} 

Octree::Octree( scalar cx, scalar cy, scalar cz, scalar dx ) {
    vec corner = {cx, cy, cz};

    this->root = new Node( corner, dx );
    this->tsize = dx;
} 
    
Octree::~Octree( ) {
    // figure out clean up later...
}

int Octree::build_tree(int n, scalar *xi, scalar *yi, scalar *zi, scalar *vxi, scalar *viy, scalar *vzi) {
    // do this later....
}