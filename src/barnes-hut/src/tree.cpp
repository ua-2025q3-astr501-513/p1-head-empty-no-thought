#include "tree.h"
#include "body.h"
#include "node.h"
#include "util.h"
#include <iostream>

Octree::Octree( ) {
    this->root = nullptr;
    this->tsize = 0;
    this->n = 0;

} 

Octree::Octree( scalar cx, scalar cy, scalar cz, scalar dx ) {
    vec corner = {cx, cy, cz};

    this->root = new Node( corner, dx );
    this->tsize = dx;
    this->n = 0;
} 
    
Octree::~Octree( ) {
    // figure out additional clean up later...
    free( this->nbody );
}

void Octree::build_tree(int n, scalar *xi, scalar *yi, scalar *zi, scalar *vxi, scalar *vyi, scalar *vzi, scalar *mass) {
    
    this->n = n; // updating the number of bodies in the simulation!
    this->nbody = (Body **) malloc( n * sizeof( Body * )); // allocating memory for the array of pointers to all bodies

    for (int i = 0; i < n; i++) {
        Body *b = new Body( xi[i], yi[i], zi[i], vxi[i], vyi[i], vzi[i], mass[i] );
        root->insert( b ); // recursion in this function will take care of the rest.
        nbody[i] = b; // adding this to our list of pointers
    }
}

void Octree::compute_forces( float theta, float t ) {
    
    for (int i = 0; i < n; i++) {
        root->get_force( nbody[i], theta);
    }

    // todo: come back to this to optimize? leapfrog maybe
    for( int i = 0; i < n; i++) {

        nbody[i]->pos += (nbody[i]->vel * t) + (nbody[i]->acc * 0.5 * (t * t) );
        nbody[i]->vel += nbody[i]->acc * t;
    }
    

}

// todo: overload this? need a better format for long simulation output files.
void Octree::print_bodies( scalar step ) {

    std::cout << "timestep ::\t" << step << "\n";

    // add some nice formatting string to this later or something...
    // also headers! headers would be good.
    for (int i = 0; i < n; i++ ) {
        std::cout << i << "\t" << nbody[i]->pos.x << "\t" << nbody[i]->pos.y << "\t" << nbody[i]->pos.z << "\n";
    }

}