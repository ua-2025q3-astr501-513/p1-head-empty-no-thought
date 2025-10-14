#include "tree.h"
#include "body.h"
#include "node.h"
#include "util.h"

#include <iostream>
#include <sstream>
#include <stdio.h>
#include <time.h>
#include <sys/stat.h>

Octree::Octree( ) {
    this->root = nullptr;
    this->tsize = 0;
    this->corner = {0, 0, 0};
    this->n = 0;

    this->kenergy = 0;
    this->penergy = 0;

} 

Octree::Octree( scalar cx, scalar cy, scalar cz, scalar dx ) {
    this->corner = {cx, cy, cz};

    this->root = new Node( corner, dx );
    this->tsize = dx;
    this->n = 0;

    this->kenergy = 0;
    this->penergy = 0;
} 
    
Octree::~Octree( ) {
    // figure out additional clean up later...
    delete[] this->nbody;
}

void Octree::build_tree(int n, scalar *xi, scalar *yi, scalar *zi, scalar *vxi, scalar *vyi, scalar *vzi, scalar *mass) {
    
    this->n = n; // updating the number of bodies in the simulation!
    this->nbody = new Body*[n];

    for (int i = 0; i < n; i++) {
        Body *b = new Body( xi[i], yi[i], zi[i], vxi[i], vyi[i], vzi[i], mass[i] * MSUN );
        nbody[i] = b; // adding this to our list of pointers
        root->insert( b ); // recursion in this function will take care of the rest.
    }
}


void Octree::delete_nodes(Node* node) {
    if (!node) return;
    for (int i = 0; i < 8; ++i) {
        if (node->children[i]) {
            delete_nodes(node->children[i]);
            // node->children[i] = nullptr;
        }
    }
    delete node;
}


void Octree::rebuild_tree() {
    delete root;
    // root = nullptr; // deleting our tree

// >>> scaling up (or down, i guess) our simulation size if needed.
    scalar farthest = 0; 
    vec origin = {0, 0, 0};
    int fidx = -1;

    for (int i = 0; i < n; i++) {
        if (distance(nbody[i]->pos, origin) > farthest) {
            fidx = i;
        }
    }

    scalar max_coord = fabs(nbody[fidx]->pos.x);
    if (fabs(nbody[fidx]->pos.y) > max_coord) max_coord = fabs(nbody[fidx]->pos.y);
    if (fabs(nbody[fidx]->pos.z) > max_coord) max_coord = fabs(nbody[fidx]->pos.z);
    
    // debugging remnant, dynamic resizing
    // std::cout << this->tsize << ", "<< max_coord << ", " << max_coord/(tsize/2) << "\n";

// >>> updating and rebuilding our tree!

    // resizing if needed
    if (max_coord > (tsize/2) * .3 ) {          // if our max coordinate is more than 30% of our simulation size
        this->tsize = max_coord * 20;           // makes the system 10^3 times larger
        this->corner =  { -this->tsize/2, -this->tsize/2, this->tsize/2};
    }
    else if ( max_coord < (tsize/2) * .10) {    // if our max coordinate is less than 15% of our simulation size
        this->tsize = max_coord / 5;            // makes the system 10^3 times smaller
        this->corner =  { -this->tsize/2, -this->tsize/2, this->tsize/2};
    }

    // creates a new root
    this->root = new Node( this->corner, this->tsize );

    // rebuilds the tree itself with the existing list of bodies.
    for (int i = 0; i < n; i++) {
        root->insert( this->nbody[i] ); // recursion in this function will take care of the rest.
    }
}

void Octree::compute_forces( scalar theta, scalar dt ) {

    // zeroing out our accelerations so they *don't* sum
    for(int i = 0; i < n; i++)
        nbody[i]->acc = {0,0,0};

    // force calculation, barnes-hut inside here!
    for (int i = 0; i < n; i++) 
        root->get_force( nbody[i], theta);

// >>> leapfrog integration here...

    // kick
    for(int i = 0; i < n; i++)
        nbody[i]->vel += nbody[i]->acc * 0.5 * dt;
    // drift
    for(int i = 0; i < n; i++)
        nbody[i]->pos += nbody[i]->vel * dt;

// >>> rebuilding our tree with updated postions
    rebuild_tree();

    // kick, again
    for(int i = 0; i < n; i++)
        nbody[i]->vel += nbody[i]->acc * 0.5 * dt;

    // fixme: update potential and kinetic energy calculation later, through simple loop
    // this->penergy = 0.0;
    // this->kenergy = 0.0;

}

// todo: overload this? need a better format for long simulation output files.
void Octree::print_bodies( int step) {

    std::cout << "timestep ::\t" << step << "\n";

    // add some nice formatting string to this later or something...
    // also headers! headers would be good.
    for (int i = 0; i < n; i++ ) {
        std::cout << i << "\t" << nbody[i]->pos.x << "\t" << nbody[i]->pos.y << "\t" << nbody[i]->pos.z << "\n";
    }

}

void Octree::save_step( int step, scalar step_time, scalar theta, const char *prefix, const char *run) {

    // we should check that our directory we want to make exists, if not we create it...
    char dname[50];
    char fname[100];
    sprintf(dname, "%s/%s", DATPATH, run);
    sprintf(fname, "%s/%s/globr_%s_%07d.dat", DATPATH, run, prefix, step);

    int status = mkdir(DATPATH, 0777);
    
    if (status <= 0 && errno != ENOENT) {
        status = mkdir(dname, 0777);

        if (status <= 0 && errno != ENOENT) {
            FILE *fout = fopen( fname, "w"); // new file to write to

            if ( fout != NULL ) { // checking that our file was made correctly

                // way too much work to get the time of our file's writing....
                time_t current_time;
                time(&current_time);
                struct tm * timeinfo;
                timeinfo = localtime(&current_time);
                char tstring[100];
                strftime( tstring, sizeof(char) * 100, "%X, %e %h %g", timeinfo);

                // a little header
                fprintf( fout, "# >>> globr_%s_%07d.dat. file written at %s.\n", prefix, step, tstring);
                // i should add something here that tells us initial conditions, perhaps.
                fprintf( fout, "# >>> timestep                  : %-15d\n", step);
                fprintf( fout, "# >>> particles                 : %-15d\n", n);
                fprintf( fout, "# >>> theta                     : %-15.3f\n", theta );
                fprintf( fout, "# >>> simulation time (yr)      : %-15.3e\n", step_time/YR);
                fprintf( fout, "# >>> simulation size (au)      : %-15.3e\n", tsize/AU);
                fprintf( fout, "# >>> kinetic energy (J)        : %-15.3e\n", kenergy);
                fprintf( fout, "# >>> potential energy (J)      : %-15.3e\n", penergy);
                fprintf( fout, "\n# >>> -------------------------------------------------------------------------------------------\n");
                fprintf( fout, "%-8s  %18s  %18s  %18s  %18s\n\n", "pID", "mass [msun]", "x [m]", "y [m]", "z [m]");
                // now onto the actual data!
                for (int i = 0; i < n; i++) {
                    fprintf(fout, "%-8d  %18.3f  %18.5e  %18.5e  %18.5e\n", i, nbody[i]->mass / MSUN, nbody[i]->pos.x, nbody[i]->pos.y, nbody[i]->pos.z);
                }

                fclose( fout );
            }
        }
    }

}