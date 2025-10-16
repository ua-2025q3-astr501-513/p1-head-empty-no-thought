#include "tree.h"
#include "body.h"
#include "node.h"
#include "util.h"

#include <iostream>
#include <sstream>
#include <stdio.h>
#include <time.h>
#include <sys/stat.h>

/** basic constructor, initalizes to zero */
Octree::Octree( ) {
    this->root = nullptr;
    this->tsize = 0;
    this->corner = {0, 0, 0};
    this->n = 0;

    this->kenergy = 0;
    this->penergy = 0;

} 

/** constructor, root node and empty tree.
 * 
 * @param cx coorner coordinate, x [m]
 * @param cy coorner coordinate, y [m]
 * @param cz coorner coordinate, z [m]
 * @param dz total simulation domain, side length [m]
 * 
*/
Octree::Octree( scalar cx, scalar cy, scalar cz, scalar dx ) {
    this->corner = {cx, cy, cz};

    this->root = new Node( corner, dx );
    this->tsize = dx;
    this->n = 0;

    this->kenergy = 0;
    this->penergy = 0;
} 
    
/**
 * destructor.
*/
Octree::~Octree( ) {
    delete[] this->nbody;
}

/**
 * recursively populates the octree given the number of bodies and their initial conditions.
 * 
 * @param n number of bodies in the system
 * @param xi initial positions, x [m]
 * @param yi initial positions, y [m]
 * @param zi initial positions, z [m]
 * @param vxi initial velocities, x [m/s]
 * @param vyi initial velocities, y [m/s]
 * @param vzi initial velocities, z [m/s]
 * @param mass masses [solar mass]
 * 
*/
void Octree::build_tree(int n, scalar *xi, scalar *yi, scalar *zi, scalar *vxi, scalar *vyi, scalar *vzi, scalar *mass) {
    
    this->n = n; // updating the number of bodies in the simulation!
    this->nbody = new Body*[n];

    for (int i = 0; i < n; i++) {
        Body *b = new Body( xi[i], yi[i], zi[i], vxi[i], vyi[i], vzi[i], mass[i] * MSUN );
        nbody[i] = b; // adding this to our list of pointers
        root->insert( b ); // recursion in this function will take care of the rest.
    }
}

/**
 * reconstructs the tree and rescales total simulation domain (size) if needed. 
 * rescaling occurs if a body's distance from the origin is greater than 30% of the
 * total domain.
 * 
 * this currently only supports upscaling, downscaling has not been debugged and 
 * implemented.
*/
void Octree::rebuild_tree( ) {
    delete root;

// >>> scaling up our simulation size if needed.
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
    // todo: fix dynamic rescaling when making simulation domain smaller, segfaulting
    // else if ( max_coord < (tsize/2) * .10) {    // if our max coordinate is less than 15% of our simulation size
    //     this->tsize /= 10;            // makes the system 10^3 times smaller
    //     this->corner =  { -this->tsize/2, -this->tsize/2, this->tsize/2};
    // }

    // creates a new root
    this->root = new Node( this->corner, this->tsize );

    // rebuilds the tree itself with the existing list of bodies.
    for (int i = 0; i < n; i++) {
        root->insert( this->nbody[i] ); // recursion in this function will take care of the rest.
    }
}

/**
 * handles high-level force computations for all bodies in the tree and updates
 * positions, velocities, and accelerations through leapfrog integration. 
 * 
 * also calculates system energies (buggy) for basic diagnostics.
 * 
 * @param theta threshold criteria for barnes-hut, ratio of node width to distance to center of mass.
 * @param dt timestep [s]
*/
void Octree::compute_forces( scalar theta, scalar dt ) {

    // zeroing out our accelerations so they *don't* sum
    for (int i = 0; i < n; i++)
        nbody[i]->acc = {0,0,0};

    // force calculation, barnes-hut inside here!
    for (int i = 0; i < n; i++) 
        root->get_force( nbody[i], theta);

// >>> leapfrog integration here...

    // kick
    for (int i = 0; i < n; i++)
        nbody[i]->vel += nbody[i]->acc * 0.5 * dt;
    // drift
    for (int i = 0; i < n; i++)
        nbody[i]->pos += nbody[i]->vel * dt;

// >>> rebuilding our tree with updated postions
    rebuild_tree();

    // kick, again
    kenergy = 0.0; // zeroing our previous kinetic
    for (int i = 0; i < n; i++) {
        nbody[i]->vel += nbody[i]->acc * 0.5 * dt;
        kenergy += 0.5 * nbody[i]->mass * nbody[i]->vel.norm(); // sneaking in a quick kinetic energy calculation
    }

    // todo: fix softening term here, should be on the scale of 0.01 pc for globular clusters.
    penergy = 0.0; // zeroing our previous potential
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if ( distance(nbody[i]->pos, nbody[j]->pos) > 1 * AU) // softening term, not sure if this is enough
                penergy += G * nbody[i]->mass * nbody[j]->mass / distance(nbody[i]->pos, nbody[j]->pos);
        }
    }

}

/**
 * Prints basic information about every particle in a system. 
 * Only good for VERY small systems, don't use unless debugging.
 * 
 * @param step integer timestep of the simulation
 */ 
void Octree::print_bodies( int step) {

    std::cout << "timestep ::\t" << step << "\n";

    for (int i = 0; i < n; i++ ) {
        std::cout << i << "\t" << nbody[i]->pos.x << "\t" << nbody[i]->pos.y << "\t" << nbody[i]->pos.z << "\n";
    }

}

/**
 * handles data output for the entire simulation and makes sure our .dat files 
 * are formatted to be pretty and machine readable.
 * 
 * all files have the naming convention globr_{run}_0000000.dat.
 * 
 * @param step integer timestep
 * @param step_time physical time of timestep [s]
 * @param theta threshold criterion
 * @param run name of simulation run, for file naming
*/
void Octree::save_step( int step, scalar step_time, scalar theta, const char *run) {

    // we should check that our directory we want to make exists, if not we create it...
    char dname[50];
    char fname[100];
    sprintf(dname, "%s/%s", DATPATH, run);
    sprintf(fname, "%s/%s/globr_%s_%07d.dat", DATPATH, run, run, step);

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
                fprintf( fout, "# >>> globr_%s_%07d.dat. file written at %s.\n", run, step, tstring);
                // i should add something here that tells us initial conditions, perhaps.
                fprintf( fout, "# >>> timestep                  : %-15d\n", step);
                fprintf( fout, "# >>> particles                 : %-15d\n", n);
                fprintf( fout, "# >>> theta                     : %-15.3f\n", theta );
                fprintf( fout, "# >>> simulation time   [yr]    : %-15.3e\n", step_time/YR);
                fprintf( fout, "# >>> simulation size   [pc]    : %-15.3e\n", tsize/PC);
                fprintf( fout, "# >>> kinetic energy    [J]     : %-15.3e\n", kenergy);
                fprintf( fout, "# >>> potential energy  [J]     : %-15.3e\n", penergy);
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