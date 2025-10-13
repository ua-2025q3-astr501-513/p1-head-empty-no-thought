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
void Octree::print_bodies( int step ) {

    std::cout << "timestep ::\t" << step << "\n";

    // add some nice formatting string to this later or something...
    // also headers! headers would be good.
    for (int i = 0; i < n; i++ ) {
        std::cout << i << "\t" << nbody[i]->pos.x << "\t" << nbody[i]->pos.y << "\t" << nbody[i]->pos.z << "\n";
    }

}

void Octree::save_step( int step, scalar step_time, const char *prefix, const char *run) {

    // we should check that our directory we want to make exists, if not we create it...
    char dname[50];
    char fname[100];
    sprintf(dname, "%s/%s", DATPATH, run);
    sprintf(fname, "%s/%s/globr_%s_%d.dat", DATPATH, run, prefix, step);

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
                fprintf( fout, "# >>> globr_%s_%d.dat. file written at %s.\n", prefix, step, tstring);
                // i should add something here that tells us initial conditions, perhaps.
                fprintf( fout, "# >>> timestep               : %-15d\n", step);
                fprintf( fout, "# >>> simulation time (unit) : %-15.3e\n", step_time);
                fprintf( fout, "\n------------------------------------------------------------------------------------------------\n");
                fprintf( fout, "%-8s  %18s  %18s  %18s  %18s\n\n", "pID", "mass [msun]", "x [unit]", "y [unit]", "z [unit]");
                // now onto the actual data!
                for (int i = 0; i < n; i++) {
                    fprintf(fout, "%-8d  %18.5e  %18.5e  %18.5e  %18.5e\n", i, nbody[i]->mass, nbody[i]->pos.x, nbody[i]->pos.y, nbody[i]->pos.z);
                }

                fclose( fout );
            }
        }
    }

}