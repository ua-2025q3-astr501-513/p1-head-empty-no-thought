#include "tree.h"
#include "body.h"
#include "node.h"
#include "util.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstring>

#define DATAPATH "../data/"
#define INITPATH "../init/"

struct config {
    int n = 0;
    int size = 10;
    scalar dt = 1;
    scalar theta = 0.5;
    int nstep = 5000;
    int fout = nstep / 1000;
    char* run = nullptr;
    char* prefix = nullptr;
    char* filename = nullptr;
};

config parse_args(int argc, char** argv) {
    config cfg;

    for (int i = 1; i < argc; ++i) {
        if ((std::strcmp(argv[i], "-N") == 0 || std::strcmp(argv[i], "--N") == 0) && i + 1 < argc) {
            cfg.n = std::atoi(argv[++i]);
        } else if (std::strcmp(argv[i], "--size") == 0 && i + 1 < argc) {
            cfg.size = std::atoi(argv[++i]);
        } else if (std::strcmp(argv[i], "--step") == 0 && i + 1 < argc) {
            cfg.dt = std::atof(argv[++i]);
        } else if (std::strcmp(argv[i], "--nstep") == 0 && i + 1 < argc) {
            cfg.nstep = std::atoi(argv[++i]);
        } else if (std::strcmp(argv[i], "--freq") == 0 && i + 1 < argc) {
            cfg.fout = std::atoi(argv[++i]);
        } else if (std::strcmp(argv[i], "--theta") == 0 && i + 1 < argc) {
            cfg.theta = std::atof(argv[++i]);
        } else if (std::strcmp(argv[i], "--run") == 0 && i + 1 < argc) {
            cfg.run = argv[++i];
        } else if (std::strcmp(argv[i], "--init") == 0 && i + 1 < argc) {
            cfg.filename = argv[++i];
        } else {
            throw std::runtime_error(std::string("Unknown or incomplete argument: ") + argv[i]);
        }
    }

    return cfg;
}

int main( int argc, char *argv[] ) {

    config cfg = parse_args( argc, argv );

    char PATH[512];
    std::snprintf(PATH, sizeof(PATH), "%s/%s", INITPATH, cfg.filename);


    scalar size = cfg.size * PC; // total simulation size
    vec c = { -size/2, -size/2, size/2};
    Octree *bhtree = new Octree( c.x, c.y, c.z, size ); // initializing our tree

    // >>> new test suite, 1000 particles for a realistic cluster

    scalar *x; scalar *y; scalar *z;
    scalar *vx; scalar *vy; scalar *vz;
    scalar *m;

    int n = cfg.n;

    x   = (scalar *) malloc( sizeof(scalar) * n);
    y   = (scalar *) malloc( sizeof(scalar) * n);
    z   = (scalar *) malloc( sizeof(scalar) * n);
    vx  = (scalar *) malloc( sizeof(scalar) * n);
    vy  = (scalar *) malloc( sizeof(scalar) * n);
    vz  = (scalar *) malloc( sizeof(scalar) * n);
    m   = (scalar *) malloc( sizeof(scalar) * n);


    // actually reading in our data
    scalar *lines[7] = { x, y, z, vx, vy, vz, m };

    FILE *fp = fopen(PATH, "r");
    // std::cout << "reading file...\n"; 

    if (!fp) {
        perror("Error opening file");
        return 1;
    }

    // int newline = 0;

    for (int row = 0; row < 7; row++) {
        for (int col = 0; col < n; col++) {
            int ret = fscanf(fp, " %f ", &lines[row][col]);
            if (ret == EOF) {
                printf("Reached EOF at row %d col %d\n", row, col);
                return 1;
            }
        }
    }

    fclose(fp);
    // std::cout << "done reading file...\n"; 

    bhtree->build_tree(n, x, y, z, vx, vy, vz, m);
    // std::cout << "tree built...\n"; 
    scalar dt = cfg.dt * YR;
    scalar simtime = 0.0;
    scalar theta = cfg.theta;

    for ( int t = 0; t < cfg.nstep; t++) {

        bhtree->compute_forces( theta, dt);
        if (t % cfg.fout == 0)
            bhtree->save_step( t, simtime, theta, cfg.run, cfg.run);
        simtime += dt;
        
    }

    // >>> memory cleanup on aisle zero.
    free(x);
    free(y);
    free(z);
    free(vx);
    free(vy);
    free(vz);
    free(m);

    return 0;
        
}