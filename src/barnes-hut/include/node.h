#ifndef NODE_H
#define NODE_H

#include "util.h"
#include "body.h"

class Node {

    public:
        scalar mass;
        scalar dx;
        int nchildren;
        vec corner;
        vec com;

        Node* parent;
        Node* children[8]; // can i do this?
        Body* particle; // single particle, if this is a leaf in the tree

        Node( vec c, scalar s);

        bool is_internal( );
        bool contains( vec v );
        void insert( Body* b);

        void update_mass( ) ;
        vec get_force( Body* b, scalar theta );

};

#endif