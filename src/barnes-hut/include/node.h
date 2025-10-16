#ifndef NODE_H
#define NODE_H

#include "util.h"
#include "body.h"

class Node {

    public:
        scalar mass; /** total mass in node [kg] */
        scalar dx; /** node size (physical) [m] */
        int nchildren; /** total number of child nodes */
        vec corner; /** coordinates of upper left corner */
        vec com; /** position of center of mass */

        Node* parent; /** parent node in the tree */
        Node* children[8]; /** list of pointers to child nodes */
        Body* particle; /** particle/body contained in node */

        Node( vec c, scalar s);
        ~Node();

        bool is_internal( );
        bool contains( vec v );
        void insert( Body* b);

        void update_mass( ) ;
        vec get_force( Body* b, scalar theta );

};

#endif