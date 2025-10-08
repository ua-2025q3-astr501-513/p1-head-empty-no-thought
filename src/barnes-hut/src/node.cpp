#include "node.h"

Node::Node( vec c, scalar s) {
    this->mass = 0.;
    this->dx = s;
    this->nchildren = 0;
    this->corner = c;
    this->com = {}; // workaround to get a zero vector, lazy :/
    
    this->parent = nullptr;
    this->particle = nullptr;
}

bool Node::is_internal( ) {
    if (this->particle != nullptr)
        return false;
    return true;
}

bool Node::contains( vec v ) {

    if ( ( v.x > this->corner.x && v.x <= this->corner.x + this->dx ) &&
            ( v.y > this->corner.y && v.y <= this->corner.y + this->dx ) &&
            ( v.z > this->corner.z && v.x <= this->corner.z + this->dx ) ) {
                return true;
            }
    return false;
}

void Node::insert( Body* b) {

    if (this->particle == nullptr ) {
        this->particle = b;
    } else {
        // do some recursion here into the child nodes...
        // update com
        // update mass
    }
}

void Node::update_mass( const Body* b) {
    // do some recursion here too? or maybe this is the only recursive area
}

void Node::get_force( Body *b, scalar theta) {
    // this is where the magic of barnes hut happens.
}