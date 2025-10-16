#include "node.h"
#include "body.h"
#include "util.h"
#include <iostream>

/** 
 * constructor, basic. all fields aside from position and size are set to 
 * zero or null pointers.
 * 
 * @param c the coordinates of the upper left corner of the node
 * @param s the physical size of the node's domain [m]
*/
Node::Node( vec c, scalar s) {
    this->mass = 0.;
    this->dx = s;
    this->nchildren = 0;
    this->corner = c;
    this->com = {0, 0, 0}; // workaround to get a zero vector, lazy :/
    
    this->parent = nullptr;
    this->particle = nullptr;
    
    for (int i = 0; i < 8; i++) {
        this->children[i] = nullptr;
    }
}

/** 
 * destructor, because i forgot to destroy the children (oops)
 * 
 * note that we DON'T destroy the particle (body) here >> 
 * that pointer is still needed by the tree to rebuild.
 */
Node::~Node( ) {
    this->parent = nullptr;
    this->particle = nullptr;
    for (int i = 0; i < 8; i++) {
        this->children[i] = nullptr;
    } 
}

/** 
 * checks if a node is internal to the tree.
 * 
 * note: an internal node is one that does have children.
 * 
 * @returns true if the node is internal (has children); false otherwise.
 */
bool Node::is_internal( ) {
    if (nchildren > 0)
        return true;
    return false;
}

/**
 * checks to see if a position vector is within the domain
 * of this node.
 * 
 * @returns true if v inside this node; false otherwise.
*/
bool Node::contains( vec v ) {

    if ( ( v.x > this->corner.x && v.x <= this->corner.x + this->dx ) &&
            ( v.y > this->corner.y && v.y <= this->corner.y + this->dx ) &&
            ( v.z > this->corner.z && v.z <= this->corner.z + this->dx ) ) {
                return true;
            }
    return false;
}

/**
 * recursively inserts a particle (body) into this node (or its children).
 * 
 * psuedocode taken from Thomas Trost's lecture slides [here](https://www.tp1.ruhr-uni-bochum.de/~grauer/lectures/compI_IIWS1819/pdfs/lec10.pdf)
 * 
 * @param b the body to be added to the node.
*/
void Node::insert( Body* b) {


    if ( particle == nullptr && !is_internal() ) {
        // creates a new particle if this node doesn't have one
        particle = b;
        return;
    } else if ( is_internal() ) {
        int q = get_quadrant(this->dx, this->corner, b->pos);
        if (children[q] == nullptr) {
            vec cnew = get_new_corner(q, this->corner, this->dx);
            children[q] = new Node( cnew, this->dx/2);
            children[q]->parent = this;
            nchildren++;
        }
        children[q]->insert( b );
        return;
    } 

    // inserting this node's particle into a new subdivision/node.
    int qold = get_quadrant( dx, corner, particle->pos);
    int qnew = get_quadrant( dx, corner, b->pos);

    if (children[qold] == nullptr) { // if there's not already a node there
        vec cnew = get_new_corner(qold, this->corner, this->dx);
        children[qold] = new Node( cnew, this->dx/2);
        children[qold]->parent = this;
        nchildren++;
    }

    children[qold]->insert(particle);
    particle = nullptr;

    if (children[qnew] == nullptr) { // if there's not already a node there
        vec cnew = get_new_corner(qnew, this->corner, this->dx);
        children[qnew] = new Node( cnew, this->dx/2);
        children[qnew]->parent = this;
        nchildren++;
    }

    children[qnew]->insert(b);

    update_mass( ); // not sure if this needs to happen here or not
    return;
    
}

/**
 * updates the total node mass (e.g. if we added a body or a child node.)
*/
void Node::update_mass( ) {
    
    // if this is a node with no children and one body:

    if (!is_internal()) {

        if (particle != nullptr) {
            mass = particle->mass;
            com = particle->pos;
            return;
        } 

        // if the node is empty for some reason...
        mass = 0;
        com = {0, 0, 0};
        return;
    }

    // otherwise, go through the children!
    scalar tmass = 0;
    vec tcom = {};

    for (int i = 0; i < 8; i++) {
        if (children[i] != nullptr) {
            children[i]->update_mass( );
            tmass += children[i]->mass;
            tcom += (children[i]->com * children[i]->mass); // weighted sum!
        }
    }

    if (tmass > 0) {
        mass = tmass;
        com = tcom * (1 / tmass); // final weighted sum
    } else { // if all children are empty.
        mass = 0;
        com = {0, 0, 0};
    }

    return;

}

/**
 * calculates the net force this node exerts on another body,
 * either using a direct calculation or a center of mass estimate
 * with the barnes-hut algorithm. recursive, to ensure we track through
 * the tree if needed.
 * 
 * psuedocode taken from Thomas Trost's lecture slides [here](https://www.tp1.ruhr-uni-bochum.de/~grauer/lectures/compI_IIWS1819/pdfs/lec10.pdf)
 * 
 * @param b the body we're calculating force on
 * @param theta threshold criteria for barnes-hut, ratio of node width to distance to center of mass.
 * 
 * @returns a force vector.
 * 
*/
vec Node::get_force( Body* b, scalar theta) {
    // this is where the magic of barnes hut happens.
    vec force = {};
    scalar F = 0;

    if (mass == 0) { return force; }

    vec rdiff = com - b->pos;
    scalar r = rdiff.norm();

    // if this node has one particle
    if (!is_internal() && particle != nullptr && particle != b) {
        // calcuate the force 
        F = G * particle->mass * b->mass / pow(r, 3); // magnitude of force
        b->acc += rdiff * F * (1/b->mass); // updating particle acceleration from force calculation.
        // i should probably update velocity and position when we handle timesteps.
        return rdiff * F; // force vector!
    }

    // barnes-hut approximation for this node
    if ( dx / distance(com, b->pos) < theta ) {
        F = G * mass * b->mass / pow(r, 3); // calcuates force using total node mass
        return rdiff * F;
    }

    // otherwise, we look at the child nodes (recursion)
    for (int i = 0; i < 8; i++) {

        if ( children[i] != nullptr ) {
            force += children[i]->get_force(b, theta);
        }
    }

    return force;

}