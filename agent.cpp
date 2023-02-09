
#include "agent.hpp"

/**
 * Create an agent with agent_degree children and link to parent
 * 
 * @param agent_depth       depth of the agent tree
 * @param agent_degree      number of childen nodes 
 * 
 * @return  agent
 */
agent::agent(size_t agent_depth, size_t agent_degree) {

    depth = agent_depth;
    degree = agent_degree;

    if( g_debug )
        cout << "Agent() dep:" << agent_depth << " deg:" << agent_degree << endl;

    // Root node has no parent
    this->parent = nullptr;

    // If not a leaf node
    if (depth > 1) {

        // Create degree number of children
        children = new agent*[degree];

        // Assign the childrens parent
        for (size_t i = 0; i < degree; ++i) {
            children[i] = new agent(depth - 1, degree);
            children[i]->parent = this;
        }

    } 
    // If not a leaf node, childen are null
    else {
        children = nullptr;
    }
}

/**
 * Deconstruct agent and its children
 * 
 * @return  void
 */
agent::~agent() {

    if( g_debug )
        cout << "~Agent() dep:" << depth << " deg:" << degree << endl;

    if (depth > 1) {
        for (size_t i = 0; i < degree; ++i) delete children[i];
        delete[] children;
    }
}

/**
 * Swap the pocket/current from agent1 with the pocket/current from agent2
 * 
 * @param agent1        first agent
 * @param i             pocket or current of first agent
 * @param agent2        second agent
 * @param j             pocket or current of second agent
 * 
 * @return  void
 */
void agent::swap(agent* agent1, int i, agent* agent2, int j) {

    if( g_debug )
        cout << "swap(agent* agent1, int i, agent* agent2, int j) i:" << i << " j:" << j << endl;

    double tmp_fit = agent1->fitness[i];
    agent1->fitness[i] = agent2->fitness[j];
    agent2->fitness[j] = tmp_fit;

    MatrixContinuedFraction tmp_frac = agent1->member[i];
    agent1->member[i] = agent2->member[j];
    agent2->member[j] = tmp_frac;
}

/**
 * Update pocket and current if the fitness of current is greater than the pocket
 * 
 * @return  bool true if swap occurs
 */
bool agent::update_pocket() {

    if( g_debug )
        cout << "update_pocket()" << endl;

    if (fitness[0] > fitness[1]) {
        agent::swap(this, 0, this, 1);
        return true;
    }
    return false;
}

/**
 * Propogate less fit pockets down the tree by the path of the fittest child
 * 
 * @return  void
 */
void agent::movedown_pocket() {

    if( g_debug )
        cout << "movedown_pocket()" << endl;

    update_pocket();

    // Dont process leafs - they have no children
    if (depth <= 1) return;

    // Determine best child
    int best_child = 0;
    double best_fit = this->children[0]->fitness[0];
    for (size_t i = 1; i < degree; ++i) {
        if (children[i]->fitness[0] < best_fit) {
            best_child = i;
            best_fit = children[i]->fitness[0];
        }
    }

    // Swap pockets of the parent and child if the parent is not fitter than the child
    if (fitness[0] > best_fit) {
        swap(this, 0, children[best_child], 0);
        children[best_child]->movedown_pocket();
    }
}

/**
 * Propogate more fit childen up the parents if they are fitter than the parent
 * 
 * @return  void
 */
void agent::propagate_pocket() {

    if( g_debug )
        cout << "propagate_pocket()" << endl;

    // Dont process root node - no parent
    if (!this->parent) return;

    // Trade child and parent pocket if the child is fitter than the parent
    while (fitness[0] < this->parent->fitness[0]) {
        agent::swap(this, 0, this->parent, 0);
        this->parent->propagate_pocket();
        this->update_pocket();
    }
}

// Heplers

void agent::debug(vector<string> names) {
    /*
     * Debug an agent
     */

    cout << "Pocket fitness:" << member[0].fitness << "\t" << g_objec << " score:" << member[0].error << endl;
    member[0].write_cfr(cout, names, Numpy);
    member[0].write_term_table(cout, names, Numpy);
    cout << endl;
   
    cout << "Current fitness:" << member[1].fitness << "\t" << g_objec << " score:" << member[1].error << endl;
    member[1].write_cfr(cout, names, Numpy);
    member[1].write_term_table(cout, names, Numpy);
    cout << endl;
}