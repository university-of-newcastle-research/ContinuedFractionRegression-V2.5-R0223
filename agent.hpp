
/*
 * File:    agent.hpp
 * Authors: Haoyuan Sun
 * 
 * Description
 * - Represent agents of the tree-based population
 * 
 * Improvements
 * -
 *
 * 210713, AC, cleanup
 * 180700, HS, initial revision
 * 
 */

#ifndef __AGENT_HPP
#define __AGENT_HPP

// Std lib
#include <stdlib.h>
#include <iostream>
// Local
#include "globals.h"
#include "cont_frac.hpp"

class agent {

    private:

        int         depth;    // The depth below this node (inclusive), leaves has depth 1
                              // This is really a 1-based inverse depth and not best 

        size_t      degree;   // The degree of the tree structure, default is 3

        double                      fitness[2];     // Fitness for pocket and current (0 and 1 respectively)
        MatrixContinuedFraction     member[2];      // Fraction for pocket and current (0 and 1 respectively)

        agent* parent;              // Parent agent
        agent** children;           // Child agent

        static void swap(agent* a, int i, agent* b, int j);

        bool    update_pocket();
        void    movedown_pocket();
        void    propagate_pocket();

    public:

        agent(size_t agent_depth, size_t agent_degree = 3);
        ~agent();

        void    update();

        void    debug(vector<string>);

        friend class population;
};

#endif // __AGENT_HPP
