/*
 * File:   pop.hpp
 * Author: Haoyuan Sun
 *
 * Description
 * - Represent an agents population
 * 
 * Improvement List
 * -
 * 
 * 210713, AC, Cleanup
 * 190705, HS, Initial version 
 */

#ifndef __POPULATION_HPP
#define __POPULATION_HPP

using namespace std;

// Std Lib
#include <vector>
#include <fstream>
#include <cassert>
#include <iostream>
#include <deque>
#include <chrono>
#include <utility>
#include <iomanip>
#include <cmath>
#include <algorithm>
// Local
#include "agent.hpp"
#include "eval.hpp"
#include "data.hpp"
#include "globals.h"
#include "local_search.hpp"
#include "debug.hpp"
#include "linear.hpp"
#include "cont_frac.hpp"

class population {

    private:

        int num_var;
        data_store train_data, test_data; //const
        evaluator train_e, test_e;
        vector<string> names;       //names of the variables (first row of data)
        string objec;               //type of objective function
        double delta;               //delta for objective function

        string solfile;

        agent* root;

        // pop_init.cpp
        void initialize(agent* current);
        void initialize(agent* current, double alfa_0);
        void initialize(agent* current, string sol_file);

        //	void diversify(agent* a);
        //	void simplify(agent* a);

        // pop_evolve.cpp
        void eval_fit(agent* a, int i);
        void recombine(agent* a, agent* b);
        void mutate(agent* a);
        void local_search(agent* a);

        // pop_recomb.cpp
        template <typename Operator>
        void variable_recomb(agent* a, agent* b, Operator op);
        void variable_intersect(agent* a, agent* b);
        void variable_union(agent* a, agent* b);
        void variable_symdiff(agent* a, agent* b);

        // pop_mutate.cpp
        void feature_toggle(MatrixContinuedFraction& frac);
        void feature_mod(MatrixContinuedFraction& frac);

        // Initialisations
        void initialize_linear_root(agent* current, string sol_file="");
        void initialize_linear(agent* current);
        void refresh(agent* current);
        
        void breed(agent* leader);
        bool check(agent* a);
        string print_agent_vars(agent *a, int sol);       
        void print_population();
        void print_population_coeff();
        inline void print_score_header();

    public:

        population(data_store train);
        population(data_store train, data_store test);
        population(data_store train, data_store test, string soln_file);


        ~population();

        void run(MatrixContinuedFraction& sol);
        pair <double, double> run(int);
        
};

#endif // __POPULATION_HPP
