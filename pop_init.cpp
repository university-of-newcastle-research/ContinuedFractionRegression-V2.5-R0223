#include "globals.h"
#include "pop.hpp"

void population::initialize(agent *current) {
    /*
     *
     */
    if( g_debug )
        cout << "initialize(agent *current)" << endl;

    current->member[0] = MatrixContinuedFraction(g_depth, num_var); 
    current->member[1] = MatrixContinuedFraction(g_depth, num_var);

    eval_fit(current, 0);
    eval_fit(current, 1);

    // Show the Fitness after evaluating the solutions in the agent
    if( g_verbose_mode > 1 )    current->debug(names);

    //    cout << current->fitness[0] << "\t" << current->fitness[1] << endl;
    if (current->depth > 1) {
        for (size_t i = 0; i < current->degree; ++i)
            initialize(current->children[i]);
    }

    current->movedown_pocket();
}

void population::initialize(agent* current, double alfa_0){
    /*
     * initialise population with the alfa_0 (avg of the targets)
     *      as the value of const in all terms of CFraction.
     */

    if( g_debug )
        cout << "initialize(agent* current, double alfa_0)" << endl;

    current->member[0] = MatrixContinuedFraction(g_depth, num_var, alfa_0); //fraction(num_var);
    current->member[1] = MatrixContinuedFraction(g_depth, num_var, alfa_0); //fraction(num_var);

    eval_fit(current, 0);
    eval_fit(current, 1);

    if( g_verbose_mode > 1 )    current->debug(names);

    if (current->depth > 2) {
      for (size_t i = 0; i < current->degree; ++i)
          initialize(current->children[i], alfa_0);
    }
    else if (current->depth > 1) {
        for (size_t i = 0; i < current->degree; ++i)
            initialize(current->children[i]);
    }

    current->movedown_pocket();

}

void population::initialize(agent* current, string sol_file){
    /*
     * initialise population with the alfa_0 (avg of the targets)
     *      as the value of const in all terms of CFraction.
     */

    solfile = sol_file;
    if( g_debug )
        cout << "initialize(agent* current, const char* sol_file)" << endl;

    // Incremental Depth :: initialise pocket of root as the exact previous solution (extra terms as 0/1)
    if(g_incremental_depth==true)
        current->member[0] = MatrixContinuedFraction(sol_file, false); 
    else
        current->member[0] = MatrixContinuedFraction(sol_file); 
    // initialise the current as normal (extra terms with random)
    current->member[1] = MatrixContinuedFraction(sol_file);

    eval_fit(current, 0);
    eval_fit(current, 1);

    if( g_verbose_mode > 1 )    current->debug(names);

    if (current->depth > 1) {
         for (size_t i = 0; i < current->degree; ++i) {
            if( g_depth_lock == -1) { // do not include previous solution if depth is locked
                initialize(current->children[i]);
            }else{                  
                //seed with previous solution in Iterative-Model inclusion, but without Depth Lock
                initialize(current->children[i], sol_file);
            }
         }  
     }

    current->movedown_pocket();
}

void population::initialize_linear_root(agent* current, string sol_file) {
    /*
     * Initialize the population with linear regression. If provided, sol_file provides the starting
     * solution, and linear regression is used to initialize the rest of the fraction. Should be
     * only for the root.
     */

    if( g_debug )
        cout << "initialize_linear_root(agent* current, string sol_file)" << endl;

    if (sol_file != "") {
        //read num_term and num_var from file
        ifstream fin(sol_file);
        if (!fin.is_open()) return;

        size_t num_terms;
        size_t num_var;
        fin >> num_terms >> num_var;
        size_t frac_depth = (num_terms - 1) / 2;

        MatrixContinuedFraction frac(sol_file);

        // if the solution file already covers the whole fraction, no
        // deeper depths for linear regression to initialize
        if (frac_depth >= frac.depth) {
            return initialize(current, sol_file);
        }

        LinearRegression linreg(frac, frac_depth + 1, train_data);

        linreg.fit();

        current->member[0] = frac;
        current->member[1] = linreg.fraction;
    }
    else {
        
        size_t frac_depth = 0;

        MatrixContinuedFraction frac(g_depth, num_var);

        LinearRegression linreg(frac, frac_depth, train_data);

        linreg.fit();

        current->member[0] = frac;
        current->member[1] = linreg.fraction;
    }

    if( g_verbose_mode > 1 )    current->debug(names);

    eval_fit(current, 0);
    eval_fit(current, 1);

    if (current->depth > 1) {
        for (size_t i = 0; i < current->degree; ++i)
            initialize_linear(current->children[i]);
    }

    current->movedown_pocket();
}

void population::initialize_linear(agent* current) {
    /* 
     * Initialize the current agent with linear regression. Meant for non-root agents.
     * Only one of the pocket/current is initialized with linear regression. The other is initialized
     * as normal.
     * First, initializes the fraction as normal. Then, randomly selects a depth
     * from 0 - fraction depth, and uses linear regression to initialize from that depth
     * onwards.
     */
     
    if( g_debug )
        cout << "initialize_linear(agent* current)" << endl;

    current->member[0] = MatrixContinuedFraction (g_depth, num_var);

    MatrixContinuedFraction frac(g_depth, num_var);

    int frac_depth = g_rint(0, g_depth);
    LinearRegression linreg(frac, frac_depth, train_data);

    linreg.fit();
    current->member[1] = linreg.fraction;

    if( g_verbose_mode > 1 )    current->debug(names);

    eval_fit(current, 0);
    eval_fit(current, 1);

    if (current->depth > 1) {
        for (size_t i = 0; i < current->degree; ++i) {
            initialize_linear(current->children[i]);
        }
    }

    current->movedown_pocket();
}

void population::refresh(agent *current) {
    /* 
     *  recalculate the fitness score for new objec while switching from Pearson to MSE
     */

    double new_fit_p = current->member[0].error * current->member[0].penalty;
    current->member[0].fitness = new_fit_p;

    double new_fit_c = current->member[1].error * current->member[1].penalty;
    current->member[1].fitness = new_fit_c;

    current->fitness[0] = new_fit_p;
    current->fitness[1] = new_fit_c;
    //current.update_pocket();    // swap pocket & current if relevant

    if (current->depth > 1) {
        for (size_t i = 0; i < current->degree; ++i)
            refresh(current->children[i]);
    }
    current->movedown_pocket();
}
