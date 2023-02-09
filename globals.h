/*
 * File:   globals.h
 * Author: Mohammad Haque
 *
 * Description
 * - Define global variables
 * 
 * Improvement List
 * -
 * 
 * 210713, AC, Cleanup
 * 190705, MH, Initial version 
 */

#ifndef GLOBALS_H
#define GLOBALS_H

#include <string>
#include "rng.hpp"

using namespace std;

extern const string VERSION;

extern char* g_trnFile;         // Training filepath
extern char* g_tstFile;         // Testing filepath
extern char* g_solnFile;        // Solution output filepath
extern char* g_logFile;         // Log filepath
extern string g_outPath;         // Location of output to created (.log, sol.json, -Train.csv, -Test.csv) for depth-lock and iterative execution
extern float g_delta;           // Fitness penalty multiplier to reduce number of used variables
extern string g_objec;          // Objective function text identifier
extern string g_numFrac;        // 
extern bool g_loocv;            // Undertake leave-on-out cross validation
extern int g_norm;              // Type of normalisation preprocessing behaviour 
extern size_t g_depth;          // Depth of graction
extern int g_numGen;            // Number of generations to run
extern uint_fast32_t g_seed;    // Seed identifier
extern float g_mu_rate;         // Mutation probability
extern bool g_run_serial;       // Serial or parallel execution
extern int g_percnt_data;       // Percentage of data to utilise in local search
extern int g_verbose_mode;      // Regulate amount of function output; 0 => Production; 1 => Exploration; 2 => Deep debug
extern int g_local_search;      // Number of generations after which local search is triggered
extern int g_reset_stuck;       // Number of generations of stagnating fitness to be stale
extern int g_nm_runs;           // Number of times the local search algorithm is run
extern int g_nm_max_iters;      // Number of iterations in the local search
extern int g_nm_reset_stag;     // Number of local search iterations of stagnating fitness to be stale
extern int g_PREC;              // Resolution of decimal calculation and output
extern double g_EPS;            // Value that we determine is small enough to be identical to another compared value
extern double g_ZERO;           // Value that we determine is small enough to be identical to zero
extern randreal g_rreal;        // Double randomisation object
extern randint g_rint;          // Integer reandomisation object
extern double EPS;              //
extern double MIN;              //
extern double g_window_size;    // Percentage of data selected at random for use in local search
extern bool g_weight_present;   // Flag determining if weight is in the datasets
extern char* g_pop_init;        // Type of initialisation
extern int g_pop_init_method;   // Method of initialisation when non-default is indicated
extern double g_linreg_chance;  // 
extern int g_type;              // Represenation type of continued fractions; 0 => traditional, 1 => continued log, 2 => power law
extern int g_base;              // Base used for the continued log type
extern bool g_cont;             //
extern bool g_origin;           // Force function to pass through the origin
extern bool g_aic_paramCount;   // AIC:: Supply the number of free parameters in each of the Meta-features
extern char* g_paramCountFile;   // ParamCount information filepath
extern int g_int_terms;         // Force terms to evaluate to integer terms
extern bool g_int_only;         // Only allow integer coefficients or constants
extern bool g_prev_model_as_var;  // To include previous model as a variable in current depth solutions
extern bool g_debug;

// !!! Temporary vars commited quickly - need to make this more robust
extern int g_depth_lock;                    // Lock terms up to and including this depth
extern int g_run_no;                        // Run when iterating the CFR with different depths
extern double g_fitness;                    // Global fitness when iterating many CFR executions
extern bool g_incremental_depth;            // Flag for iterative learning of the CFR
// extern vector<int> aic_param_count;         // AIC: parameter cout for feature/meta-feature
#endif /* GLOBALS_H */
