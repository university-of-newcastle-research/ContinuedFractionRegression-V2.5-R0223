
/*
 * File:    cont_frac.cpp
 * Authors: mnhaque
 * 
 * Description
 * - Matrix implementation of the CFR for increased efficiency
 * 
 * Improvements
 * -
 * 
 * 210930, AC, Cleanup and documetation
 * 210714, AC, Create rand_var_at such that all coefficient changes occur in singular func
 * 210713, AC, Utilise CCL fraction structure for a power-law representation,
 *             General cleanup and documentation
 * 210000, MH, Nil constant in first term of CFR to force origin intercept
 * 200526, MH, Created first version 11:07
 */

#ifndef MATRIXCONTINUEDFRACTION_HPP
#define MATRIXCONTINUEDFRACTION_HPP

// Std Lib
#include <random>
#include <string>
#include <chrono> 
#include <ctime>
#include <sstream>
#include <iomanip>
#include <cerrno>
#include <unordered_map>
#include <map>
// Packages
#include "nlohmann/json.hpp"

// Local code
#include "globals.h"
#include "helper.hpp"

using namespace std;
using namespace chrono;
using json = nlohmann::json;


/**
 * Fraction outputs
 * 
 * @return          void
 */
enum ptype {
    Latex,
    Excel,
    Numpy
};

class MatrixContinuedFraction {
    
    /*
     * Represent a continued fractions in a matrix form
     *
     * The original reference to 'term' is from the CFR paper (or CFs) which has been
     * denoted as g_i(x) and h_i(x) in the form of f(x) = g_0(x)+h_0(x)/(g_1(x)+h_2(x)/(g_2(x)+...))
     * where each g and h is a function of the IVs in the form g_i(x) = (c1x1+c2x2+...cNxN+constant)
     * 
     * This concept has become clouded when our CF representation has evolved to have multiple functions
     * of the IVs for each term. Consider the continued logarithms form 
     *  g_i(x) = (c1x1+c2x2+...cNxN+constant)*base^(c1x1+c2x2+...cNxN+constant) or the power law representation
     *  g_i(x) = (c1x1+c2x2+...cNxN+constant)^(c1x1+c2x2+...cNxN+constant) 
     * Complexity will increase with the user of more advanced forms such as fourier
     * 
     * For this reason we refer to "base term" as the inner terms and have a baseTermsPerTerm given the specific 
     * representation e.g. 1 for normal, 2 for continued logarithm and power law
     */

    // We create a matrix M x N matrix to map our continued fractions coefficients.
    // currently hardcoded for a maximum depth of 20 when utilising one base term per
    // continued fractions term. When we utilise the continued log or power law version
    // we reduce the depth that we can achieve with the same static matrix.
    // !!! we should explore the use of the eigen library in the matrix represenation
    
    private:

        void    init_terms();                       // Common setup for all constructors; number of base terms, depth etc.
        void    randomise(double constant = 0);     // Initial selection of coeffiencts and constants
        void    force_origin();                     // Hide the constant terms to force passage through origin

        size_t  term_depth(int);                    // Determine depth for base term
        size_t  term_no(int);                       // Determine term number for base term

        double  eval_term(const vector<double>& vals, int term);        // Calculate a specific term

        void    write_base_term(ostream& , vector<string> names, int term, ptype, bool print_rational=false);      // Print CFR base term
        void    write_term(ostream&, vector<string> names, int term, ptype, bool print_rational=false);            // Print CFR full term


    public:

        static const int M = 41;        // Hard-coded columns representing variables/constants
        static const int N = 1002;      // Hard-coded rows representing each term

        double  coeffs[M][N];           // Coefficients for each independent variable and the constant term
        bool    active[M][N];           // Local mask for each coeff or constant term

        
        //Depei's Iterative Depth and MWRAE objective function uses these lists
        unordered_map<int, double> residual_weight;
        unordered_map<int, double> temp_weight;

        double  penalty;                // Penalty value for the fraction due to the number of variables used
        double  error;                  // Error of the fraction given the objective function (g_objec)
        double  fitness;                // Fraction fitness as a combination of error and penalty
        int     freeParams;             // number of free parameters in the model, the value of K in AIC
        // bool    isRationalised;         // Flag to track the conversion of Continued Fraction into Rational format

        size_t  depth;                  // Number of convergents for the continued fractions
                                        // this is 1 for the first full term and +1 every two full terms there after e.g.
                                        //   f(x) = g_0(x)                              => Depth 0
                                        //   f(x) = g_0(x)+h_0(x)/g1(x)                 => Depth 1
                                        //   f(x) = g_0(x)+h_0(x)/(g1(x)+h1(x)/g2(x))   => Depth 2
                                        //   etc.

        size_t  vars;                   // Number of variables per base term where item at maximum position vars is the constant
        size_t  base_terms;             // Number of base terms 
        size_t  base_terms_per_term;

        vector<bool> feature_active;    // Global mask for each independent varaible down the term

        // Constructors

        MatrixContinuedFraction(size_t frac_depth = g_depth, size_t frac_vars = 0, double frac_alpha = 0);
        MatrixContinuedFraction(string filename, bool rand_extra_depth = true);         
        MatrixContinuedFraction(const MatrixContinuedFraction &cp);      
        ~MatrixContinuedFraction();
    
        int     used_features();

        void    rand_var_int(size_t, size_t, int min = 1, int max = 3 );                
        void    rand_var_real(size_t, size_t, double min = 0, double max = 1);          
        
        int     sum_param_count(vector<int> paramCounter);

        double  eval(const vector<double>& vals, int to_depth = g_depth);

        bool    import_cfr(string solfile, bool rand_extra_depth = true);                           
        bool    export_cfr(string, int, vector<string>);                   
               
        // void    Ratonalise();
        // Debug

        void    write_cfr(ostream&, vector<string>, ptype, bool isSummary = false, bool print_rational = false);    
        void    write_term_table(ostream&, vector<string>, ptype);       
};

#endif /* MATRIXCONTINUEDFRACTION_HPP */

