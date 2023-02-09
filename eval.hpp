
/*
 * File:    eval.hpp
 * Authors: Haoyuan Sun
 * 
 * Description
 * - 
 * 
 * Improvements
 * -
 *
 * 210713, AC, cleanup
 * 180000, HS, initial revision
 * 
 */

#ifndef __EVAL_HPP
#define __EVAL_HPP

// Std Lib
#include <unordered_set>
#include <string>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <limits>
#include <vector>       
#include <algorithm>    // Sort
#include <iomanip>
// Local
#include "globals.h"
#include "data.hpp"
#include "cont_frac.hpp"
#include "objective.hpp"

using namespace std;

class evaluator {
    
    private:

        const data_store& data;         // const
        unordered_set<int> selection;

    public:
        
        evaluator(const data_store& td);
        evaluator(const data_store& td, int num);
        evaluator(const data_store& td, double window);

        data_store getData();
        double eval_fit(MatrixContinuedFraction& frac) const;
        double eval_fit_full(MatrixContinuedFraction& frac) const;

        // moh: added functions to compute additional metrics
        double compute_MSE(MatrixContinuedFraction & frac) const;
        double compute_ChiSquared(MatrixContinuedFraction & frac) const;
        double compute_ChiSquaredPct(MatrixContinuedFraction& frac) const;
        double compute_MSET(MatrixContinuedFraction & frac, const data_store& td) const;
        double compute_RMSE(MatrixContinuedFraction & frac) const;
        double compute_NMSE(MatrixContinuedFraction & frac) const;
        double compute_PEARSON(MatrixContinuedFraction& frac) const;
        void compute_Scores(MatrixContinuedFraction& frac) const;
        void print_val(const char* filename, MatrixContinuedFraction& frac) const;
        void print_val(const char* filename, MatrixContinuedFraction& frac, const data_store& td) const;
        void print_val_nd(const char* filename, MatrixContinuedFraction& frac, const data_store& td) const;

};

#endif // __EVAL_HPP
