/*
 * File:   objective.hpp
 * Author: Haoyuan Sun, Mohammed Haque
 *
 * Description
 * - Define different objective functions for use
 * 
 * Improvement List
 * -
 * 
 * 210713, AC, Cleanup
 * 200819, AC, Migrate original work from eval.cpp developed by HS & MH 
 * 190000, MH, Initial work
 * 180800, HS, Initial revision
 */

#ifndef OBJECTIVE_HPP
#define OBJECTIVE_HPP

using namespace std;

// Std Lib
#include <unordered_set>
#include <cmath> /* log */
// Local
#include "data.hpp"
#include "globals.h"
#include "cont_frac.hpp"
#include "debug.hpp"

class ObjectiveFuncs {
private:

    // Define the series of objective functions used by the memetic algorithm
    ///   The following are common inputs to the evaulation
    ///   frac: the current continued fraction used to evaluate
    ///     data: the entire input dataset
    ///     subset: the elements within data to consider, which may be the entire set
    ///     multi_constant: ???
    ///

    static double mse(MatrixContinuedFraction& frac,
            const data_store& data,
            const unordered_set<int>& subset,
            const double multi_const);

    static double mwe(MatrixContinuedFraction& frac,
            const data_store& data,
            const unordered_set<int>& subset,
            const double multi_const);
    
    static double uwae(MatrixContinuedFraction& frac,
            const data_store& data,
            const unordered_set<int>& subset,
            const double multi_const);
    
    static double abserr(MatrixContinuedFraction& frac,
            const data_store& data,
            const unordered_set<int>& subset,
            const double multi_const);
    
    static double mase(MatrixContinuedFraction& frac,
            const data_store& data,
            const unordered_set<int>& subset);

    static double medae(MatrixContinuedFraction& frac,
            const data_store& data,
            const unordered_set<int>& subset);

    static double mederr(MatrixContinuedFraction& frac,
            const data_store& data,
            const unordered_set<int>& subset);

    static double emc(MatrixContinuedFraction& frac,
            const data_store& data,
            const unordered_set<int>& subset);

    static double theil(MatrixContinuedFraction& frac,
            const data_store& data,
            const unordered_set<int>& subset);

    static double pearson(MatrixContinuedFraction& frac,
            const data_store& data,
            const unordered_set<int>& subset,
            const double multi_const);

    static double nmse(MatrixContinuedFraction& frac,
            const data_store& data,
            const unordered_set<int>& subset,
            const double multi_const);

    static double ub(MatrixContinuedFraction& frac,
            const data_store& data,
            const unordered_set<int>& subset,
            const double multi_const);

    static double lb(MatrixContinuedFraction& frac,
            const data_store& data,
            const unordered_set<int>& subset,
            const double multi_const);

    static double cw(MatrixContinuedFraction& frac,
            const data_store& data,
            const unordered_set<int>& subset,
            const double multi_const);
    
    static double cwy(MatrixContinuedFraction& frac,
            const data_store& data,
            const unordered_set<int>& subset,
            const double multi_const);

    static double aic(MatrixContinuedFraction& frac,
            const data_store& data,
            const unordered_set<int>& subset);
    static double mre(MatrixContinuedFraction& frac,
            const data_store& data,
            const unordered_set<int>& subset,
            const double multi_const);

    static double mse_for_id(MatrixContinuedFraction& frac,
            const data_store& data,
            const unordered_set<int>& subset,
            const double multi_const);

    static double mwrae(MatrixContinuedFraction& frac,
            const data_store& data,
            const unordered_set<int>& subset,
            const double multi_const); 
public:

    // Define the evaluation function

    static double evaluate(MatrixContinuedFraction& frac,
            const data_store& dateset,
            const unordered_set<int>& subset,
            string objectiveFunc);

};

#endif