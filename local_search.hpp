/*
 * File:   local_search.hpp
 * Author: Haoyuan Sun
 *
 * Description
 * - Local search functions
 * 
 * Improvement List
 * -
 * 
 * 210713, AC, Cleanup
 * 180700, HS, Initial revision
 */

#ifndef __OPTIMIZE_HPP__
#define __OPTIMIZE_HPP__

using namespace std;

// Std Lib
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <thread>
#include <atomic>
#include <mutex>
#include <iostream>
#include <cassert>
#include <limits>
#include <stdlib.h>     // exit, EXIT_FAILURE
#include <omp.h>
// Local
#include "data.hpp"
#include "eval.hpp"
#include "cont_frac.hpp"
#include "globals.h"
#include "linear.hpp"
#include "cont_frac.hpp"
#include "eval.hpp"
#include "eval.hpp"
#include "helper.hpp"
#include "data.hpp"

class optimize {

	private:

		const data_store& train_data;
		MatrixContinuedFraction& frac;

		string objec; //type of objective function
		double delta; //delta for objective function


		using pii = pair<int, int>;
		int ndim;
		vector<pii> var_map;

		double eval_fit (vector<double>& vec, MatrixContinuedFraction& buf,
		const evaluator& e) const;
		double nelder_mead (MatrixContinuedFraction& frac, const evaluator& e) const;
		double runSerial();
		double runParallel();


	public:
		
		optimize(const data_store& td, MatrixContinuedFraction& f, string objective, double del);
		double run ();

};

#endif // __OPTIMIZE_HPP__
