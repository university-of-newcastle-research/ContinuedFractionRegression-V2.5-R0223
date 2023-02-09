/*
 * File:   linear.h
 * Author: Kevin Huange
 *
 * Description
 * - This class is a wrapper for MatrixContinuedFraction that performs a linear regression
 * 	 as an initialization method for the coefficients. 
 *   A LinearRegression object should be instantiated whenever you wish to use linear
 *   regression to fit certain depths for an existing MatrixContinuedFraction object.
 * 
 * Improvement List
 * -
 * 
 * 210713, AC, Cleanup
 * 200815, KH, Initial version 
 */

#ifndef __LINEAR_HPP
#define __LINEAR_HPP

// Std Lib
#include <Eigen/Dense>
#include <fstream>
// Local
#include "data.hpp"
#include "cont_frac.hpp"

class LinearRegression {

	private:
		
		int starting_depth;			// depth to start linear regression from
		int max_depth;				// depth of the fraction

		// X is a n x m data matrix, y is a size n target vector, where
		// n is the number of data points, and m is the number of features, plus 1
		// for the constant
		Eigen::MatrixXd X;
		Eigen::VectorXd y;

		Eigen::VectorXd var_used;		// boolean vector of whether a variable is included or not
		Eigen::VectorXd yfit;			// yfit is the current target to fit for a certain depth for linear regression
		Eigen::VectorXd C;				// C contains the constants that are added to make sure that spurious poles do not occur

		// corresponds to MatrixContinuedFraction::coeffs, but transposed
		// for better locality since Eigen matrices are stored column by column
		// by default
		Eigen::MatrixXd coef;

		// Helper functions to convert from term indices to depth indices
		inline int term_to_depth(int term) { return (term + 1) / 2; }
		inline int num_index(int depth) { return 2 * depth - 1; }
		inline int den_index(int depth) { return 2 * depth; }

		void load_data(data_store& data);
		void compute_starting_residual();
		void fit_depth(int d);
		void update_fraction();
		
		Eigen::VectorXd eval_term(int i, Eigen::MatrixXd& X_eval);
		Eigen::VectorXd eval(Eigen::MatrixXd& X_eval);
		Eigen::VectorXd eval_rec(Eigen::MatrixXd& X_eval, int d);
		Eigen::VectorXd linear_least_squares(Eigen::MatrixXd& X, Eigen::VectorXd& y);

	public:

		MatrixContinuedFraction fraction;		// The fraction to perform linear regression on

		LinearRegression(MatrixContinuedFraction fraction, int starting_depth, data_store& train);
			
		void fit();
		void log_debug();

};

#endif // __LINEAR_HPP
