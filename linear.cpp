/*
 * File:   linear.cpp
 * Author: Kevin Huang
 *
 * Created on 15 August 2020
 */

#include "linear.hpp"

LinearRegression::LinearRegression(MatrixContinuedFraction fraction, int starting_depth, data_store& train) {
    this->fraction = fraction;
    this->starting_depth = starting_depth;

    load_data(train);
    
    max_depth = fraction.depth;
    
    C = Eigen::VectorXd::Zero(max_depth + 1);

    compute_starting_residual();

    // Only use variables which have already been used in the fraction
    var_used = Eigen::VectorXd::Zero(fraction.vars + 1);
    for (size_t r = 0; r < fraction.base_terms; r++) {
        for (size_t c = 0; c < fraction.vars; c++) {
            if (fraction.coeffs[r][c] != 0) { //moh: need to use a proper floating point comparison function
                var_used(c) = true;
            }
        }
    }
    // constant is always used
    var_used(fraction.vars) = true;

    // If variable isn't used, then set the corresponding column
    // in X to 0
    for (size_t i = 0; i < (size_t)X.cols(); i++) {
        X.col(i) = X.col(i).array() * var_used(i);
    }
}

void LinearRegression::load_data(data_store& data) {
    /* 
     * load information into Eigen matrices for easy linear algebra
     */
    vector<vector<double>> dataX = data.input;
    vector<double> datay = data.expected;

    X = Eigen::MatrixXd::Zero(data.num_entry, data.num_var + 1);
    y = Eigen::VectorXd::Zero(data.num_entry);

    for (size_t i = 0; i < dataX.size(); i++) {
        for (size_t j = 0; j < dataX[0].size(); j++) {
            X(i, j) = dataX[i][j];
        }
        // constant term multiplies by 1
        X(i, dataX[0].size()) = 1;
        y(i) = datay[i];
    }
    yfit = y;

    //load fraction coefficients
    size_t starting_term = num_index(starting_depth);
    if (starting_depth == 0) {
        starting_term = 0;
    }
    coef = Eigen::MatrixXd::Zero(fraction.vars + 1, fraction.base_terms);
    for (size_t r = 0; r < fraction.base_terms && r < starting_term; r++) {
        for (size_t c = 0; c <= fraction.vars; c++) {
            if (!fraction.active[r][c]) {
                coef(c, r) = 0;
            } else {
                coef(c, r) = fraction.coeffs[r][c];
            }
        }
    }
}

void LinearRegression::compute_starting_residual() {
     /* 
     * Compute the starting value for yfit, so that we start fitting at starting_depth
     */
    for (int d = 0; d < starting_depth; d++) {
        Eigen::VectorXd residual = Eigen::VectorXd::Zero(yfit.rows());

        // Compute the residual of d-th linear regression
        if (d == 0) {
            residual = yfit - eval_term(0, X);
        } else {
            residual = yfit.array() * eval_term(num_index(d), X).array() - eval_term(den_index(d), X).array();
        }

        if (d == starting_depth - 1) {
            // Ensure residual is positive to avoid poles
            C(d) = abs(residual.minCoeff()) + 1;
            residual = residual.array() + C(d);
        }
        yfit = 1 / residual.array();
    }
}

Eigen::VectorXd LinearRegression::eval_term(int i, Eigen::MatrixXd& X_eval) {
    /* 
     * Evaluate the term at X_eval 
     */
    return X_eval * coef.col(i);
}

Eigen::VectorXd LinearRegression::eval(Eigen::MatrixXd& X_eval) {
    /* 
     * Evaluate the entire fraction at X_eval
     */
    if (max_depth == 0) {
        return eval_term(0, X_eval);
    }
    
    return eval_term(0, X_eval) + eval_rec(X_eval, 1);
}

Eigen::VectorXd LinearRegression::eval_rec(Eigen::MatrixXd& X_eval, int d) {
    /*
     * Recursive helper function for eval
     */
    if (d == max_depth) {
        return eval_term(num_index(d), X_eval).array() / eval_term(den_index(d), X_eval).array();
    }

    return eval_term(num_index(d), X_eval).array() / (eval_term(den_index(d), X_eval).array() + eval_rec(X_eval, d + 1).array());
}

void LinearRegression::fit() {
    for (int d = starting_depth; d <= max_depth; d++) {
        fit_depth(d);
    }
    update_fraction();
    log_debug();
}

void LinearRegression::fit_depth(int d) {
    /*
     * A helper method for fit(). This just does linear regression for a specific
	 * depth d, and then sets up yfit for the next depth.
     */
    if (d > 0) {
        coef(fraction.vars, num_index(d)) = 1;
        coef(fraction.vars, den_index(d - 1)) = coef(fraction.vars, den_index(d - 1)) - C(d - 1);
    }

    Eigen::VectorXd den_solution = linear_least_squares(X, yfit);

    // Get rid of coeff for unused variables
    den_solution = den_solution.array() * var_used.array();

    coef.col(den_index(d)) = den_solution;

    Eigen::VectorXd residual = yfit - eval_term(den_index(d), X);

    C(d) = abs(residual.minCoeff()) + 1;
    residual = residual.array() + C(d);

    yfit = 1 / residual.array();
}

void LinearRegression::update_fraction() {
    /*
     * Updates MatrixContinuedFraction with linear regression results. 
	 * Should be called only after fitting
     */
    for (int c = 0; c < coef.cols(); c++) {
        for (int r = 0; r < coef.rows(); r++) {
            //if (almost_equal(coef(r, c), g_ZERO, g_PREC))
            if (coef(r, c) == 0) {
                fraction.active[c][r] = false;
            } else {
                fraction.active[c][r] = true;
            }

            fraction.coeffs[c][r] = coef(r, c);
        }
    }
}

Eigen::VectorXd LinearRegression::linear_least_squares(Eigen::MatrixXd& X, Eigen::VectorXd& y) {
    /* 
     * Solve a linear least squares problem on X and y. Returns the weights
     */
    return X.colPivHouseholderQr().solve(y);
}

void LinearRegression::log_debug() {

    ofstream fout("debug.txt");

    Eigen::VectorXd y_pred = eval(X);
    for (int i = 0; i < X.rows(); i++) {
        fout << X(i, 0) << " " << y_pred(i) << " " << y(i) << endl;
    }
}