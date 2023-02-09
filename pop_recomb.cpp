
#include "pop.hpp"

template <typename Operator>
void population::variable_recomb(agent* a, agent* b, Operator op) {
    /*
     *
     */

    // cout << "---Start variable_recomb" << endl;
    auto move_val = [&, this] (double a, double b) {

        double coeff = a + g_rint(-1, 4) * (b - a) / 3;
        if( g_int_only ) 
            coeff = round(coeff);

        return coeff;
    };

    MatrixContinuedFraction& p1 = a->member[0];
    MatrixContinuedFraction& p2 = b->member[1];
    MatrixContinuedFraction& ch = a->member[1];

    // !!! With g_type == 1 CCL we are not considering the exponent IVs for recombination !!!
    // With g_type == 2 we do not want to recombine exponent IVs
    for (int i = 0; i < num_var; ++i) {
        ch.feature_active[i] = op(p1.feature_active[i], p2.feature_active[i]);
    }
    
    size_t max_dep = min(ch.depth, min(p1.depth, p2.depth));
    
    size_t start_term = 0;
    if( g_depth_lock != -1)
        start_term = 2*g_depth_lock+1;

    for (size_t i = start_term; i < 2 * max_dep + 1; ++i) {

        vector<double> acoeff(p1.coeffs[i], p1.coeffs[i]+(num_var + 1));
        vector<bool> afeature(p1.active[i], p1.active[i]+(num_var + 1));

       vector<double> bcoeff(p2.coeffs[i], p2.coeffs[i]+(num_var + 1));
        vector<bool> bfeature(p2.active[i], p2.active[i]+(num_var + 1));
       
       
        // find the shared variables
        vector<double> new_coeff(num_var + 1);
        for (int j = 0; j < num_var; ++j) {
            if (!ch.feature_active[j]) {;
                continue;
            }

            if (afeature[j] && bfeature[j]) {
                // six equally spaced points around the parents
                new_coeff[j] = move_val(acoeff[j], bcoeff[j]);
            } else if (afeature[j])
                new_coeff[j] = acoeff[j];
            else if (bfeature[j])
                new_coeff[j] = bcoeff[j];
        }

        new_coeff.back() = move_val(acoeff.back(), bcoeff.back());
        copy(new_coeff.begin(), new_coeff.end(), ch.coeffs[i] );

        // Unmask non-zero coeffs
        for (size_t j = 0; j < ch.vars; ++j) {

            if (abs(ch.coeffs[i][j]) < g_EPS) {
                ch.coeffs[i][j] = g_ZERO;
                ch.active[i][j] = false;
            } 
            else 
                ch.active[i][j] = true;
        }

    }

    for (size_t i = 2 * max_dep + 1; i < ch.depth; ++i) {
        for (int j = 0; j < num_var; ++j) {
            if (!ch.feature_active[j]) {
                ch.active[i][j] = false;
            }
        }
    }

}

void population::variable_intersect(agent* a, agent* b) {
    /*
     *
     */
    variable_recomb(a, b, bit_and<bool>());
}

void population::variable_union(agent* a, agent* b) {
    /*
     *
     */
    variable_recomb(a, b, bit_or<bool>());
}

void population::variable_symdiff(agent* a, agent* b) {
    /*
     *
     */
    variable_recomb(a, b, bit_xor<bool>());
}
