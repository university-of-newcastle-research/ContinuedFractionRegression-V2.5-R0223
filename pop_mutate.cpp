
#include "pop.hpp"

void population::feature_toggle(MatrixContinuedFraction& frac) {
    /*
     *  Fewer feature -> high chance of turning on a feature
     */
    
    // Select IV at random
    int iv = g_rint(0, frac.vars - 1);

    // For all terms
    for (size_t term = 0; term < frac.base_terms; term++) {

        // If we have the constant exponent CFR and we are in the exponent
        // term then do not toggle
        if( g_type == 3 && ((term+1) % (int)frac.base_terms_per_term) == 0)
            continue;

        // If variable is unmasked
        if (frac.active[term][iv]) {

            // Turn varaible off
            frac.coeffs[term][iv] = 0;
            frac.active[term][iv] = false;
        }

        else if (!frac.active[term][iv]) {

            // Randomise coefficient with 50% chance of being set to 0
            if( g_rreal() < 0.5 )
                frac.rand_var_int(term, iv);
            else {
                // Turn varaible on but set to zero
                frac.coeffs[term][iv] = 0;
                frac.active[term][iv] = true;
            }
        }
    }

    // Toggle IV on/off flag
    frac.feature_active[iv] = !frac.feature_active[iv];

}

void population::feature_mod(MatrixContinuedFraction& frac) {
    /*
     * Randomly turn on/off a variable from an existing feature
     */

    if( g_debug )
        cout << "feature_mod(MatrixContinuedFraction& frac)" << endl;

    // Determine features that are on, take no action if all are off
    // !! Uncertain of why we return if there are no features on, don't we still want
    // !! to turn one on?
    vector<int> arr;
    for (size_t i = 0; i < frac.vars; ++i) {
        if (frac.feature_active[i])
            arr.push_back(i);
    }
    if (arr.empty()) return;

    // Select random feature that is ON 
    int iv = g_rint(0, arr.size() - 1);

    size_t start_term = 0;
    if( g_depth_lock != -1) 
        start_term = 2*g_depth_lock+1;

    // Select random term from the fraction
    int term = g_rint(start_term, frac.base_terms - 1);  

    // If we have the constant exponent CFR and we are in the exponent
    // term then do not toggle
    if( g_type == 3 && ((term+1) % (int)frac.base_terms_per_term) == 0)
        return;

    // When the variable is on
    if (frac.active[term][iv]) {

        // Turn it off
        frac.coeffs[term][iv] = 0;
        frac.active[term][iv] = false;

    } else {

        // Randomise coefficient with 50% chance of being set to 0
        if( g_rreal() < 0.5 )
            frac.rand_var_int(term, iv);
        else {
            // Turn varaible on but set to zero
            frac.coeffs[term][iv] = 0;
            frac.active[term][iv] = true;
        }
    }
}
