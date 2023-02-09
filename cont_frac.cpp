
#include "cont_frac.hpp"
#include <stdexcept>

using namespace std;

/**
 * Force the fraction through the origin
 * Required for physical applications where this property is known
 * 
 * @return          void
 */
void MatrixContinuedFraction::force_origin() {

    // This should be achieved by system of equations or another formal process 
    // however we crudely solve the issue by zeroing the appropriate constants such that
    // the when the ivs evaluate at 0 the fraction collapses to 0

    // In all CFR cases we remove the constant from g0(x)
    active[0][vars] = false;
    coeffs[0][vars] = 0.0;

    // We may need to consider other cases for different fractions
}

/**
 * Randomise the continued fraction terms
 * 
 * @param constant  Variable to initialise the term constants to
 *                  Defaults to 0 which indicates random
 *
 * @return          void
 */
void MatrixContinuedFraction::randomise(double constant) {

    // Randomise the global feature mask for the fraction with equal probability
    for (size_t iv = 0; iv < vars; ++iv) {

        if( g_rreal() < 0.5 )
            feature_active[iv] = true;
        else 
            feature_active[iv] = false;        
    
    }
  
    // Randomise coefficients for each variable of each term, including the constant (iv == vars)
    for (size_t term = 0; term < base_terms; term++) { 

        for (size_t var = 0; var <= vars; var++) {

            // Assume the mask and variable are off
            active[term][var] = false;
            coeffs[term][var] = 0;

            // Leave iv off if it is off globally and not a constant
            if(!feature_active[var] && var != vars)
                continue;

            // Set to forced constant when provided
            if( var == vars &&           // If the constant and
                constant != 0               // Constant value provided
            ) {
                coeffs[term][var] = constant;
                active[term][var] = true;
                feature_active[var] = true; //turn the Constant ON
                continue;
            }

            // Leave iv off if constant exponent CFR and we are considering the ivs in the exponent
            if( g_type == 3 &&                                  // If we are in constant exponent CFR and
                ((term+1) % (int)base_terms_per_term) == 0 &&   // We are looking at the exponent and
                var != vars                                  // We are considering the independent variables
            )
                continue;

            // Force constant to start at 1 for constant exponent CFR 
            if( g_type == 3 &&                                  // If we are in constant exponent CFR and
                ((term+1) % (int)base_terms_per_term) == 0 &&   // We are looking at the exponent and
                var != vars                                  // We are considering the independent variables
            ) {
                coeffs[term][var] = 1;
                active[term][var] = true;
            }
                
            // Otherwise initialise as integer
            rand_var_int(term, var);
        }

    }
}

/**
 * Initalise term settings based on the fraction depth and the type of fraction
 * 
 * @return          void
 */
void MatrixContinuedFraction::init_terms() {

    // For depth = 0 there is only one term g0(x) then for
    // subsequent depths we have two terms per depth
    // i.e. for depth 1 h1(x)/(g1(x)/...)
    //      for depth 2 h1(x)/(g1(x)/(h2(x)/g2(x))
    //      and so on
    base_terms = 2 * depth + 1;

    // List of features or indepedent variables x1,...,xN based 
    // from the dataset that is being represented by this CFR
    // this is essentially the column header of the matrix and there
    // exist a row of coefficients for each term.
    // 
    // We +1 for the constant term added to the end of the CFR
    feature_active = vector<bool>(vars + 1, false);
    feature_active[vars] = true; //turn the Constant ON
    
    // In the event of the continued logarithm form (g_type == 1)
    // or the power law form (g_type == 2) we have more than one place
    // in each term that is a function of the independent variables.
    // For instance gi(x) = fi(x1,..xn,c)*b^(fi(x1,...,xn,c))
    // In both cases, we have double the terms (i.e. in each term we 
    // have functions of the IVs in two places). Thus, we times by 2
    // but we do not neccesarily limit ourselves to two places.
    // g_type = 4 is the hybrid term proposed by Depei
    if (g_type == 1 || g_type == 2 || g_type == 3 || g_type == 4) {

        // Double terms as there is a coefficient and exponent to consider
        base_terms = 2 * base_terms;

    }

    // Determine the number of base terms in each term
    base_terms_per_term = 1;
    if(g_type == 1 or g_type == 2 || g_type == 3)
        base_terms_per_term = 2;
    // g_type=4 by Depei has following 3 terms per base term: f(x)* cos(f(x)) * e^(f(x))
    else if (g_type == 4) base_terms_per_term = 3; 
    
}

/**
 * Identify the term number that the base term belongs to
 * 
 * @param base_term     Base term to determine the term of
 * 
 * @return              size_t
 */
size_t MatrixContinuedFraction::term_no(int base_term) {
    
    size_t baseTermDepth = (base_term)/(base_terms_per_term);
    return baseTermDepth;

}

/**
 * Identify the depth that a base term belongs to
 * 
 * @param base_term     Base term to determine the term of
 * 
 * @return              size_t
 */
size_t MatrixContinuedFraction::term_depth(int base_term) {

    // Each depth of the continued fraction    
    size_t termsPerDepth = 2;

    size_t depth = 0;
    if( base_term>=2 ) {
        depth = 1+floor( (base_term-termsPerDepth)/(termsPerDepth*2) );
    }
    return depth;
}

/**
 * Caclulate a term given a set of variable values
 * 
 * @param term      Term within the fraction to evaluate
 * 
 * @return          double
 */
double MatrixContinuedFraction::eval_term(const vector<double>& vals, int term) {
    
    double ret = 0;

    // For all variables in ther term that are not masked
    for (size_t i = 0; i < vars; ++i) {

        // Multiply with the coefficient
        if (active[term][i]) {
            ret = Helper::add(ret, Helper::times(vals[i], coeffs[term][i]));
        }

    }
    
    // Add the constant if its provided
    if(active[term][vars])
        ret = Helper::add(ret, coeffs[term][vars]);

    // Force terms to evaluate as ceil/round/floor as required
    if(g_int_terms == 1)    return floor(ret);
    if(g_int_terms == 2)    return round(ret);
    if(g_int_terms == 3)    return ceil(ret);

    return ret;
}

/**
 * Return number of variables used and active at any point in the fraction
 * 
 * @return          int
 */
int MatrixContinuedFraction::used_features() {
    /*
     * Determine the number of used features used in at least one base term
     */

    // Vector for the independent variables set to false by default
    vector<bool> indpVars(vars, false);

    // For all base terms and varaibles
    for (size_t baseTerm = 0; baseTerm < base_terms; baseTerm++)
        for (size_t iv = 0; iv < vars; ++iv)

            // If the independent varaible is not masked and the coefficient is non-zero 
            if (Helper::almost_equal(coeffs[baseTerm][iv], g_ZERO, g_PREC) == false && active[baseTerm][iv])

                //Ensure the variable is included
                indpVars[iv] = true;

    // Return the number of variables on at any depth
    int ret = 0;
    for( size_t i = 0; i < indpVars.size(); i++) {

        if( indpVars[i] )
            ret++;
    }
    return ret;
}

/**
 * Evaluate a fraction to a specific depth
 * 
 * @param depth     Depth to evaluate to
 *                  Defaults to depth
 * 
 * @return          double
 */
double MatrixContinuedFraction::eval(const vector<double>& vals, int to_depth) {

    double ret = 0;

    try {

        // Override depth if fraction is not that deep
        if ( (size_t) to_depth > depth) 
            to_depth = depth;

        // Identify the number of terms to compute
        unsigned terms = 2 * depth;
        if(g_type == 1 || g_type == 2 || g_type == 3) 
            terms = 4 * depth + 1;

        // Generate an array of terms
        double term_vals[terms+1];
        for (unsigned term = 0; term <= terms; term++) {

            term_vals[term] = eval_term(vals, term);
            
        }

        // Calculate Traditional CFR
        if (g_type == 0) {

            ret = term_vals[terms];
            for (int term = terms - 2; term >= 0; term = term - 2) {
                ret = Helper::add(term_vals[term], Helper::divide(term_vals[term + 1], ret));
            }
        }
        // Caclulate Continued Logarithm
        else if (g_type == 1) {
            
            ret = Helper::times(term_vals[terms-1], Helper::pow(g_base, term_vals[terms]));
            for (int term = terms - 5; term >= 0; term = term - 4) {

                ret = Helper::times(term_vals[term], Helper::pow(g_base, term_vals[term + 1]))+
                        Helper::divide( 
                            Helper::times(term_vals[term + 2], Helper::pow(g_base, term_vals[term + 3])), ret);

            }
        } 
        // Calculate power law
        else if (g_type == 2 || g_type == 3) {

            ret =  Helper::pow(term_vals[terms-1], term_vals[terms]);
            for (int term = terms - 5; term >= 0; term = term - 4) {
                ret = Helper::pow(term_vals[terms], term_vals[terms+1]) +
                        Helper::divide( 
                            Helper::pow(term_vals[terms+2], term_vals[terms+3]), ret);
            }
        }
    } catch (...) {

        return numeric_limits<double>::max();

    }
    return ret; 
}



/**
 * Construct a continued fraction of the form
 *  cfr(x) = g0(x)+h0(x)/(h1(x)+g1(x)/h2(x)...)
 * 
 * With randomised coefficients
 * 
 * For the traditional CFR, each term gi(x), hi(x) is a linear function of x
 *  f(x) = c1x1+c2x2+c3x3+...+cNxN+c
 * 
 * Where {c1, c2, ..., cN} are the tunable parameters (coefficients) for each independent varaible x
 *
 * @param depth     Fraction depth (0 for a single term, incrementing 1 for each additional 2 terms). 
 *                  Defaults to the global depth value g_depth
 * 
 * @param vars      Number of independent variables in each linear function. 
 *                  Defaults to 1
 * 
 * @param alpha     Initalisation value for the constant of each term
 *                  Defaults to 0 which indicates to randomise the value
 * 
 * @return          MatrixContinuedFraction
 */
MatrixContinuedFraction::MatrixContinuedFraction(size_t frac_depth, size_t frac_vars, double frac_alpha) {

    if( g_debug )
        cout << "MatrixContinuedFraction(size_t depth, size_t vars, double alpha) dep:" << frac_depth << " vars:" << frac_vars << " alpha:" << frac_alpha << endl;

    depth = frac_depth;
    vars = frac_vars;

    error    = 0.0;
    fitness  = 0.0;
    penalty  = 1.0;
    freeParams = 0;
    // isRationalised = false; //turn OFF the flag of rationalisation
    
    init_terms();
    randomise(frac_alpha);
    
    if( g_origin ) force_origin();

}

/**
 * Construct a continued fraction from a solution file
 *
 * @param filename  Location of well-formatted JSON solution file 
 * 
 * @return          MatrixContinuedFraction
 */
MatrixContinuedFraction::MatrixContinuedFraction(string filename, bool rand_extra_depth) {
    // isRationalised = false; //turn OFF the flag of rationalisation at the beginning
    import_cfr(filename, rand_extra_depth);
}

/***
 * Copy constructor to Create a Copy of the MatrixContinuedFraction 
 * 
 */ 
MatrixContinuedFraction::MatrixContinuedFraction(const MatrixContinuedFraction &cp){
         // This will copy entire artay to new array
        memcpy(coeffs, cp.coeffs, sizeof(cp.coeffs));        
        memcpy(active, cp.active, sizeof(cp.active)); 

        // Use STL Copy Constructor
        residual_weight      = cp.residual_weight;
        temp_weight          = cp.temp_weight;
        feature_active       = cp.feature_active;

        // Copy variable values
        penalty              = cp.penalty;    
        error                = cp.error;              
        fitness              = cp.fitness;            
        freeParams           = cp.freeParams;         

        depth                =  cp.depth;      
        vars                 =  cp.vars;       
        base_terms           =  cp.base_terms; 
        base_terms_per_term  =  cp.base_terms_per_term;
    }

/**
 * Deconstruct a continued fraction and deallocate memory
 *
 * @return          void
 */
MatrixContinuedFraction::~MatrixContinuedFraction() {
    
    /* !!! Issues with deconstructor performing correctly

    // deallocate memory
    
    // delete coeffs
    for (int i = 0; i < base_terms; ++i) {
        free(coeffs[i]);
    }
    free(coeffs);
    
    //delete feature mask
    for (int i = 0; i < base_terms; ++i) {
        free(active[i]);
    }
    free(active);

    */
}


/**
 * Randomise an integer coefficient variable inclusively between (min, max) or (-max, -min)
 * 
 * @param term      Term to randomise
 * 
 * @param var       Variable within the term to randomise
 * 
 * @param min       Absolute minimum value
 *                  Defaults to 1
 * 
 * @param max       Absolute maximum value
 *                  Defaults to 3
 * 
 * @param chance    
 * 
 * @return          MatrixContinuedFraction
 */
void MatrixContinuedFraction::rand_var_int(size_t term, size_t var, int min, int max) {

    // determine coeff
    int coeff = g_rint(min, max);

    // Equal chance of negativity
    coeff *= g_rint(0, 1) ? 1 : -1;

    // Set the mask to show coefficient or constant
    active[term][var] = true;
    coeffs[term][var] = coeff;

    // We are unable to have 0 to a negative exponent, we should protect against this in the power law
    // This becomes -Inf (https://www.wolframalpha.com/input/?i=0%5E-0.5)
    if( g_type == 2) {

        // !!! Protect against negative exponent here
        return;

    }
}

/**
 * Randomise a real coefficient variable inclusively between (min, max)
 * 
 * @param term      Term to randomise
 * 
 * @param var       Variable within the term to randomise
 * 
 * @param min       Absolute minimum value
 *                  Defaults to 0
 * 
 * @param min       Absolute maximum value
 *                  Defaults to 1
 * 
 * @param chance    
 * 
 * @return          MatrixContinuedFraction
 */
void MatrixContinuedFraction::rand_var_real(size_t term, size_t var, double min, double max) {

    // Determine coeff
    double coeff = g_rreal(min, max);

    // We select only as positive

    // Set the mask to show coefficient or constant
    active[term][var] = true;
    coeffs[term][var] = coeff;
}

// Debug

/**
 * Return representation for power string based on the print type
 * 
 * @param print_type    Type of output desired
 * 
 * @return              string power symbol based on the print type 
 *                      ^ for Excel or Latex
 *                      ** for Numpy
 */
inline string power_str(ptype print_type) {

    if(print_type == Excel ||
       print_type == Latex)  return "^";
    
    if(print_type == Numpy)  return "**";

    return "";
}

/**
 * Return Excel cell name for a given column and row number
 * 
 * @param col       Excel column number
 * @param row       Excel row number
 * 
 * @return          string Excel name e.g.
 *                  (1,1) = A1, (2,2) = B2, (5,4) = E4 
 */
inline string excel_cell_name(int col, int row) {

    // Determine character given the column
    string str = "";
    while (col > 0) {
        str = (char) (65 + (col - 1) % 26) + str;
        col = (col - 1) / 26;
    }

    // Add row number
    str = str + to_string(row);

    return str;
}

/**
 * Output one base term of a fraction term given a print type, for example given a term
 *  gi(x) = f1(x)*e^*(f2(x)) we may be outputting f1(x) that may be a linear equation 
 *  in the form f1(x)=c1*x1+c2*x2+...cN*xN+c
 * 
 * @param out               Output stream to print to
 * @param names             Independent variable names
 * @param base_term         Base term to print
 * @param print_type        Type of output
 * @param do_print_rational Print coefficients as Rational format
 *
 * @return              void
 *                      
 */
void MatrixContinuedFraction::write_base_term(ostream& out, vector<string> names, int base_term, ptype print_type, bool do_print_rational) {

    bool needPlus = false;

    // Add IVs for the terbase_termm
    for(size_t i = 0; i <= vars; i++) {
        if(active[base_term][i]) {

            // Where coeffs are zero, don't bother displaying
            if(coeffs[base_term][i] != 0) {
            
                // Add plus when we are not the first coefficient. We only handle positive
                // coeffs as the output of the negative is handled by ostream
                if( needPlus ) {

                    if( coeffs[base_term][i] > 0 )
                        out << "+";
                    
                    // We omitt where abs(coeff) == 1 below, so in the case of -1 coeffient, we dont
                    // generate the negative sign. We cover this edge case here for the ivs but not the 
                    // constant at vars
                    if( (coeffs[base_term][i] == -1 && i != vars) || (coeffs[base_term][i]<0 && i != vars))
                        out << "-";
                }
                else if(coeffs[base_term][i] < 0 && i != vars){
                    out << "-";
                }

                // For the IVs times by the varaible, don't do for constant at last cell
                if(i != vars) {

                    if(do_print_rational) {
                            long long numer = 0, denom =0;
                            Helper::decimalToRational(abs(coeffs[base_term][i]), numer, denom);
                            cout << "(" << numer << "/" << denom  << ")";
                            if(print_type != Latex)      out << "*";
                    }else if( abs(coeffs[base_term][i]) != 1 ) {  // No need to show the '1' in 1*x1 
                        out << abs(coeffs[base_term][i]);
                        if(print_type != Latex)      out << "*";
                    }
                    
                    // Based on the print type generate the symbol
                    if( print_type == Numpy ) {
                        out << names[i+1] ;
                    }
                    else if( print_type == Latex ) {
                        out << "\\:\\text{" << names[i+1] << "}";
                    }
                    else if( print_type == Excel ) {

                        // We want the column as a 1-based number and to skipp the
                        // DV in column 1 therefore +2. We want row 2 as we assume headers
                        // in the first
                        out << excel_cell_name(2+i, 2);
                    }

                } else {
                    if(do_print_rational) {
                            long long numer = 0, denom =0;
                            Helper::decimalToRational(abs(coeffs[base_term][i]), numer, denom);
                            if (numer<0 || denom<0)
                                cout << "-(" << abs(numer) << "/" << abs(denom)  << ")";
                            else
                                cout << "(" << numer << "/" << denom  << ")";
                    }
                    else
                        out << coeffs[base_term][i];
                }
                
                needPlus = true;
            }      
        }
    }
}

/**
 * Output an entire term that may consist of multiple base terms to a specific print type
 * 
 * @param out               Output stream to print to
 * @param names             Independent variable names
 * @param base_term         Base term to print
 * @param print_type        Type of output
 * @param do_print_rational Print coefficients as Rational format
 * 
 * @return              void
 *                      
 */
void MatrixContinuedFraction::write_term(ostream& out,vector<string> names, int term, ptype print_type, bool do_print_rational) {

    if(print_type == Excel) {
        if( g_int_terms == 1)   out << "FLOOR(";
        if( g_int_terms == 2)   out << "ROUND(";
        if( g_int_terms == 3)   out << "CEILING(";
    }
    
    out << "(";
    write_base_term(out, names,term,print_type, do_print_rational);
    out << ")";

    // Continued log uses global base argument
    if(g_type == 1) {
        out << "*";
        out << g_base;

        if(print_type  == Latex)     out << "{";
        out << power_str(print_type);
        out << "(";
        write_base_term(out, names,term+1,print_type, do_print_rational);
        out << ")";
        if(print_type  == Latex)     out << "}";
    }

    // Power law has no base
    if(g_type == 2 || g_type == 3) {
        out << power_str(print_type);
        if(print_type  == Latex)     out << "{";
        out << "(";
        write_base_term(out, names,term+1,print_type, do_print_rational);
        out << ")";
        if(print_type  == Latex)     out << "}";
    }

    if(print_type == Excel) {
        if( g_int_terms == 1)   out << ",1)";
        if( g_int_terms == 2)   out << ",0)";
        if( g_int_terms == 3)   out << ",1)";
    }

}

/**
 * Output a fraction in table format
 * 
 * @param out           Output stream to print to
 * @param names         Independent variable names
 * @param base_term     Base term to print
 * @param print_type    Type of output
 * 
 * @return              void
 *                      
 */
void MatrixContinuedFraction::write_term_table(ostream& out, vector<string> names, ptype print_type) {

    int pad = 10;
    
    size_t lastTerm = 0;
    size_t firstBaseTerm = 0;

    // Debug CFR data
    out << endl << setw(pad) << "depth" << setw(pad) << "term" << setw(pad) << "base" << setw(pad) << "varNo" << setw(pad+g_PREC) << "coeff" << setw(pad) << "used?"  << setw(pad) << "feat" << endl;
    for(size_t baseTerm = 0; baseTerm < base_terms; baseTerm++) {

        size_t depth = term_depth(baseTerm);
        size_t term = term_no(baseTerm);

        // Print out each g_0, h_0, g_1,... etc. for debugging
        if(term != lastTerm) {
            write_term(out, names,firstBaseTerm,print_type);
            firstBaseTerm = baseTerm;
            lastTerm = term;
            out << endl;
        }

        for(size_t iv = 0; iv <= vars; iv++) {

            out << setw(pad) << depth;      // Depth

            // Term name
            string sTerm = "";
            if( term == 0)              sTerm += "g0(x)";
            else if( term == 1)         sTerm += "h0(x)";

            // Here we use g and h not gi, g(i+1) so we need to /2 again
            else if( term % 2 == 0)     sTerm = "g" + to_string(term/2)+"(x)";
            else                        sTerm = "h" + to_string(term/2)+"(x)";
            out << setw(pad) << sTerm;

            out << setw(pad) << baseTerm;                           // Base term
            out << setw(pad) << baseTerm*(1+vars)+iv;            // Var no
            
            out << setw(pad+g_PREC) << coeffs[baseTerm][iv];    // Coeff
            out << setw(pad) << active[baseTerm][iv];     // Mask

            // Feature name or constant
            if( iv == vars ){
                out << "\t" << "alpha" << endl; 
                break;   
            }                  
            else
                out << "\t" << names[iv+1] << endl;

        }
    }
    write_term(out, names,firstBaseTerm,print_type);
    out << endl;

}

/**
 * Output a fraction in evaluation form
 * 
 * @param out               Output stream to print to
 * @param names             Independent variable names
 * @param base_term         Base term to print
 * @param print_type        Type of output
 * @param do_summary        Print fitness, error, penalty etc
 * @param do_print_rational Print coefficients as Rational format
 * 
 * @return              void
 *                      
 */
void MatrixContinuedFraction::write_cfr(ostream& out, vector<string> names, ptype print_type, bool do_summary, bool do_print_rational) {

    size_t prec = cout.precision();
    
    //if( !do_summary )
    setprecision(g_PREC);

    size_t fullTerms = base_terms/base_terms_per_term;

    if( print_type == Latex)
        out << "$ ";

    // Dispaly CFR form
    for(size_t i = 0; i < fullTerms; i++) {

        // Print the single full term for g_0(x)
        if(i == 0) {
            write_term(out, names, i, print_type, do_print_rational);
            continue;
        }

        // for g_i terms excluding g_0 (i.e. even numbered full terms)
        if( i % 2 != 0) {

            out << "+";
            if( print_type == Latex)     cout << "\\cfrac{";
            write_term(out, names, i*base_terms_per_term, print_type, do_print_rational);
        }
        // for h_i terms (i.e. odd numbered full terms)
        else {

            if( print_type == Latex)     cout << "}";
            else                        out << "/";
            if( print_type == Latex)     cout << "{";
            out << "(";
            write_term(out, names, i*base_terms_per_term, print_type, do_print_rational);

        }

    }

    // Add braces to the end
    for(size_t i = 0; i < (fullTerms-1)/2; i++) {

        out << ")";
        if( print_type == Latex)     cout << "}";

    }
    // if( print_type == Latex)     cout << "}";
    if( print_type == Latex)
        out << " $";

    if( do_summary ) {
        out << endl;        
        out << "fitness:" << fitness << " error:" << error << " penalty:" << penalty << " used_feats:" << used_features();
        out << endl;
    }
    
    //if( !do_summary )
    setprecision(prec);

}


/**
 * Export fraction to a standardised JSON format usbale by import_cfr
 * 
 * @param solfile       Location to ouput JSON file
 * @param prefix        File prefix should multiple fractions be output
 * @param names         Independent variable names
 * 
 * @return              void
 *                      
 */
bool MatrixContinuedFraction::export_cfr(string solfile, int prefix, vector<string> names) {

    // We ignore file at the moment and set the location based on the code
    // For use within MTBTM we may want to modify this

    string soutFile = g_outPath + "Run"+to_string(prefix)+".Sol.json";  //Creating the solFile location as per execution number

    const char *sOut = soutFile.c_str();

    ofstream fout(sOut);

    if (!fout.is_open()) 
        return false;

    ostringstream oss;
    time_t ct = std::time(0);
    char* cc = ctime(&ct);
    oss << cc;
    string now_s = oss.str();

    string typeText;
    if( g_type == 0 )
        typeText = "Continued Fraction";
    if( g_type == 1 )
        typeText = "Continued Logarithm";
    if( g_type == 2 )
        typeText = "Power Law Continued Fraction";
    if( g_type == 2 )
        typeText = "Power Law Constant Exponent Continued Fraction";


    json all =  {
                    {"date", now_s},
                    {"type", g_type},
                    {"typeText", typeText},
                    {"seed", g_seed},
                    {"solution", {}},
                    {"equations", {}},
                };

    all["solution"]["objective"] = g_objec;
    all["solution"]["fitness"] = fitness;
    all["solution"]["error"] = error;
    all["solution"]["penalty"] = g_delta;
    if(g_objec=="aic"){
        all["solution"]["freeParams"] = freeParams;
    }
    all["solution"]["depth"] = depth;
    all["solution"]["terms"] = base_terms;
    all["solution"]["base_terms"] = base_terms_per_term;
    all["solution"]["vars"] = vars;
    all["solution"]["features"] = {};
    all["solution"]["variables"] = {};
    

    all["solution"]["global_mask"] = {};


    size_t var = 0;
    feature_active[vars] = true; //trun ON the Constant used in global_mask

    for(size_t term = 0; term < base_terms; term++) {
        for(size_t iv = 0; iv <= vars; iv++) {

            string name = iv == vars ? "constant" : names[iv+1];

            if( term == 0 ) {
                bool feat_mask = feature_active[iv] ? true : false;
                all["solution"]["global_mask"][iv] = feat_mask;
                all["solution"]["features"][iv] = name;
            }

            bool mask = active[term][iv] ? true : false;

            all["solution"]["variables"][var] = {
                {"baseTerm", term},
                {"varNo", iv},
                {"name", name},
                {"value", coeffs[term][iv]},
                {"active", mask}
            };
            var += 1;
        }
    }

    stringstream cfr_s;
    cfr_s.precision(g_PREC);

    write_cfr(cfr_s, names, Excel, false);
    all["equations"]["Excel"] = cfr_s.str();
    cfr_s.str("");

    write_cfr(cfr_s, names, Numpy, false);
    all["equations"]["Numpy"] = cfr_s.str();
    cfr_s.str("");

    write_cfr(cfr_s, names, Latex, false);
    all["equations"]["Latex"] = cfr_s.str();

    fout << all.dump(4) << endl;

    cout << endl << "Solution is written at: " << soutFile << endl;

    return true;
}

/**
 * Import fraction from a standardised JSON format created by export_cfr
 * 
 * @param solfile       Location to ouput JSON file
 * 
 * @return              void
 *                      
 */
bool MatrixContinuedFraction::import_cfr(string solfile, bool rand_extra_depth) {

    for( size_t i = 0; i < 41; i++) {
        for( size_t j = 0; j < 1002; j++) {
            coeffs[i][j] = 0;
            active[i][j] = false;
        }
    }

    ifstream ifs(solfile);
    json all = json::parse(ifs);

    if(g_verbose_mode > 0) {
        cout << "Importing file " << solfile << " of type " << all["typeText"] << endl;
        cout << "Excel: " << all["equations"]["Excel"] << endl;
        cout << "Numpy: " << all["equations"]["Numpy"] << endl;
        cout << "Latex: " << all["equations"]["Latex"] << endl;
    }

    g_type = all["type"];
    if(g_verbose_mode > 0)
        cout << "Type override to " << g_type << endl;

    base_terms_per_term = all["solution"]["base_terms"];
    vars = all["solution"]["vars"];

    if( g_incremental_depth==true && g_prev_model_as_var==true )
        vars++;

    if(g_objec=="aic")
        freeParams = all["solution"]["freeParams"];

    for (auto el : all["solution"]["variables"].items())
    {
        size_t term = el.value()["baseTerm"];
        size_t iv = el.value()["varNo"];
        string name = el.value()["name"];
        double coeff = el.value()["value"];
        bool is_active = el.value()["active"];
                
        coeffs[term][iv] = coeff; 
        active[term][iv] = is_active;        
    }

    feature_active = vector<bool>(vars + 1, false);
    feature_active[vars] = true;  // turn the Constant ON
    for( size_t i = 0; i < vars; i++) {
        feature_active[i] = all["solution"]["global_mask"][i];
    }

    depth = g_depth;    //update the depth value form cont_frac structure
    base_terms = 2 * depth + 1;
    if (g_type == 1 || g_type == 2 || g_type == 3) {
        base_terms = 2 * base_terms;
    }

    base_terms_per_term = 1;
    if(g_type == 1 or g_type == 2 || g_type == 3)
        base_terms_per_term = 2;
   
    // Process the extra terms
    for( size_t i = all["solution"]["terms"]; i < base_terms; i++) {
        for( size_t j = 0; j <= vars; j++ )  {

            if( rand_extra_depth ) {

                if(!feature_active[j] && j != vars)
                    continue;

                if( g_rreal() < 0.5 )   // randomly initialise terms for extra depth
                    rand_var_int(i,j); 

            } else {    // keep the same fitness score as solution from Previous Depth: by setting extra terms as (0)/(1)

                if( j == vars && i != all["solution"]["terms"]) {   //making denominator term=(1) by setting constant=1
                    coeffs[i][j] = 1; 
                    active[i][j] = true;
                } else {    // make numerator term=(0) by setting all coefficeints to zero
                    coeffs[i][j] = 0; 
                    active[i][j] = true;
                }
            }
        }
    }
    

    return true;
}


int MatrixContinuedFraction::sum_param_count(vector<int> paramCounter ) {
    /*
     * Determine the sum of the free paraemeters used in the entire fraction
     *  free parameters: non zero coefficients and teh constant
     *  paramCounter :  // keep the number of parameters [coefficeints with feature/meta-feature plus constant] in the model
     */

    int paramCount = 0;
    // find the cost of meta-features
    vector<bool> metaConsidered = vector<bool>(vars+1, false);

    // For all base terms and varaibles
    for (size_t baseTerm = 0; baseTerm < base_terms; baseTerm++)
        for (size_t iv = 0; iv <= vars; ++iv)
        {
            // If the independent varaible is not masked and the coefficient is non-zero 
            if (Helper::almost_equal(coeffs[baseTerm][iv], g_ZERO, g_PREC) == false && active[baseTerm][iv])
            {
                //Increament the parameter count as per value in the vector
                if(paramCounter[iv]>0){
                    // consider the cost of the meta feature once
                    if(metaConsidered[iv]==false)
                    {
                        paramCount = paramCount+ paramCounter[iv];
                        metaConsidered[iv] = true;
                    }
                    paramCount +=1;  //cost of associated coefficent of meta-feature
                }
                else  paramCount +=1; //cost of associated coefficent/const
             }
        }
    return paramCount;
}
