
#include "local_search.hpp"

optimize::optimize(const data_store& td, MatrixContinuedFraction& f, string objective, double del) :
    train_data(td), frac(f), objec(objective), delta(del) {
    /*
     *
     */

     //run linreg first
    if (g_rreal(0.0, 1.0) < g_linreg_chance) {
        MatrixContinuedFraction frac(f);
        data_store data(td);
        LinearRegression linreg(frac, (int) f.depth, data);
        linreg.fit();
        f = linreg.fraction;
    }
    ndim = 0;

    size_t start_term = 0;
    if( g_depth_lock != -1)
        start_term = 2*g_depth_lock+1;

    for (size_t i = start_term; i < frac.base_terms; ++i) {
        for (size_t j = 0; j < frac.vars +1; ++j) { //frac.repr[i].coeff.size() is var+1
            if (frac.active[i][j] == true) { // feature is on
                ++ndim;
                var_map.push_back({i, j});
            }
        }
    }
}

double optimize::run() {
    /*
     *
     */

    double best = -1.0;

    if (g_run_serial == true) {
        best = runSerial();
        if (g_verbose_mode == 3) {
            cout << endl;
        }
    } else
        best = runParallel();
    return best;
}

struct Compare {
    double val;
    MatrixContinuedFraction res;
};

#pragma omp declare reduction(minimum : struct Compare : omp_out = omp_in.val < omp_out.val ? omp_in : omp_out) initializer (omp_priv = omp_orig)

double optimize::runParallel() {
     /*
     * Local search parallel execution
     */
    
    double best;
    struct Compare minres = {numeric_limits<double>::max(), frac};
    int c;

    //non-window version
    if (g_window_size < g_ZERO) {
#pragma omp parallel for shared(train_data, g_percnt_data) firstprivate(frac) reduction(minimum:minres)
        for (c = 0; c < g_nm_runs; c++) {

            MatrixContinuedFraction buf;
            double tmp;

            buf = frac;
            evaluator e(train_data, g_percnt_data);

            nelder_mead(buf, e);
            tmp = e.eval_fit_full(buf);

            if (tmp < minres.val) {
                minres.val = tmp;
                minres.res = buf;
                frac = minres.res;
            }
        }
        best = minres.val;
        frac = minres.res;

        return best;
    } else { // window version
#pragma omp parallel for shared(train_data, g_window_size) firstprivate(frac) reduction(minimum:minres)
        for (c = 0; c < g_nm_runs; c++) {

            MatrixContinuedFraction buf;
            double tmp;

            buf = frac;
            evaluator e(train_data, g_window_size);

            nelder_mead(buf, e);
            tmp = e.eval_fit_full(buf);

            if (tmp < minres.val) {
                minres.val = tmp;
                minres.res = buf;
                frac = minres.res;
            }
        }
        best = minres.val;
        frac = minres.res;

        return best;
    }
}

double optimize::runSerial() {
    /*
     * Local search serial execution
     */

    double best = numeric_limits<double>::max(); // was -1; Mohammad changed it

    //non-window 
    if (g_window_size < g_ZERO) {
        for (int c = 0; c < g_nm_runs; c++) {
            MatrixContinuedFraction buf = frac;
            evaluator e(train_data, g_percnt_data); //percentage of data
            nelder_mead(buf, e);

            double tmp = e.eval_fit_full(buf);
            if (tmp < best) { //not necessary to check best < g_ZERO ||
                best = tmp;
                frac = buf;
            }
        }
    } else { //window of x
        for (int c = 0; c < g_nm_runs; c++) {
            MatrixContinuedFraction buf = frac;
            evaluator e(train_data, g_window_size); // initialise the evaluator for the window
            nelder_mead(buf, e);
            double tmp = e.eval_fit_full(buf);
            if (tmp < best) { //not necessary to check best < g_ZERO ||
                best = tmp;
                frac = buf;
            }
        }
    }
    if (g_verbose_mode == 3) {
        cout << endl;
    }
    return best;
}

double optimize::eval_fit(vector<double>& vec, MatrixContinuedFraction& buf, const evaluator& e) const {
    /*
     * Evaluation of CFR in local search
     */
    
    for (int i = 0; i < ndim; ++i) {

        const pii& pos = var_map[i];
       
        // When using the power law representation, we can't have a negative exponent without a constant, i.e. where we have 0^(c) where c is negative
        // For ease we ensure that c is positive in the exponent everywhere to easily avoid the violation. NB this is especially troublesome when
        // forcing the function through the origin
    
        // force integer
        if( g_int_only )        
            vec[i] = round(vec[i]);

        if( (g_type == 2 || g_type == 3) &&        // Power-law CFR and
            (size_t)pos.second == buf.vars &&    // Is the consant term and
            vec[i] < 0 ) {                          // Is proposed to be negative by the local search
           
            // For the exponent terms for, enforce a non-negative constant exponent
            if( ((pos.first+1) % (int)buf.base_terms_per_term) == 0 )
                vec[i] = 0;

            // For the constants in the leading coefficients anywhere, enforce a positive
            // constant so we do not have -a^(c) where c is a decimal; i.e. a complex number
            // !!! Can't replicate this issue e.g. the follow does not have a complex value https://www.wolframalpha.com/input/?i=-5%5E%283.24%29
            //if( ((pos.first) % (int)buf.base_terms_per_term) == 0 )
            //    vec[i] = 0;

        }
   
        buf.coeffs[pos.first][pos.second] = vec[i];
    }
    return e.eval_fit(buf);
}

double optimize::nelder_mead(MatrixContinuedFraction& buf, const evaluator& e) const {
    /*
     * Speed-focused simplified nelder mead as defined in Algorithm 3 of 
     * https://direct.mit.edu/evco/article/25/3/351/1046/Evolving-a-Nelder-Mead-Algorithm-for-Optimization
     * A secondary reference exists if needed http://www.scholarpedia.org/article/Nelder-Mead_algorithm#Initial_simplex
     * Implements reflection (project the worst point in between the two best)
     * Expansion (extend reflection projection if we result in a better point than the two best)
     * Contraction (shrink reflection projection if we result in a worst point than the two best)
     */

    using coord = vector<double>;

    // Reflection
    auto refl = [this](const coord& a, const coord & b) {

        // a and b represent the centroid between the two best points
        // therefore we multiply twice their difference to cover the distance
        // from the worst point to the centroid, then further the same distance
        // beyond. See below where
        //
        // - xW is the worst point
        // - x1, x2 are the best two points
        // - c is the centroid
        // - r is the reflected point, 2x the distance of xW to c
        // 
        //              x1
        //  
        //   xW         c          r
        // 
        //              x2

        coord ret(ndim);
        for (int i = 0; i < ndim; ++i) {
            ret[i] = 2 * a[i] - b[i];
        }
        return ret;
    };

    // Expansion
    auto expa = [this](const coord& a, const coord & b) {
        coord ret(ndim);
        for (int i = 0; i < ndim; ++i) {
            ret[i] = 3 * a[i] - 2 * b[i];
        }
        return ret;
    };

    // Contraction
    auto contr = [this](const coord& a, const coord & b) {
        coord ret(ndim);
        for (int i = 0; i < ndim; ++i) {
            ret[i] = 0.5 * (a[i] + b[i]);
        }
        return ret;
    };

    // Simplex object with fitness value and the associated values for each dimension
    // greater is the sorting mechanism 
    multimap<double, coord, greater<double>> simplex;
    
    // Simplex step size
    double step = 0.5;

    //if(g_int_terms > 0 || g_int_only)
    //    step = 5.0;

    // Centroid
    coord cent(ndim);      
    double cent_fit;

    // Extract values from the current continued fraction
    coord tmp(ndim);
    for (int i = 0; i < ndim; ++i) {
        const pii& pos = var_map[i];
        tmp[i] = frac.coeffs[pos.first][pos.second];
    }

    // Evaluation of the simplex at the current point
    simplex.insert({eval_fit(tmp, buf, e), tmp});

    // Create simplex points by steping forward in each dimension a length of 'step' and appending
    // to the simplex. After this operation our simplex will be ndim+1 in length, as we have the original
    // fitness and values, then 10 modifications where each dimension is stepped
    for (int i = 0; i < ndim; ++i) {
        tmp[i] += step;
        simplex.insert({eval_fit(tmp, buf, e), tmp});
        tmp[i] -= step;
    }

    // Given our array of simplex points, calculate the centroid point and its fitness
    for (int i = 0; i < ndim; ++i) {
        cent[i] = 0;
        for (auto it = ++simplex.begin(); it != simplex.end(); ++it)

            // Sum the values from each dimension
            // Caclulate the average value on the fly to minimise risk of overflow
            cent[i] += it->second[i] /= ndim;
        
    }
    cent_fit = eval_fit(cent, buf, e);

    int iter = 0;       // either converge, or reach max number of iterations, or stagnate for too long
    int stag = 0;       // stagnation is when the new vertex is still the worst

    while ( simplex.begin()->first - (--simplex.end())->first > g_EPS  &&   // The range in the simplex points has a difference > some epsilon   
            iter < g_nm_max_iters &&                                        // We have not reached the max iterations
            stag < g_nm_reset_stag                                          // We have not stagnated
        ) {
        
        // Get the worst point i.e. highest fiteness 
        coord& vw = simplex.begin()->second;
        double vw_fit = simplex.begin()->first;
        
        // Get the second worst point
        double vsw_fit = (++simplex.begin())->first;

        // Get the best point
        coord& vb = (--simplex.end())->second;

        // Get the fitness after reflection
        coord vr = refl(cent, vw);
        double vr_fit = eval_fit(vr, buf, e);

        // Get the point of expansion from the word and centroid points
        coord ve = expa(cent, vw);

        coord vtmp;
        // If reflection is better than our worst point, try expanding and
        // settle on the best result
        if (vr_fit < vw_fit) {

            if (eval_fit(ve, buf, e) < cent_fit) {      
                if (vr_fit < cent_fit)
                    vtmp = contr(vr, ve);
                else
                    vtmp = contr(refl(cent, vr), cent);
            } else
                vtmp = cent;

            vtmp = contr(contr(vtmp, cent), ve);

        } else {

            if (vsw_fit < cent_fit)
                vtmp = vb;
            else
                vtmp = cent;

            vector<double> ncvec = contr(refl(cent, ve), contr(vtmp, cent));
            if (eval_fit(ncvec, buf, e)
                    < cent_fit)
                vtmp = refl(cent, vr);
            else
                vtmp = cent;

            //moh: floating point comparison is used here
            if (eval_fit(vtmp, buf, e) < cent_fit)
                vtmp = vr;
            else
                vtmp = contr(cent, vw);

            vtmp = contr(vw, contr(vtmp, cent));
        }

        // replace the worst
        simplex.erase(simplex.begin());
        
        //moh: floating point comparison is used here
        double vtmp_fit = eval_fit(vtmp, buf, e);
        if (vtmp_fit >= simplex.begin()->first) ++stag;
        else stag = 0;

        simplex.insert({eval_fit(vtmp, buf, e), vtmp});

        // recompute the centroid
        for (int i = 0; i < ndim; ++i) {
            cent[i] = 0;
            for (auto it = ++simplex.begin(); it != simplex.end(); ++it)
                cent[i] += it->second[i] / ndim;
        }
        cent_fit = eval_fit(cent, buf, e);

        ++iter;
    }

    double ret_fit_full = e.eval_fit_full(buf);
    if (g_verbose_mode == 3) {
        cout << " --> " << ret_fit_full; //moh: temp print
    }
    return ret_fit_full;
}
