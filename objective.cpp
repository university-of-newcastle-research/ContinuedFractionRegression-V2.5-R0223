
#include "objective.hpp"

double ObjectiveFuncs::mse(MatrixContinuedFraction& frac, const data_store& data, const unordered_set<int>& subset, const double multi_const) {

    /// Perform MSE evaluation of frac on a subset of points within the data

    double squareErrorSum = 0.0,
            weightedErrorSum = 0.0;

    // for the subset of elements in data
    for (int i : subset) {

        double err = (data.expected[i] - frac.eval(data.input[i]));
        double squareError = err*err;
        squareErrorSum += (squareError); // * data.expected[i]); // for the time being
        if (g_weight_present) weightedErrorSum += (squareError * data.samp_weight[i]);

    }

    if (!g_weight_present) weightedErrorSum = squareErrorSum;

    frac.penalty = multi_const;
    frac.error = squareErrorSum / subset.size();
    frac.fitness = weightedErrorSum / subset.size() * multi_const;

    return frac.fitness;

}

double ObjectiveFuncs::mwe(MatrixContinuedFraction& frac, const data_store& data, const unordered_set<int>& subset, const double multi_const) {

    /// Perform Mean Weighted Absolute Error evaluation of cfrac on a subset of points within the data
    double squareErrorSum = 0.0,
            weightedErrorSum = 0.0;

    // for the subset of elements in data
    for (int i : subset) {
        double pred = frac.eval(data.input[i]);
        double err = fabs(data.expected[i] - pred);
        weightedErrorSum += (err * data.samp_weight[i]);

        double squareError = err*err;
        squareErrorSum += (squareError);

    }

    frac.penalty = multi_const;
    frac.error = squareErrorSum / subset.size();
    frac.fitness = weightedErrorSum / subset.size() * multi_const;

    return frac.fitness;

}

double ObjectiveFuncs::uwae(MatrixContinuedFraction& frac, const data_store& data, const unordered_set<int>& subset, const double multi_const) {

    /// Perform Absolute Error weighted for Universe data
    //  Pablo, 2:48 PM on 3/11/2020:
    //      abs = |y_obs - y_fit|
    //      sum_i=1^N y_obs(i) / dy (uncertainty)(i) * abs(i)
    //

    double squareErrorSum = 0.0, weightedErrorSum = 0.0;

    // for the subset of elements in data
    for (int i : subset) {
        double pred = frac.eval(data.input[i]);
        double abserr = fabs(data.expected[i] - pred);
        double score = (data.expected[i] / data.dy[i]); // dy(uncertainty) should be passed as weight
        weightedErrorSum += (score * abserr); // y_obs(i) / dy (uncertainty)(i) * abs(i)
        squareErrorSum += (abserr * abserr);
    }

    frac.penalty = multi_const;
    frac.error = squareErrorSum / subset.size();
    frac.fitness = weightedErrorSum * multi_const;

    return frac.fitness;

}

double ObjectiveFuncs::abserr(MatrixContinuedFraction& frac, const data_store& data, const unordered_set<int>& subset, const double multi_const) {

    /// Perform Absolute Error for Universe data
    //  Pablo, 10:45 AM on 6/11/2020:
    //      abs = |y_obs - y_fit|
    //      sum_i=1^N abs(i)
    //

    double squareErrorSum = 0.0, errorSum = 0.0;

    // for the subset of elements in data
    for (int i : subset) {
        double pred = frac.eval(data.input[i]);
        double abserr = fabs(data.expected[i] - pred);
        errorSum += (abserr);
        squareErrorSum += (abserr * abserr);
    }

    frac.penalty = multi_const;
    frac.error = squareErrorSum / subset.size();
    frac.fitness = errorSum * multi_const;

    return frac.fitness;
}

double ObjectiveFuncs::mase(MatrixContinuedFraction& frac, const data_store& data, const unordered_set<int>& subset) {

    /// Perform 'mean absolute scaled error' evaluation of frac on a subset of points within the data

    //compute value of numerator
    double sum_forecast_err = 0.0;
    for (int i : subset) {
        double actual_val = data.expected[i]; // Y_i
        double forecast_val = frac.eval(data.input[i]); // F_i
        // e_i = Y_i - F_i
        double forecast_err = g_ZERO;
        if (actual_val < g_ZERO && forecast_val < g_ZERO)
            forecast_err = fabs(fabs(actual_val) - fabs(forecast_val));
        else if (actual_val < g_ZERO && forecast_val > g_ZERO)
            forecast_err = fabs(actual_val) + forecast_val;
        else if (actual_val > g_ZERO && forecast_val < g_ZERO)
            forecast_err = actual_val + fabs(forecast_val);
        else forecast_err = fabs(actual_val - forecast_val);

        sum_forecast_err += (forecast_err);
    }
    double nom = sum_forecast_err / subset.size();

    //compute value of denominator
    double sum_ae = 0.0;
    bool first = true;
    int tp;
    for (int i : subset) {
        if (first == false) {
            int t = i; // t_i
            double Y_t = data.expected[t];
            double Y_prev = data.expected[tp];
            double abs_err = g_ZERO;
            if (Y_t < g_ZERO && Y_prev < g_ZERO)
                abs_err = fabs(fabs(Y_t) - fabs(Y_prev));
            else if (Y_t < g_ZERO && Y_prev > g_ZERO)
                abs_err = fabs(Y_t) + Y_prev;
            else if (Y_t > g_ZERO && Y_prev < g_ZERO)
                abs_err = Y_t + fabs(Y_prev);
            else abs_err = fabs(Y_t - Y_prev);

            sum_ae += abs_err;
            tp = i;
        } else {
            tp = i;
            first = false;
        }
    }
    double denom = sum_ae / (subset.size() - 1);

    // compute MASE
    return nom / denom;

}

double ObjectiveFuncs::medae(MatrixContinuedFraction& frac, const data_store& data, const unordered_set<int>& subset) {

    /// Perform 'median absolute error' evaluation of frac on a subset of points within the data

    // median_absolute_error
    double med_abs_error = 0.0;
    vector<double> abs_errors;
    for (int i : subset) {
        // residual = Yi - Y'i
        double pred_val = frac.eval(data.input[i]);
        double meas_val = data.expected[i];
        double abs_err = g_ZERO;
        if (meas_val < g_ZERO && pred_val < g_ZERO)
            abs_err = fabs(fabs(meas_val) - fabs(pred_val));
        else if (meas_val < g_ZERO && pred_val > g_ZERO)
            abs_err = fabs(meas_val) + pred_val;
        else if (meas_val > g_ZERO && pred_val < g_ZERO)
            abs_err = meas_val + fabs(pred_val);
        else abs_err = fabs(meas_val - pred_val);

        abs_errors.push_back(abs_err);
    }
    sort(abs_errors.begin(), abs_errors.end());
    if (abs_errors.size() % 2 == 0)
        med_abs_error = (abs_errors[abs_errors.size() / 2 - 1] + abs_errors[abs_errors.size() / 2]) / 2;
    else
        med_abs_error = abs_errors[abs_errors.size() / 2];

    return med_abs_error;

}

double ObjectiveFuncs::mederr(MatrixContinuedFraction& frac, const data_store& data, const unordered_set<int>& subset) {

    /// Perform 'median absolute error' evaluation of frac on a subset of points within the data

    double med_error = 0.0;
    vector<double> errors;
    for (int i : subset) {

        double pred_val = frac.eval(data.input[i]);
        double meas_val = data.expected[i];
        double err = (meas_val - pred_val);

        errors.push_back(err);
    }
    sort(errors.begin(), errors.end());
    if (errors.size() % 2 == 0)
        med_error = (errors[errors.size() / 2 - 1] + errors[errors.size() / 2]) / 2;
    else
        med_error = errors[errors.size() / 2];

    return med_error;

}

double ObjectiveFuncs::emc(MatrixContinuedFraction& frac, const data_store& data, const unordered_set<int>& subset) {

    /// Perform 'excursion matching criteria' evaluation of frac on a subset of points within the data
    ///  https://contextearth.com/2017/10/25/improved-solver-target-error-metric/ 

    double emc;
    double nom_sum = 0.0, denom_sum = 0.0;
    for (int i : subset) {
        double pred_val = frac.eval(data.input[i]);
        double meas_val = data.expected[i];
        nom_sum = nom_sum + (meas_val * pred_val);
        denom_sum = denom_sum + (fabs(meas_val) * fabs(pred_val));
    }
    emc = nom_sum / denom_sum;
    return 1.0 - emc;

}

double ObjectiveFuncs::theil(MatrixContinuedFraction& frac, const data_store& data, const unordered_set<int>& subset) {

    /// Perform 'Minimum Theil index' evaluation of frac on a subset of points within the data

    vector<double> slopes;
    for (int i : subset) {
        for (int j : subset) {
            if (i != j) {
                double pred_val_i = frac.eval(data.input[i]);
                double meas_val_i = data.expected[i];
                double pred_val_j = frac.eval(data.input[j]);
                double meas_val_j = data.expected[j];

                double slope = (pred_val_j - pred_val_i) / (meas_val_j - meas_val_i);
                slopes.push_back(slope);
            }
        }
    }
    size_t size = slopes.size();
    sort(slopes.begin(), slopes.end());

    if (size % 2 == 0) return (slopes[size / 2 - 1] + slopes[size / 2]) / 2;
    else return slopes[size / 2];

}

double ObjectiveFuncs::pearson(MatrixContinuedFraction& frac, const data_store& data, const unordered_set<int>& subset, const double multi_const) {

    /// Perform 'Minimum Theil index' evaluation of frac on a subset of points within the data

    double sum_sq_err = 0.0, p_sum = 0.0, m_sum = 0.0, ret = 0.0, pred_avg = 0.0, meas_avg = 0.0, mse = 0.0;
    for (int i : subset) {
        double pred_val = frac.eval(data.input[i]);
        double meas_val = data.expected[i];
        double err = (meas_val - pred_val);
        double sq_err = err * err;
        sum_sq_err += sq_err;
        p_sum += pred_val;
        m_sum += meas_val;
    }

    mse = sum_sq_err / subset.size();
    pred_avg = p_sum / subset.size();
    meas_avg = m_sum / subset.size();

    double pearson = 0.0, sum_num = 0.0, meas_var = 0.0, pred_var = 0.0;

    for (int i : subset) {

        double pred_val = frac.eval(data.input[i]);
        double meas_val = data.expected[i];
        double meas_diff = (meas_val - meas_avg);
        double pred_diff = (pred_val - pred_avg);

        sum_num += pred_diff*meas_diff;
        meas_var += meas_diff*meas_diff;
        pred_var += pred_diff*pred_diff;
    }
    pearson = sum_num / pow(meas_var*pred_var, 0.5);
    ret = 1.0 - pearson; //multi_const * mse * (1.0 - pearson); // with Terry

    // saving fraction member value
    frac.penalty = multi_const;
    frac.error = mse;
    frac.fitness = ret;

    return ret;

}

double ObjectiveFuncs::nmse(MatrixContinuedFraction& frac, const data_store& data, const unordered_set<int>& subset, const double multi_const) {

    /// Perform 'normalised minimum square error' evaluation of frac on a subset of points within the data

    double sum_sq_err = 0.0, p_sum = 0.0, m_sum = 0.0, ret = 0.0, meas_avg = 0.0, mse = 0.0, nmse = 0.0;

    for (int i : subset) {

        double pred_val = frac.eval(data.input[i]);
        double meas_val = data.expected[i];
        double err = (meas_val - pred_val);
        double sq_err = err * err;
        sum_sq_err += sq_err;
        p_sum += pred_val;
        m_sum += meas_val;
    }
    mse = sum_sq_err / subset.size();
    meas_avg = m_sum / subset.size();

    double m_sum_var = 0.0;

    for (int i : subset) {
        double var_i = (data.expected[i] - meas_avg);
        m_sum_var += var_i * var_i;
    }
    double m_variance = m_sum_var / (subset.size() - 1);
    nmse = mse / m_variance;
    ret = nmse * multi_const;

    // saving fraction member value
    frac.penalty = multi_const;
    frac.error = mse;
    frac.fitness = ret;

    return ret;

}

double ObjectiveFuncs::ub(MatrixContinuedFraction& frac, const data_store& data, const unordered_set<int>& subset, const double multi_const) {

    /// Perform 'upper bound' evaluation of frac on a subset of points within the data

    double sum_sym_err = 0.0, sum_sq_err = 0.0, ret = 0.0;
    for (int i : subset) {
        // residual = Yi - Y'i
        double pred_val = frac.eval(data.input[i]);
        double meas_val = data.expected[i];
        double err = (meas_val - pred_val);
        sum_sq_err += (err * err);

        if (err <= g_ZERO) {
            err = 1.0 + fabs(err);
        } else {
            double errval = (1.0 + fabs(err));
            err = pow(errval, 8);
        }
        sum_sym_err = sum_sym_err + err;
    }

    double sum_err = sum_sym_err / subset.size();
    ret = sum_err * multi_const;
    // saving fraction member value
    frac.penalty = multi_const;
    frac.fitness = ret;
    frac.error = sum_sq_err / subset.size();

    return ret;

}

double ObjectiveFuncs::lb(MatrixContinuedFraction& frac, const data_store& data, const unordered_set<int>& subset, const double multi_const) {

    /// Perform 'lower bound' evaluation of frac on a subset of points within the data

    double sum_sym_err = 0.0, sum_sq_err = 0.0, ret = 0.0;
    for (int i : subset) {

        double pred_val = frac.eval(data.input[i]);
        double meas_val = data.expected[i];
        double err = (meas_val - pred_val);
        sum_sq_err += (err * err);

        if (err >= g_ZERO) {
            err = 1.0 + fabs(err);
        } else {
            double errval = (1.0 + fabs(err));
            err = pow(errval, 8);
        }
        sum_sym_err = sum_sym_err + err;
    }

    double sum_err = sum_sym_err / subset.size();
    ret = sum_err * multi_const;
    // saving fraction member value
    frac.penalty = multi_const;
    frac.fitness = ret;
    frac.error = sum_sq_err / subset.size();

    return ret;

}

double ObjectiveFuncs::cw(MatrixContinuedFraction& frac, const data_store& data, const unordered_set<int>& subset, const double multi_const) {

    // Perform Chi-Squared Error Metric with Uncertainity 
    // chi-squared = 1/n * \sum_{i=1}^{n} \frac{1}{sigma_i^2} (obs_i - pred_i)^2

    double squareErrorSum = 0.0,
            weightedErrorSum = 0.0;

    // for the subset of elements in data
    for (int i : subset) {

        double err = (data.expected[i] - frac.eval(data.input[i]));
        double squareError = err*err;
        squareErrorSum += (squareError);
        double norm_dy = 1/(data.dy[i] * data.dy[i]); // 1/sigma_i^2
        weightedErrorSum += (squareError * norm_dy);
    }

    frac.penalty = multi_const;
    frac.error = squareErrorSum / subset.size();
    frac.fitness = weightedErrorSum / subset.size() * multi_const;

    return frac.fitness;
}

double ObjectiveFuncs::cwy(MatrixContinuedFraction& frac, const data_store& data, const unordered_set<int>& subset, const double multi_const) {

    // Perform Chi-Squared Error Metric with Uncertainity/Weight as: (dy/y)^2
    // chi-squared = 1/n * \sum_{i=1}^{n} \frac{(obs_i - pred_i)^2} {(\frac{sigma_i}{obs_i})^2 }
    //             = 1/n * \sum_{i=1}^{n} {(obs_i - pred_i)^2} * {(obs_i / sigma_i)^2 }

    double squareErrorSum = 0.0,
            weightedErrorSum = 0.0;

    // for the subset of elements in data
    for (int i : subset) {

        double err = (data.expected[i] - frac.eval(data.input[i]));
        double squareError = err*err;
        squareErrorSum += (squareError);
        double norm_dy = (data.expected[i] / data.dy[i]);  // y/dy
        norm_dy = norm_dy * norm_dy; // (y/sigma_i)^2
        weightedErrorSum += (squareError * norm_dy); //  sq_err * (y/dy)^2 = sq_err /  (dy/y)^2
    }

    frac.penalty = multi_const;
    frac.error = squareErrorSum / subset.size();
    frac.fitness = weightedErrorSum / subset.size() * multi_const;

    return frac.fitness;
}


double ObjectiveFuncs::aic(MatrixContinuedFraction& frac, const data_store& data, const unordered_set<int>& subset) {

    /// Perform AIC evaluation of frac on a subset of points within the data

    double squareErrorSum = 0.0;
    double weightedErrorSum = 0.0;

    // for the subset of elements in data
    if (g_weight_present){
        for (int i : subset) {
            double err = (data.expected[i] - frac.eval(data.input[i]));
            double squareError = err*err;
            squareErrorSum += (squareError); // * data.expected[i]); // for the time being
            weightedErrorSum += (squareError * data.samp_weight[i]);
        }
    }else{
        for (int i : subset) {
            double err = (data.expected[i] - frac.eval(data.input[i]));
            double squareError = err*err;
            squareErrorSum += (squareError); // * data.expected[i]); // for the time being
        }
        weightedErrorSum = squareErrorSum;
    }
    
    // if (!g_weight_present) weightedErrorSum = squareErrorSum;

    // calculate the AIC by given formula 2.2.3 in http://faculty.wwu.edu/jmcl/Biostat/l_AIC4.pdf
    // where mse = RSS/n
    // AIC = n ln(mse) + 2K
    // --- Calculate the Parameter Counter
    // frac.write_cfr(cout, data.names, Excel, false);
    int param_count = frac.sum_param_count(data.aic_param_count);
    // cout << " :" << param_count << endl;
    
    int n = subset.size();
    double fit = n* log( weightedErrorSum / n ) + 2*param_count;

    if (param_count>0 && (n / param_count < 40))
        fit = fit + ((2 * param_count * (param_count + 1)) / (n - param_count - 1));

    // frac.write_cfr(cout, data.names, Numpy, false);
    // cout <<endl << "SqErr: " << weightedErrorSum << "\tn: " << n << "\tK: " << param_count << "\tAIC: " << fit << endl;


    frac.penalty = param_count;
    frac.error = squareErrorSum / subset.size();
    frac.fitness = fit;
    frac.freeParams = param_count;
    return frac.fitness;

}


double ObjectiveFuncs::mre(MatrixContinuedFraction& frac, const data_store& data, const unordered_set<int>& subset, const double multi_const) {

    /// Perform MSE evaluation of frac on a subset of points within the data

    double errSum = 0.0;

    // For all examples
    for (int i : subset) {

        // If zero, we sum nothing
        if(data.expected[i] != 0)
            errSum += abs((data.expected[i] - frac.eval(data.input[i]))/data.expected[i]);

    }

    frac.fitness = errSum / subset.size();

    return frac.fitness;

}


double ObjectiveFuncs::mwrae(MatrixContinuedFraction& frac, const data_store& data, const unordered_set<int>& subset, const double multi_const) {
    /// Perform MWRAE evaluation of frac on a subset of points within the data
    double squareErrorSum = 0.0,
            weightedErrorSum = 0.0;
    double error_no_weight_at_all = 0.0;
    
    const double EPS = 0.00000001;
   
    //ADDED THIS IN CASE INITIALIZATION DIDN'T HAPPEN FAST ENOUGH
    if (frac.residual_weight.size() == 0) {
        for (unsigned int i = 0; i < data.num_entry; i++) {
            frac.residual_weight[i] = 1;
            frac.temp_weight[i] = 1;
        }
    }
    int cur_weight = 1;
   
    //{abs(residual), datapoint index}
    vector<pair<double, int> > data_weights;
    for (auto it = frac.residual_weight.begin(); it != frac.residual_weight.end(); it++) {
        data_weights.push_back({it->second, it->first});
    }
    sort(data_weights.begin(), data_weights.end());
    double last_weight = 1.0;
    for(pair<double, int> curP : data_weights) {
        //Added this for generalization
        if (!subset.count(curP.second)) {
            continue;
        }
        double err = (data.expected[curP.second] - frac.eval(data.input[curP.second]));
        //BAD NAME vvvvv
        double squareError = fabs(err) / fabs(data.expected[curP.second]);
        error_no_weight_at_all += err*err;
        squareErrorSum += (squareError * cur_weight); // * data.expected[i]); // for the time being
        //cout << err << "      " << squareError << "      " << squareError * cur_weight << "        " << squareErrorSum << endl;
        if (g_weight_present) weightedErrorSum += (squareError * data.samp_weight[curP.second] * cur_weight);
        //ADDED THIS SINCE LAST MEETING WITH ANDREW
        if (fabs(last_weight - curP.first) > EPS) {
            cur_weight++;
        }
        //ADDED THIS
        last_weight = curP.first;
    }

    if (!g_weight_present) weightedErrorSum = squareErrorSum;

    frac.penalty = multi_const;
    frac.error = error_no_weight_at_all / subset.size();
    frac.fitness = weightedErrorSum / subset.size() * multi_const;

    return frac.fitness;

}


// Evaluation point

double ObjectiveFuncs::evaluate(MatrixContinuedFraction& frac, const data_store& dataset, const unordered_set<int>& subset, string objectiveFunc) {

    double multi_const;
    if (g_numFrac == "true") multi_const = 1.0 / (1.0 - double(frac.used_features()) / double(frac.vars));
    else multi_const = 1 + frac.used_features() * g_delta;

    double ret = 0;

    //dbg(fileSummary, "Evaluating with:"+objectiveFunc);

    if (objectiveFunc == "mse") ret = ObjectiveFuncs::mse(frac, dataset, subset, multi_const);
    else if (objectiveFunc == "mwe") ret = ObjectiveFuncs::mwe(frac, dataset, subset, multi_const);
    else if (objectiveFunc == "uwae") ret = ObjectiveFuncs::uwae(frac, dataset, subset, multi_const);
    else if (objectiveFunc == "abserr") ret = ObjectiveFuncs::abserr(frac, dataset, subset, multi_const);
    else if (objectiveFunc == "mase") ret = ObjectiveFuncs::mase(frac, dataset, subset);
    else if (objectiveFunc == "med.ae") ret = ObjectiveFuncs::medae(frac, dataset, subset);
    else if (objectiveFunc == "med.err") ret = ObjectiveFuncs::mederr(frac, dataset, subset);
    else if (objectiveFunc == "emc") ret = ObjectiveFuncs::emc(frac, dataset, subset);
    else if (objectiveFunc == "theil") ret = ObjectiveFuncs::theil(frac, dataset, subset);
    else if (objectiveFunc == "pearson") ret = ObjectiveFuncs::pearson(frac, dataset, subset, multi_const);
    else if (objectiveFunc == "nmse") ret = ObjectiveFuncs::nmse(frac, dataset, subset, multi_const);
    else if (objectiveFunc == "ub") ret = ObjectiveFuncs::ub(frac, dataset, subset, multi_const);
    else if (objectiveFunc == "lb") ret = ObjectiveFuncs::lb(frac, dataset, subset, multi_const);
    else if (objectiveFunc == "cw") ret = ObjectiveFuncs::cw(frac, dataset, subset, multi_const);
    else if (objectiveFunc == "cwy") ret = ObjectiveFuncs::cwy(frac, dataset, subset, multi_const);
    else if (objectiveFunc == "aic") ret = ObjectiveFuncs::aic(frac, dataset, subset);
    else if (objectiveFunc == "mre") ret = ObjectiveFuncs::mre(frac, dataset, subset, multi_const);
    // else if (objectiveFunc == "id") ret = ObjectiveFuncs::id(frac, dataset, subset, multi_const);
    else if (objectiveFunc == "mwrae") ret = ObjectiveFuncs::mwrae(frac, dataset, subset, multi_const);

    else dbg(fileSummary, "g_objec not defined");

    ret = isfinite(ret) ? ret : numeric_limits<double>::max() / (subset.size() + 1);
    return ret;

}