
#include "eval.hpp"

static double mse;
static double rmse;
static double nmse;
static double pred_avg;
static double meas_avg;
static double MYMAX;

evaluator::evaluator(const data_store& td) : data(td) {
    selection = unordered_set<int>();
    for(unsigned int i = 0; i < data.num_entry; i++)
        selection.insert(i);
    MYMAX = numeric_limits<double>::max() / (data.num_entry + 1); //moh: define new big number
}

evaluator::evaluator(const data_store& td, int num) : data(td) {
    selection = unordered_set<int>();
    unsigned int num_samples = (data.num_entry * num) / 100;

    MYMAX = numeric_limits<double>::max() / (data.num_entry + 1); //moh: define new big number

    // if too many being selected, might just as choose the whole
    //	if (num *2 >= data.num_entry) return;
    // moh:
    if (num_samples == (unsigned int) data.num_entry) {
        for(unsigned int i = 0; i < data.num_entry; i++)
            selection.insert(i);
        return;
    }

    //randint r_int;
    unsigned int iter = 0; // guard against super bad RNG
    //CCL
    if (!g_cont) {
        while (num_samples > 0 && iter++ < data.num_entry) {
            auto result = selection.insert(g_rint(0, data.num_entry - 1));
            if (result.second) --num_samples;
        }
    }
    else {
        unsigned int start_sample = g_rint(0, data.num_entry - 1 - num_samples);
        while (num_samples > 0 && iter++ < data.num_entry) {
            auto result = selection.insert(start_sample);
            start_sample++;
            if (result.second) --num_samples;
        }
    }

    // Sorting the int vector
    //sort(selection.begin(), selection.end());
}

evaluator::evaluator(const data_store& td, double window) : data(td) {
    /*
     * initialise the evaluator with user-defined size of the window
     * @param full training data and size of the window in percentage
     *
     */

    MYMAX = numeric_limits<double>::max() / (data.num_entry + 1); //moh: define new big number

    selection = unordered_set<int>();
    unsigned int num_samples = (data.num_entry * window);
    //    if(num_samples>=data.num_entry/2)   return;

    unsigned int start_pos = g_rint(0, (data.num_entry - num_samples)); //pick a windows starting point uniformly at random
    unsigned int end_pos = start_pos + num_samples;

    for (auto i = start_pos; i < end_pos; i++)
        selection.insert(i);
}

static inline double eval_pt(const data_store& td, MatrixContinuedFraction& frac, int i) {
    auto sqr = [] (double x) -> double {
        return x*x;
    };
    return sqr(td.expected[i] - frac.eval(td.input[i])) / td.num_entry;
}

static inline double process(double val) {
    return isfinite(val) ? val : MYMAX;
}

data_store evaluator::getData() {
    return data;
}

double evaluator::eval_fit(MatrixContinuedFraction& frac) const {
    /*
     * When we are evaluating only on a part of data
    */

    return ObjectiveFuncs::evaluate(frac, data, selection, g_objec);
}

double evaluator::eval_fit_full(MatrixContinuedFraction& frac) const {

    if(selection.size() != data.num_entry) {
        unordered_set<int> temp = unordered_set<int>();
        for(unsigned int i = 0; i < data.num_entry; i++) {
            temp.insert(i);
        }
        return ObjectiveFuncs::evaluate(frac, data, temp, g_objec);
    } else
        return eval_fit(frac);

}

void evaluator::print_val(const char* filename, MatrixContinuedFraction& frac) const {
    ofstream fout(filename);

    for (size_t i = 0; i < data.num_entry; ++i) {
        fout << data.expected[i] << ',';
        fout << frac.eval(data.input[i]) << '\n';
    }

    fout.close();
}

void evaluator::print_val_nd(const char* filename, MatrixContinuedFraction& frac, const data_store& trnData) const {
    ofstream fout(filename);

    //    cout<<"Test Dataset" <<endl;
    //    data.print();
    for (size_t i = 0; i < data.num_entry; ++i) {
        double nrmTrgt = (data.expected[i]*(trnData.trgtMax - trnData.trgtMin) + trnData.trgtMin);
        //        double nrmTrgt = data.expected[i];
        fout << nrmTrgt << ',';

        vector<double> nrm_data = data.input[i];
        for (size_t di = 0; di < nrm_data.size(); di++) {
            nrm_data[di] = (nrm_data[di]*(trnData.featMax[di] - trnData.featMin[di]) + trnData.featMin[di]);
        }
        //        double nrmPredict = frac.eval(data.input[i]);
        double nrmPredict = frac.eval(data.input[i])*((trnData.trgtMax - trnData.trgtMin) + trnData.trgtMin);
        fout << nrmPredict << '\n';

        cout << nrmTrgt << " " << nrmPredict << endl;
    }

    fout.close();
}

double evaluator::compute_MSET(MatrixContinuedFraction& frac, const data_store& trnData) const {

    /*	Mohammad added this function
    * http://lasa.epfl.ch/teaching/lectures/ML_MSc_Advanced/Practical/TP4b.pdf
    * 2.1 Regression Metrics
    * Mean Square Error (MSE) of the predictor
    * MSE = 1/n * sum_{i=1 to n} {(Yi-Y'i)^2}
    * Y' = predictor, Y=measured/observed
    */

    double sum_sq_err = 0.0, p_sum = 0.0, m_sum = 0.0;

    for (size_t i = 0; i < data.num_entry; ++i) {
        vector<double> nrm_data = data.input[i];
        for (size_t di = 0; di < nrm_data.size(); di++) {
            nrm_data[di] = (nrm_data[di]*(trnData.featMax[di] - trnData.featMin[di]) + trnData.featMin[di]);
        }
        double factor = (trnData.trgtMax - trnData.trgtMin);
        double pred_val = frac.eval(nrm_data) * factor + trnData.trgtMin; //reverse-normalize target

        double meas_val = data.expected[i];
        double err = (meas_val - pred_val);
        double sq_err = err * err; //residual^2
        sum_sq_err += sq_err;
        p_sum += pred_val;
        m_sum += meas_val;
    }
    mse = sum_sq_err / data.num_entry;
    pred_avg = p_sum / data.num_entry;
    meas_avg = m_sum / data.num_entry;
    //cout<< "Num Samples: " << data.num_entry << endl;
    return mse;
}

double evaluator::compute_MSE(MatrixContinuedFraction& frac) const {
    double sum_sq_err = 0.0, p_sum = 0.0, m_sum = 0.0;

    for (size_t i = 0; i < data.num_entry; i++) {
        // residual = Yi - Y'i
        double pred_val = frac.eval(data.input[i]);
        double meas_val = data.expected[i];
        double err = (meas_val - pred_val);
        double sq_err = err * err; //residual^2
        sum_sq_err += sq_err;
        p_sum += pred_val;
        m_sum += meas_val;
    }
    mse = sum_sq_err / data.num_entry;
    pred_avg = p_sum / data.num_entry;
    meas_avg = m_sum / data.num_entry;
    //    cout<< "Num Samples: " << data.num_entry << endl;
    return mse;
}

double evaluator::compute_ChiSquared(MatrixContinuedFraction& frac) const {
    
    // Perform Chi-Squared Error Metric with Uncertainity/Weight 
    // chi-squared = 1/n * \sum_{i=1}^{n} \frac{1}{sigma_i^2} (obs_i - pred_i)^2

    if(!data.has_uncertainty)   return -1;

    double squareErrorSum = 0.0,
            weightedErrorSum = 0.0;

    // for the subset of elements in data
    for (size_t i = 0; i < data.num_entry; i++) {

        double err = (data.expected[i] - frac.eval(data.input[i]));
        double squareError = err*err;
        squareErrorSum += (squareError);
        double norm_dy = 1/(data.dy[i] * data.dy[i]); // 1/sigma_i^2
        weightedErrorSum += (squareError * norm_dy);
    }

    return weightedErrorSum / (data.num_entry);
    
}

double evaluator::compute_ChiSquaredPct(MatrixContinuedFraction& frac) const {
    
    // Perform Chi-Squared Error Metric with Uncertainity/Weight 
    // chi-squared = 1/n * \sum_{i=1}^{n} \frac{1}{sigma_i^2} (obs_i - pred_i)^2

    if(!data.has_uncertainty)   return -1;

    int good = 0.0,
        bad = 0.0;

    // for the subset of elements in data
    for (size_t i = 0; i < data.num_entry; i++) {

        double err = (data.expected[i] - frac.eval(data.input[i]));
        if (abs(err) <= data.dy[i] )
            good += 1;
        else
            bad += 1;

    }

    double pct = (good/(double)data.num_entry)*100;
    return pct;
}

double evaluator::compute_RMSE(MatrixContinuedFraction& frac) const {
    rmse = sqrt(mse);
    return rmse;
}

double evaluator::compute_PEARSON(MatrixContinuedFraction& frac) const {
    double pearson = 0.0, sum_num = 0.0, meas_var = 0.0, pred_var = 0.0;
    //            double pearson_type = -1; //-1 to prefer positive coefficient
    for (size_t i = 0; i < data.num_entry; i++) {
        // residual = Yi - Y'i
        double pred_val = frac.eval(data.input[i]);
        double meas_val = data.expected[i];
        double meas_diff = (meas_val - meas_avg);
        double pred_diff = (pred_val - pred_avg);

        sum_num += pred_diff*meas_diff;
        meas_var += meas_diff*meas_diff;
        pred_var += pred_diff*pred_diff;
    }
    pearson = sum_num / pow(meas_var*pred_var, 0.5);
    return pearson;
}

double evaluator::compute_NMSE(MatrixContinuedFraction& frac) const {

    /*
     * Mohammad added this function
     * http://lasa.epfl.ch/teaching/lectures/ML_MSc_Advanced/Practical/TP4b.pdf
     * 2.1 Regression Metrics
    * The Normalized Mean Square Error (NMSE) is simply the MSE normalized by the
    * variance of the observed values:
    *  NMSE = MSE / VAR(Y)
    *  where VAR(Y) = 1/(n-1) * sum_{i=1 to n} {(Yi-mu)^2} and mu = avg(Y)
    */

    double m_sum_var = 0.0;

    for (size_t i = 0; i < data.num_entry; i++) {
        double var_i = (data.expected[i] - meas_avg);
        m_sum_var += (var_i * var_i); //residual^2
    }
    double m_variance = m_sum_var / (data.num_entry - 1);
    nmse = mse / m_variance;
    return nmse;
}

void evaluator::compute_Scores(MatrixContinuedFraction& frac) const {

    double mse, nmse, pearson, theil, med_err, med_ae, mase;
    //    double sum_num = 0.0, meas_var = 0.0, pred_var = 0.0;
    vector<double> errors; //med.err
    vector<double> abs_errors;
    double sum_forecast_err = 0.0; //mase

    /* --------- MSE --------- */
    double sum_sq_err = 0.0, p_sum = 0.0, m_sum = 0.0;
    for (size_t i = 0; i < data.num_entry; i++) {
        // residual = Yi - Y'i
        double pred_val = frac.eval(data.input[i]);
        double meas_val = data.expected[i];
        double err = (meas_val - pred_val);
        double sq_err = err * err; //residual^2
        sum_sq_err += sq_err;
        p_sum += pred_val;
        m_sum += meas_val;

        //med.err
        if (meas_val < pred_val)
            err = (pred_val - meas_val);
        errors.push_back(err);

        //med.ae
        double abs_err = g_ZERO;
        if (meas_val < g_ZERO && pred_val < g_ZERO)
            abs_err = fabs(fabs(meas_val) - fabs(pred_val));
        else if (meas_val < g_ZERO && pred_val > g_ZERO)
            abs_err = fabs(meas_val) + pred_val;
        else if (meas_val > g_ZERO && pred_val < g_ZERO)
            abs_err = meas_val + fabs(pred_val);
        else abs_err = fabs(meas_val - pred_val);

        abs_errors.push_back(abs_err);
        sum_forecast_err += abs_err;
    }
    mse = sum_sq_err / data.num_entry;
    pred_avg = p_sum / data.num_entry;
    meas_avg = m_sum / data.num_entry;

    /* --------- PEARSON --------- */
    double m_sum_var = 0.0; // NMSE
    double sum_num = 0.0, meas_var = 0.0, pred_var = 0.0;

    for (size_t i = 0; i < data.num_entry; i++) {
        // residual = Yi - Y'i
        double pred_val = frac.eval(data.input[i]);
        double meas_val = data.expected[i];
        double meas_diff = (meas_val - meas_avg); // variance
        double pred_diff = (pred_val - pred_avg);

        sum_num += pred_diff*meas_diff;
        meas_var += meas_diff*meas_diff;
        pred_var += pred_diff*pred_diff;

        //NMSE
        m_sum_var += (meas_diff * meas_diff); //residual^2
    }
    pearson = sum_num / pow(meas_var*pred_var, 0.5);

    /* --------- NMSE --------- */
    double m_variance = m_sum_var / (data.num_entry - 1);
    nmse = mse / m_variance;

    /* --------- MED.ERR --------- */
    sort(errors.begin(), errors.end());
    if (errors.size() % 2 == 0)
        med_err = (errors[errors.size() / 2 - 1] + errors[errors.size() / 2]) / 2;
    else
        med_err = errors[errors.size() / 2];

    /* --------- MED.AE --------- */
    sort(abs_errors.begin(), abs_errors.end());
    if (abs_errors.size() % 2 == 0)
        med_ae = (abs_errors[abs_errors.size() / 2 - 1] + abs_errors[abs_errors.size() / 2]) / 2;
    else
        med_ae = abs_errors[abs_errors.size() / 2];

    /* --------- MASE --------- */
    // mean absolute scaled error
    //compute value of numerator
    double nom = sum_forecast_err / data.num_entry;

    //compute value of denominator
    double sum_ae = 0.0;
    for (long t = 1; t < data.num_entry; ++t) {
        double Y_t = data.expected[t];
        double Y_prev = data.expected[t - 1];
        double abs_err = g_ZERO;
        if (Y_t < g_ZERO && Y_prev < g_ZERO)
            abs_err = fabs(fabs(Y_t) - fabs(Y_prev));
        else if (Y_t < g_ZERO && Y_prev > g_ZERO)
            abs_err = fabs(Y_t) + Y_prev;
        else if (Y_t > g_ZERO && Y_prev < g_ZERO)
            abs_err = Y_t + fabs(Y_prev);
        else abs_err = fabs(Y_t - Y_prev);

        sum_ae += abs_err;
    }
    double denom = sum_ae / (data.num_entry - 1);

    // compute MASE
    mase = nom / denom;


    /* --------- THEIL --------- */
    vector<double> slopes;
    for (size_t i = 0; i < data.num_entry; ++i) {
        for (size_t j = 0; j < data.num_entry; ++j) {
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
    if(slopes.size() > 0 ) {
        size_t size = slopes.size();
        sort(slopes.begin(), slopes.end());
        if (size % 2 == 0) {
            theil = (slopes[size / 2 - 1] + slopes[size / 2]) / 2;
        } else {
            theil = slopes[size / 2];
        }
    }


    // size_t pad = g_PREC+5;

    cout.precision(g_PREC);
    cout << mse 
         << "\t" << nmse 
         << "\t" << pearson 
         << "\t" << theil 
         << "\t" << med_err 
         << "\t" << med_ae 
         << "\t" << mase;

    if( data.has_uncertainty ) {
        cout << "\t" << compute_ChiSquared(frac);
        cout << "\t" << compute_ChiSquaredPct(frac) << "%";
    } else {
        cout << "\t" << "N/A";
        cout << "\t" << "N/A";
    }
        
    cout << endl;

}
