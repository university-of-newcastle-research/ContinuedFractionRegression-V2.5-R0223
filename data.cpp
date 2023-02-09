
#include "data.hpp"

data_store::data_store() {
    /*
     * Construct datastore with default values and empty array
     */

    num_entry = 0;
    num_var = 0;
    num_train = 0;
    num_test = 0;

    isSplit = false;
    isNormalised = false;

    expected = vector<double>(); //target
    input = vector<vector<double>>(); //data
    names = vector<string>();
    aic_param_count = vector<int>(num_var+1, 0);    
}

void data_store::clean() {
    /*
     * Clean the datastore with default values and empty array
     */

    num_entry = 0;
    num_var = 0;
    num_train = 0;
    num_test = 0;

    isSplit = false;
    isNormalised = false;

    expected = vector<double>(); //target
    input = vector<vector<double>>(); //data
    names = vector<string>();
    aic_param_count = vector<int>(num_var+1, 0);    
}

void data_store::assign(vector<vector<double>> data, vector<double> target) {
    /*
     * 
     */

    num_entry = data.size();
    num_var = data[0].size();
    num_train = 0;
    num_test = 0;

    isSplit = false;
    isNormalised = false;

    expected = target; //target
    input = data; //data
    names = vector<string>();
    aic_param_count = vector<int>(num_var+1, 0);  
}

void data_store::assign(vector<vector<double>> data, vector<double> target, vector<string> var_names) {
    /*
     * 
     */

    num_entry = data.size();
    num_var = data[0].size();
    num_train = 0;
    num_test = 0;

    isSplit = false;
    isNormalised = false;

    expected = target; //target
    input = data; //data
    names = var_names;
    aic_param_count = vector<int>(num_var+1, 0);
}

// Pre-process functions

void data_store::split(double trnSize) {
    /*
     * 
     */

    num_train = (long) (trnSize * num_entry);
    num_test = num_entry - num_train;

    // initialise a vector with the sequential index of samples
    vector<long> indexVector;
    back_insert_iterator<vector<long> > p =
            back_inserter(indexVector);

    for (unsigned int i = 0; i < num_entry; i++)
        *p = i;

    // random shuffle those indices
    random_shuffle(indexVector.begin(), indexVector.end());

    //    cout << "Start Train Split: " ;
    // split first trnSize percentage of samples into train set
    vector<long>::iterator it = indexVector.begin();
    double sum = 0.0;
    for (unsigned int count = 0; count < num_train; ++it, count++) {
        unsigned int index = *it;
        //        cout << index << " ";
        train_targets.push_back(expected.at(index));
        sum += expected.at(index);
        vector<double> data_at_index = input.at(index);
        train_input.push_back(data_at_index);
    }
    trnTrgtAvg = sum / num_train;


    //remaining data should be placed in test set
    sum = 0.0;
    while (it != indexVector.end()) {
        unsigned int index = *it;
        //        cout << index << " ";

        test_targets.push_back(expected.at(index));
        sum += expected.at(index);
        vector<double> data_at_index = input.at(index);
        test_input.push_back(data_at_index);
        it++;
    }
    tstTrgtAvg = sum / (num_entry - num_train);

    isSplit = true;
}

vector<size_t> data_store::sort_indexes(const vector<vector<double>> &v) {

    // initialize original index locations
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values 
    stable_sort(idx.begin(), idx.end(),
        [&v](size_t i1, size_t i2) {return abs(v[i1][0]) < abs(v[i2][0]);});

    return idx;

}

data_store data_store::filter(data_store &new_data, double pct) {
    /*
     * 
     */   

    new_data.names = names;
    new_data.num_var = input[0].size();
    
    new_data.num_train = (long) (pct * num_entry);
    new_data.num_test = new_data.num_train;
    new_data.num_entry = new_data.num_train;

    new_data.has_uncertainty = has_uncertainty;
    new_data.has_weight = has_weight;
    new_data.isSplit = false;
    new_data.isNormalised = false;

    vector<size_t> idx = sort_indexes(input);

    size_t count = 0;
    for(int i : idx) {

        new_data.input.push_back(input[i]);
        new_data.expected.push_back(expected[i]);

        if( has_uncertainty )
            new_data.dy.push_back(dy[i]);
        
        if( has_weight )
            new_data.samp_weight.push_back(samp_weight[i]);

        count++;

        if( count == new_data.num_train )
            break;

    }

    return new_data;
}

void data_store::decimalScalingNormalise() {
    /*
     * Normalise the data using decimal scaling
     */

    maxValue.resize(num_var);

    // Find max values  
    for (size_t i = 0; i < input.size(); i++) {
        for (size_t j = 0; j < input[0].size(); j++) {
            maxValue[j] = max(maxValue[j], input[i][j]);
        }
    }
  
    // Max of target
    trgtMaxVal = *max_element(expected.begin(), expected.end());

    // Normalise each variable/feature based on max value
    for (size_t i = 0; i < input.size(); i++) {
        for (size_t j = 0; j < input[0].size(); j++) {
            if (Helper::almost_equal(maxValue[j], g_ZERO, g_PREC))
                input[i][j] = g_ZERO;
            else {
                int logVal = (int) log10(maxValue[j]) + 1;
                input[i][j] = (input[i][j] / (pow(10, logVal)));
            }
        }
    }
    
    isNormalised = true;
    normType = 3;
}

void data_store::minmaxNormalise() {
    /*
     * Normalise the data using Min-Max
     */ 
    
    // Set size of vectors
    featMin.resize(num_var);
    featMax.resize(num_var);

    // Find min and max values
    for (size_t i = 0; i < input.size(); i++) {
        for (size_t j = 0; j < input[0].size(); j++) {
            featMin[j] = min(featMin[j], input[i][j]);
            featMax[j] = max(featMax[j], input[i][j]);
        }
    }

    // Min and max of target
    trgtMin = *min_element(expected.begin(), expected.end());
    trgtMax = *max_element(expected.begin(), expected.end());

    // Normalise the matrix values: min-max normalisation
    for (size_t i = 0; i < input.size(); i++) {
        for (size_t j = 0; j < input[0].size(); j++) {
            double denom = (featMax[j] - featMin[j]);
            if (Helper::almost_equal(denom, g_ZERO, g_PREC))
                input[i][j] = g_ZERO;
            else
                input[i][j] = ((input[i][j] - featMin[j])*1) / denom;
        }
    }

    // Normalize target
    double denom = (trgtMax - trgtMin);
    for (size_t j = 0; j < expected.size(); j++) {
        if (Helper::almost_equal(denom, g_ZERO, g_PREC))
            expected[j] = g_ZERO;
        else
            expected[j] = (expected[j] - trgtMin) / denom;
    }

    isNormalised = true;
    normType = 1;
}

void data_store::reverseMinMax() {
    /*
     * Reverse normalised matrix values: min-max de-normalisation
     */

    for (size_t i = 0; i < input.size(); i++) {
        for (size_t j = 0; j < input[0].size(); j++) {
            input[i][j] = (input[i][j] * (featMax[j] - featMin[j]) + featMin[j]);
        }
    }

    // De-normalize target
    double factor = (trgtMax - trgtMin);
    for (size_t j = 0; j < expected.size(); j++) {
        expected[j] = (expected[j] * factor + trgtMin);
    }
}

void data_store::reverseZScore() {
    /*
     * Reverse normalised matrix values: z-score de-normalisation
     */

    for (size_t i = 0; i < input.size(); i++) {
        for (size_t j = 0; j < input[0].size(); j++) {
            input[i][j] = ((input[i][j] * featStd[j]) + featMu[j]);
        }
    }

    // De-normalize target
    for (size_t j = 0; j < expected.size(); j++) {
        expected[j] = ((expected[j] * trgtStd) + trgtMu);
    }
}

void data_store::zScoreNormalise() {
    /*
     * z-score normalisation of data
     */

    featMu.resize(num_var);
    featStd.resize(num_var);

    // Transpose the input matrix so that we can have the values per features in a single vector
    vector<vector<double>> dataTrans(input[0].size(), vector<double>(input.size()));
    for (size_t i = 0; i < input.size(); i++)
        for (size_t j = 0; j < input[0].size(); j++)
            dataTrans[j][i] = input[i][j];

    // Find the standard deviation for each features/variables    
    for (size_t i = 0; i < dataTrans.size(); i++) {
        
        vector<double> v = dataTrans[i];

        // Compute mean score of the feature/variable
        double sum = accumulate(begin(v), end(v), 0.0);
        double mean = sum / v.size();
        featMu[i] = mean; // save mean of feature

        // z-score uses "Sample Standard Deviation" sx = sqrt( sum((x-mean)^2)/(n-1)) )
        // "sample standard deviation" formula is different than the "population standard deviation"
        double accum = 0.0;
        for_each(begin(v), end(v), [&](const double d) {
            accum += (d - mean) * (d - mean);
        });
        double stdev = sqrt(accum / (v.size() - 1));
        featStd[i] = stdev; // save std of feature score

        // Put 0 as z-score when standard dev is 0
        if (stdev == 0.0) {
            fill(dataTrans[i].begin(), dataTrans[i].end(), 0.0);
        } else {
            transform(v.begin(), v.end(), dataTrans[i].begin(), [mean, stdev](double x) {
                return (x - mean) / stdev; });
        }

    }

    // Put the normalized data back into input
    for (size_t i = 0; i < input.size(); i++) {
        for (size_t j = 0; j < input[0].size(); j++) {
            input[i][j] = dataTrans[j][i];
        }
    }

    // Normalize target
    double sum = accumulate(begin(expected), end(expected), 0.0);
    double mean = sum / expected.size();
    trgtMu = mean;      // save mean of target

    double accum = 0.0;
    for_each(begin(expected), end(expected), [&](const double d) {
        accum += (d - mean) * (d - mean);
    });
    double stdev = sqrt(accum / (expected.size() - 1));
    trgtStd = stdev;    // save std of target

    // Put 0 as z-score when standard dev is 0
    if (stdev == 0.0) {
        fill(expected.begin(), expected.end(), 0.0);
    } else {
        transform(expected.begin(), expected.end(), expected.begin(), [mean, stdev](double x) {
            return (x - mean) / stdev; });
    }

    isNormalised = true;
    normType = 2;
}

// Accessor function

int data_store::getNormType() {
    return normType;
}

vector<double> data_store::getFeatMin() {
    return featMin;
}

vector<double> data_store::getFeatMax() {
    return featMax;
}

vector<double> data_store::getFeatMean() {
    return featMu;
}

vector<double> data_store::getFeatStd() {
    return featStd;
}

// Helper functions

bool checkSuffix(const string& s, const string& suffix) {
    return s.rfind(suffix) == (s.size() - suffix.size());
}

void data_store::read(const char* filename) {
    if( g_verbose_mode > 0 )    cout << "Data Loading From CSV: " << filename << endl;

    string targetHeader = "y";
    string weightHeader = "w";
    string uncertantyHeader = "dy";

     if (checkSuffix(filename,".csv")) {

        // Read in CSV file
        ifstream fin(filename);
        
        if (!fin.is_open()) {
            cout << "Unable to read file: " << filename << endl;
            exit(1);
        }

        bool isFirst = true;
        if (fin.is_open()) {
            
            int weightCol = -1;
            int dyCol = -1;

            regex words_regex("[^,]+");
            
            // Read each line
            string line;
            while ( getline(fin, line) ) {

                // Process header row                
                if (isFirst) {

                    // Setup regex
                    auto words_begin = sregex_iterator(line.begin(), line.end(), words_regex);
                    auto words_end = sregex_iterator();

                    // For each variable
                    num_var = 0;
                    size_t colNo = 0;
                    for (sregex_iterator i = words_begin; i != words_end; ++i) {

                        // Extract String and replace end of line on header
                        string text = Helper::trim((*i).str());
                        size_t index = text.find("\r", 0);
                        if (index != string::npos)
                            text.replace(index, 1, "");                    
        
                        if ( colNo == 0 ) {
                            if( g_verbose_mode > 0 )        cout << "Target forced at " << colNo << endl;
                            names.push_back(text);
                        }
                        // Exclude weight column
                        else if ( text.compare(weightHeader) == 0 ) {

                            if ( weightCol != -1 ) {
                                cout << "Data has TWO COLUMNS with w as header" << endl;
                                exit(1);
                            }

                            weightCol = colNo;
                            has_weight = true;
                            if( g_verbose_mode > 0 )    cout << "Weight detected at " << colNo << endl;
                        }
                        // Exclude uncertainy column
                        else if ( text.compare(uncertantyHeader) == 0 ) {

                            if ( dyCol != -1 ) {
                                cout << "Data has TWO COLUMNS with dy as header" << endl;
                                exit(1);
                            }

                            dyCol = colNo;
                            has_uncertainty = true;
                            if( g_verbose_mode > 0 ) cout << "Uncertainty detected at " << colNo << endl;
                        } 
                        // Process all other cols
                        else {
                            names.push_back(text); 
                            num_var += 1;
                            if( g_verbose_mode > 0 ) cout << "IV detected " << text << " at " << colNo << endl;
                        }
                        colNo += 1;

                    }

                    if( g_verbose_mode > 0 && dyCol == -1 )       cout << "No uncertainty column detected" << endl;
                    if( g_verbose_mode > 0 && weightCol == -1 )   cout << "No weight column detected" << endl;

                    if( dyCol == -1 && (g_objec == "cw" || g_objec == "cwy" || g_objec == "uwae") ) {
                        cout << "No uncertainty when required for objective function " << g_objec << endl;
                        cout << "Exiting" << endl;
                        exit(1);
                    }
                        

                    isFirst = false;

                } else {

                    // Setup regex
                    vector<string> split_line;
                    auto words_begin = sregex_iterator(line.begin(), line.end(), words_regex);
                    auto words_end = sregex_iterator();

                    // Extract row
                    for (sregex_iterator i = words_begin; i != words_end; ++i)
                        split_line.push_back((*i).str());                           
                    
                    // Increment row count
                    num_entry += 1;
                    for (size_t j = 0; j < split_line.size(); j++) {
                                            
                        double val = stod(split_line[j]);
                        
                        if ( j == 0 ) {
                            expected.push_back(val);
                        }
                        else if ( (int)j == weightCol) {
                            samp_weight.push_back(val);  
                        }
                        else if ( (int)j == dyCol) {
                            dy.push_back(val);
                        } else {
                            // if first IV in sample, create a new vector
                            if( j == 1)
                                input.push_back(vector<double>());
                            input.back().push_back(val);
                        }

                    }

                    // Handle no weight col
                    if (weightCol == -1)
                        samp_weight.push_back(1.0); 
                }
            }
            fin.close();
        }

        // read file if AIC parameter count is present
        aic_param_count = vector<int>(num_var+1, 0);  

        if( g_aic_paramCount==true){
            // Read in CSV file
            ifstream finp(g_paramCountFile);
            
            if (!finp.is_open()) {
                cout << "Unable to read file: " << g_paramCountFile << endl;
                exit(1);
            }

            regex words_regex("[^,]+");

            // Read each line
            string pline;
            while ( getline(finp, pline) ) {
                vector<string> split_line;
                auto words_begin = sregex_iterator(pline.begin(), pline.end(), words_regex);
                auto words_end = sregex_iterator();
                // Extract row
                for (sregex_iterator i = words_begin; i != words_end; ++i)
                    split_line.push_back((*i).str());                           
                    
                // Read Value
                for (size_t j = 0; j < split_line.size(); j++) {
                                        
                    int val = stoi(split_line[j]);
                    aic_param_count[j] = val;
                }  
                break; 
            }
        }

    } else if (checkSuffix(filename,".in")) {
        cout << ".in Filetype no longer used. Please use CFR" << endl;
        exit(1);
    } else {
        cout << "Invalid file type specified" << endl;
        exit(1);
    }

    if( g_verbose_mode > 1 ) {
        debug();
        cout << endl;
    }   
    

}

void data_store::debug() {

    /*
     * Debug data store class
     */

    cout << "Data debug samples:" << num_entry << " feats:" << num_var << " weight?:" << has_weight << " uncertainty?:" << has_uncertainty << endl;

    for (unsigned int i = 0; i < num_entry; i++) {

        cout << "y:" << expected[i] << ", ";
        for (unsigned int j = 0; j < num_var; j++) {

            cout << names[j+1] << ": " << input[i][j] << " ";
        }

        if( has_weight )
            cout << ", w: " << samp_weight[i];
        
        if( has_uncertainty )
            cout << ", dy: " << dy[i];

        cout << endl;
    }
    cout << endl;

}

void data_store::readAndNormData(data_store &ds, char* fileName, bool isTest) {
    /*
    * Function to read data and apply chosen Normalisation method only on training data
    * @param data-storage to save, data file name, indicate if test file
    */
    ds.read(fileName);
    if (isTest == false) { //apply normalisation if it is Training data
        if (g_norm == 1) {
            ds.minmaxNormalise(); //normalise train data
            //        ds.reverseMinMax(); // de-normlize min-max
            cout << "Min Max Normalisation" << endl;
        } else if (g_norm == 2) {
            ds.zScoreNormalise(); //z-score normalise train data
            //        ds.reverseZScore();
            cout << "z-Score Normalisation" << endl;
        } else if (g_norm == 3) {
            ds.decimalScalingNormalise(); //decimnal scaling normalise train data
            cout << "Decimal Scaling Normalisation" << endl;
        }
    }
    for(unsigned int i = 0; i < ds.num_entry; i++) {
        dbg(fileSolution,to_string(ds.expected[i])+","+to_string(ds.input[i][0]),1,0);
    }
}

void data_store::printData(vector<vector<double>> data, vector<double> target) {
    /*
     * Show the data as <sample> : <target>
     * @param sample and target
     */

    for (size_t i = 0; i < data.size(); ++i) {
        for (size_t j = 0; j < data[0].size(); ++j) {
            cout << data[i][j] << " ";
        }
        cout << ": " << target[i] << endl;
    }
}

void data_store::applyOnTestZScore(vector<vector<double>> &tstData, vector<double> &tstTarget, vector<double> fmu, vector<double> fstd, double trgtMu, double trgtStd) {
    /*
    * Apply Z-Score Normalisation on the Test Set
    * @param data <sample><target> and required parameters for the z-score computed on training data
    */

    //	we can have the values per features in a single vector
    vector<vector<double>> dataTrans(tstData[0].size(),
            vector<double>(tstData.size()));
    for (size_t i = 0; i < tstData.size(); i++)
        for (size_t j = 0; j < tstData[0].size(); j++)
            dataTrans[j][i] = tstData[i][j];

    for (size_t i = 0; i < dataTrans.size(); i++) {
        vector<double> v = dataTrans[i];
        double mean = fmu[i]; // save mean of feature
        double stdev = fstd[i]; // save std of feature score

        // put 0 as z-score when standard dev is 0
        if (stdev == 0.0) {
            fill(dataTrans[i].begin(), dataTrans[i].end(), 0.0);
        } else {
            transform(v.begin(), v.end(), dataTrans[i].begin(), [mean, stdev](double x) {
                return (x - mean) / stdev; });
        }
    }

    // put the normalized data back into input
    for (size_t i = 0; i < tstData.size(); i++) {
        for (size_t j = 0; j < tstData[0].size(); j++) {
            tstData[i][j] = dataTrans[j][i];
        }
    }

    // put 0 as z-score when standard dev is 0
    if (trgtStd == 0.0) {
        fill(tstTarget.begin(), tstTarget.end(), 0.0);
    } else {
        transform(tstTarget.begin(), tstTarget.end(), tstTarget.begin(), [trgtMu, trgtStd](double x) {
            return (x - trgtMu) / trgtStd; });
    }
}

void data_store::applyOnTestMinMax(vector<vector<double>> &tstData, vector<double> &tstTarget, vector<double> featMax, vector<double> featMin, double trgtMax, double trgtMin) {
    /*
     * Apply MinMax Normalisation on the Test Set
     * @param data <sample><target> and required parameters for the MinMax computed on training data
     */

    // normalise the matrix values: min-max normalisation
    for (size_t i = 0; i < tstData.size(); i++) {
        for (size_t j = 0; j < tstData[0].size(); j++) {
            double denom = (featMax[j] - featMin[j]);
            if (Helper::almost_equal(denom, g_ZERO, g_PREC))
                tstData[i][j] = g_ZERO;
            else
                tstData[i][j] = ((tstData[i][j] - featMin[j])*1) / denom;
        }
    }

    //normalize target
    double denom = (trgtMax - trgtMin);
    for (size_t j = 0; j < tstTarget.size(); j++) {
        if (Helper::almost_equal(denom, g_ZERO, g_PREC))
            tstTarget[j] = g_ZERO;
        else
            tstTarget[j] = (tstTarget[j] - trgtMin) / denom;
    }
}

void data_store::write_csv(string file, int run, bool isTrain) {

    // We ignore file at the moment
    string soutFile = g_logFile;

    // Setup name
    if( isTrain )
        soutFile.replace(soutFile.length()-7,7,"Output/Run"+to_string(run)+".Train.csv");
    else 
        soutFile.replace(soutFile.length()-7,7,"Output/Run"+to_string(run)+".Test.csv");
    
    // Open File
    ofstream fout(soutFile);
    if (!fout.is_open()) {
        cout << "Unable to open file " << soutFile << endl;
        return ;
    }

    // Output header and all rows
    for (unsigned int i = 0; i < num_entry; i++) {

        // Output header
        if( i == 0 ) {

            fout << "y";
            for (unsigned int j = 0; j < num_var; j++) {
                fout << "," << names[j+1];
            }

            if( has_weight )
                fout << "," << samp_weight[i];
    
            if( has_uncertainty )
                fout << "," << dy[i];

            fout << endl;
        }
        
        // Output rows
        fout << expected[i];                                // y
        for (unsigned int j = 0; j < num_var; j++) {        // all vars
            double val = input[i][j];
            fout << "," << val;
        }

        if( has_weight )
            fout << "," << samp_weight[i];
        
        if( has_uncertainty )
            fout << "," << dy[i];

        fout << endl;

    }
}

void data_store::write_csv_model(string file, int run, bool isTrain, MatrixContinuedFraction& frac) {

    // We ignore file at the moment
    string soutFile = g_logFile;

    // Setup name
    if( isTrain )
        soutFile = g_outPath + "Run"+to_string(run)+".Train.csv";  //Creating the Train output file as per execution number
    else
        soutFile = g_outPath + "Run"+to_string(run)+".Test.csv";  //Creating the Test output file as per execution number   

    if(g_objec=="aic"){
        string soutParamCounter = g_outPath + "AICParamCounter"+to_string(run)+".csv";
        ofstream fpout(soutParamCounter);
        if (!fpout.is_open()) {
            cout << "Unable to open file " << soutParamCounter << endl;
            return ;
        }
        // Output header and all rows
        fpout << aic_param_count[0];
        for (unsigned int i = 1; i < num_var; i++) {
            fpout << "," << aic_param_count[i];
        }
        if( has_weight )
            fpout << ",0";

        if( has_uncertainty )
            fpout << ",0";

        fpout << "," << frac.freeParams << endl;
        // fpout.close();
    }

    // Open File
    ofstream fout(soutFile);
    if (!fout.is_open()) {
        cout << "Unable to open file " << soutFile << endl;
        return ;
    }

    // Output header and all rows
    for (unsigned int i = 0; i < num_entry; i++) {

        // Output header
        if( i == 0 ) {

            fout << "y";
            for (unsigned int j = 0; j < num_var; j++) {
                fout << "," << names[j+1];
            }

            if( has_weight )
                fout << "," << samp_weight[i];
    
            if( has_uncertainty )
                fout << "," << dy[i];

            fout << ",m" << run << endl;
        }
        
        // Output rows
        fout << expected[i];                                // y
        for (unsigned int j = 0; j < num_var; j++) {        // all vars
            double val = input[i][j];
            fout << "," << val;
        }

        if( has_weight )
            fout << "," << samp_weight[i];
        
        if( has_uncertainty )
            fout << "," << dy[i];

        fout << "," << frac.eval(input[i]) << endl;

    }
}
