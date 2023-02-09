
#ifdef __linux__
#include <signal.h>
#endif

using namespace std;


// Std Libs
#include <cmath>
#include <string>
#include <string>
#include <climits>
#include <iostream>
#include <cstdlib>
#include <chrono>
#include <fstream>
#include <ctime>
#include <cstdio>
// Local Libs
#include "data.hpp"
#include "globals.h"
#include "pop.hpp"
#include "helper.hpp"
#include "args.hpp"
#include "debug.hpp"

#define atoa(x) #x

// Default Global Variables (see "globals.h")
char*       g_trnFile = (char *)    "";
char*       g_tstFile = (char *)    "";
char*       g_solnFile = (char *)   "";
char*       g_paramCountFile=(char *)   "";
char*       g_logFile = (char *)    "Output/Run.log";
string      g_outPath =             "Output/";
float       g_delta =               0.35;
float       g_mu_rate =             0.2;
string      g_objec =               "mse";
string      g_numFrac =             "false";
bool        g_loocv =               false;
int         g_norm =                0;
int         g_numGen =              200;
int         g_current_run;
int         g_pop_init_method =     0;
bool        g_weight_present =      false;
size_t      g_depth =               4;
uint_fast32_t   g_seed =            0;
bool        g_run_serial =          true;
int         g_percnt_data =         100;
int         g_verbose_mode =        2;
int         g_local_search =        1;
int         g_reset_stuck =         5;
int         g_nm_runs =             4;          // Default as per the initial papers 2, 3, 4 on CFR
int         g_nm_max_iters =        250;        //  ""
int         g_nm_reset_stag =       10;         //  ""
double      g_window_size =         -1.0;
double      g_linreg_chance =       0;
int         g_PREC =                18;
char*       g_pop_init;
int         g_type;
int         g_base;
bool        g_origin;
bool        g_cont;
bool        g_aic_paramCount=       false;
int         g_int_terms =           0;
bool        g_int_only =            false;

// !!! Temporary vars commited quickly - need to make this more robust
bool        g_debug =               false;
int         g_depth_lock =          -1;                 
int         g_run_no =              1;                     
double      g_fitness =             numeric_limits<double>::max();
bool        g_prev_model_as_var =   false;  // true = if want to use previous model as meta feature
bool        g_incremental_depth =   false;  // Use the Incremental Depth Approach

// Other defaults
double      g_EPS = 1e-10;
double      g_ZERO = 0.000000000000000000;
double      MIN = numeric_limits<double>::min();
double      MAX = numeric_limits<double>::max();
randreal    g_rreal;
randint     g_rint;

const string VERSION =
    R"(VERSION:
    Continued Fraction Regression with Matrix Representation.
    Last update: 05-05-2021
    CFR version: 2.5.1
)";

int main(int argc, char* argv[]) {

    // Version header
    cout << VERSION << endl;

    // Read command line arguments
    Args::readCmdArguments(argc, argv);

    // Initialise random generator with device random
    g_rreal = randreal(g_seed); 
    g_rint = randint(g_seed);

    // Train/test datas
    data_store test_data, train_data;

    // Stamp start time
    auto start = chrono::system_clock::now();

    MatrixContinuedFraction sol;

    //NOT Leave-one-out Cross-validation
    if (g_loocv == false) {

        // read train data
        data_store::readAndNormData(train_data, g_trnFile, false);        

        //No test file supplied:: do 75-25 split
        if (strlen(g_tstFile) == 0) // split 75-25
        {
            if( g_verbose_mode > 0 )
                cout << "Data Separating into Train-Test: 75-25 split" << endl;

            train_data.split(0.75);
            test_data.assign(train_data.test_input, train_data.test_targets, train_data.names);
            train_data.assign(train_data.train_input, train_data.train_targets, train_data.names);
            
            if( g_verbose_mode > 1 ) {
                cout << "Train data" << endl;
                train_data.debug();
                cout << endl;
                cout << "Test data" << endl;
                test_data.debug();
            }
            
        } else {
            // read supplied Test data
            data_store::readAndNormData(test_data, g_tstFile, true);
        }

        //If Data Normalisation Selected
        // MinMax:: Apply normalisation on Test with the same scale of training data
        if (g_norm == 1) { // MinMax
            vector<double> fMin = train_data.getFeatMin();
            vector<double> fMax = train_data.getFeatMax();
            double tMin = train_data.trgtMin;
            double tMax = train_data.trgtMax;
            cout << "Target Min: " << tMin << "\tTarget Max: " << tMax << endl;
            data_store::applyOnTestMinMax(test_data.input, test_data.expected, fMax, fMin, tMax, tMin);
        } else if (g_norm == 2) { // z-Score
            vector<double> fMean = train_data.getFeatMean();
            vector<double> fStd = train_data.getFeatStd();
            double tmu = train_data.trgtMu;
            double tstd = train_data.trgtStd;
            data_store::applyOnTestZScore(test_data.input, test_data.expected, fMean, fStd, tmu, tstd);
        }

        /*
        *       Print Program Parameters before starting the Algorithm
        */
        Args::printProgramParameters();   
        
        /*
         * start the algorithm
         */

        start = chrono::system_clock::now();
        // Population Initialisation Method Check
        if (g_pop_init_method == 0) { //pop-init :: random (default)
           
            if( g_incremental_depth==false)  // for a Fixed Depth
            {
                if(g_objec=="aic")   // export the variable parameter count for AIC objective Fucntion 
                {
                    string soutParamCounter = g_outPath + "AICParamCounter"+to_string(g_depth)+".csv";
                    g_paramCountFile = strdup(soutParamCounter.c_str());
                    g_aic_paramCount = true;
                }
                population ant(train_data, test_data);
                ant.run(sol);
            }
            else{    // If Incremental Depth Approach is enabled
                // prepare out directory path from the log Directory Path
                g_outPath = g_logFile;
                int pos = g_outPath.find_last_of("/")+1;
                if (pos<0)
                    pos = g_outPath.find_last_of('\\')+1; // for Windows File System
                g_outPath.erase(pos);

                size_t iter_max_depth = g_depth;
                double best_fitness = MAX;
                for(size_t i = 0; i <= iter_max_depth; i++ ) {
                    g_depth = i;
                    // Prepare log file for each depth
                    // string soutFile = g_outPath +  "Run"+to_string(i)+".Sol.json";

                    if( i != 0 ) { 
                        // load previous model as the variable
                        if(g_prev_model_as_var == true)
                        {
                            if(g_objec=="aic") {
                                string soutParamCounter = g_outPath + "AICParamCounter"+to_string(i-1)+".csv";
                                g_paramCountFile = strdup(soutParamCounter.c_str());
                                g_aic_paramCount = true;
                            }

                            // read train data
                            train_data.clean();
                            string trainDataFile = g_outPath + "Run"+to_string(i-1)+".Train.csv";
                            char *dataTrn = strdup(trainDataFile.c_str());
                            data_store::readAndNormData(train_data, dataTrn, false);    

                            // read test data
                            test_data.clean();
                            string testDataFile = g_outPath + "Run"+to_string(i-1)+".Test.csv";
                            char *dataTst = strdup(testDataFile.c_str());
                            data_store::readAndNormData(test_data, dataTst, true);
                        }

                        // initialise the popluation with solution form previous depth
                        string prevSolnFile = g_outPath + "Run"+to_string(i-1)+".Sol.json";    
                        cout <<endl<< "Initialising Population with Previous Solution: " << prevSolnFile << endl;
                        population ant(train_data, test_data, prevSolnFile.c_str());
                        ant.run(sol);
                    }
                    else {
                        population ant(train_data, test_data);
                        ant.run(sol);
                    }
                    if(best_fitness > g_fitness) 
                    {   
                        best_fitness = g_fitness;
                        // only for Incremental Depth CFR
                        cout << endl << endl << "============== Depth Increased to " << g_depth+1 << " =================" << endl << endl;
                    }
                    else {
                        cout << endl << endl << "============== <Terminated> :: Fitness Stopped Improving at Depth=" << g_depth << " =================" << endl;
                        cout << "<Best Solution> found at Depth=" << g_depth-1 << endl;
                        break;
                    }
                    g_run_no++; // increament the global run counter for iterative execution
                }
            }

        }

    } else { // LOOCV option is selected

        
        data_store mydata;
        //read training data
        data_store::readAndNormData(mydata, g_trnFile, true);

        double avg_mse_trn = 0.0, avg_mse_tst = 0.0;
        char *fout = strcat(g_trnFile, ".LOOCV.csv");
        ofstream outfile(fout);
        cout << "Run\tMSE_Train\tMSE_Test" << endl;

        // Split the data into Train-Test using LOOCV method
        for (unsigned int i = 0; i < mydata.num_entry; i++) {
            vector<vector<double>> data = mydata.input;
            vector<double> targets = mydata.expected;

            //Create Test Data
            vector<vector<double>> tst_data;
            tst_data.push_back(data[i]);
            vector<double> tst_trgt;
            tst_trgt.push_back(mydata.expected[i]);
            test_data.assign(tst_data, tst_trgt);

            //Create train split of the Data
            targets.erase(targets.begin() + i);
            data.erase(data.begin() + i);
            train_data.assign(data, targets);
            train_data.num_entry = mydata.num_entry - 1;
            train_data.names = mydata.names;
            train_data.num_var = mydata.num_var;


            // Apply selected Normalisation Method on Training data
            if (g_norm == 1) {
                train_data.minmaxNormalise(); //Min-Max
            } else if (g_norm == 2) {
                train_data.zScoreNormalise(); //z-score
            } else if (g_norm == 3) {
                train_data.decimalScalingNormalise(); //decimnal scaling
            }

            //Apply normalisation on Test with the same scale of training data
            if (g_norm == 1) {
                vector<double> fMin = train_data.getFeatMin();
                vector<double> fMax = train_data.getFeatMax();
                double tMin = train_data.trgtMin;
                double tMax = train_data.trgtMax;
                data_store::applyOnTestMinMax(test_data.input, test_data.expected, fMax, fMin, tMax, tMin);

            } else if (g_norm == 2) {
                vector<double> fMean = train_data.getFeatMean();
                vector<double> fStd = train_data.getFeatStd();
                double tmu = train_data.trgtMu;
                double tstd = train_data.trgtStd;
                data_store::applyOnTestZScore(test_data.input, test_data.expected, fMean, fStd, tmu, tstd);
            }

            /*
             *  start the algorithm with LOOCV
             */

            if (g_pop_init_method == 0) { //pop-init :: random (default)
                population ant(train_data, test_data);
                ant.run(sol);
                pair <double, double> res_mse = ant.run(i);

                double trn_mse = res_mse.first;
                double tst_mse = res_mse.second;

                avg_mse_trn = avg_mse_trn + trn_mse;
                avg_mse_tst = avg_mse_tst + tst_mse;

                outfile << i << "\t" << trn_mse << "\t" << tst_mse << endl;
                cout << i << "\t" << trn_mse << "\t" << tst_mse << endl;
            } 
            
        }
        avg_mse_trn = avg_mse_trn / mydata.num_entry;
        avg_mse_tst = avg_mse_tst / mydata.num_entry;

        outfile << "a.LOOCV\t" << avg_mse_trn << "\t" << avg_mse_tst << endl;
        cout << "a.LOOCV\t" << avg_mse_trn << "\t" << avg_mse_tst << endl;
        outfile.close();
    }
    
    auto end = chrono::system_clock::now();
    chrono::duration<double> diff = end - start;


    if (g_verbose_mode > 0) {
        cout.precision(3);
        if (g_verbose_mode == 2) {
            cout << "Time in seconds: " << diff.count() << endl;
            cout << "Time in minutes: " << diff.count() / 60 << endl;
            cout << endl;
        } else
            cout << "Time in seconds: " << diff.count() << endl;
    }
    cout << "Done!" << endl;

    return 0;
}

