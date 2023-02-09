
/*
 * File:    data.hpp
 * Author: Haoyuan Sun
 *
 * * Description
 * - Train/test data ingestion & preprocess procedures
 * 
 * Improvement List
 * - Cache fraction evaluation
 * 
 * 210713, AC, Cleanup
 * 210000, MH, Improvements/additions
 * 180700, HS, Initial version 
 */
#include "cont_frac.hpp" 

#ifndef __DATA_STORE_HPP
#define __DATA_STORE_HPP

// Std Libs
#include <iosfwd>
#include <sstream>
#include <fstream>
#include <string>
#include <cassert>
#include <regex>
#include <numeric> 
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <iomanip>
#include <limits>
#include <type_traits>
// Local
#include "helper.hpp"
#include "globals.h"
#include "debug.hpp"


using namespace std;

class data_store {

    public: 

        vector<double> featMu;
        vector<double> featStd;
        double trgtMu, trgtStd;

        //min-max
        vector<double> featMin;
        vector<double> featMax;
        double trgtMin, trgtMax;

        //decimal scaling
        vector<double> maxValue;
        double trgtMaxVal;
        double trgtAvg;                 // Average of target value
        unsigned int num_entry;
        unsigned int num_var;
        unsigned int num_train = 0;
        unsigned int num_test = 0;
        bool has_weight = false;
        bool has_uncertainty = false;

        int normType;                   // Normalisation 1=min max, 2=z-score 3= Decimal scaling

        bool isSplit = false;
        bool isNormalised = false;

        vector<double> expected;        // target colum
        vector<double> samp_weight;     // weight column
        vector<double> dy;              // uncertainty column

        vector<size_t> sort_indexes(const vector<vector<double>> &v);
        
        vector<vector<double>> input;   // Data
        vector<string> names;           // Names of the variables (first row of data)
        vector<int> aic_param_count;    // AIC: parameter cout for feature/meta-feature

        vector<double> train_targets;
        vector<double> test_targets;
        vector<vector<double>> train_input;
        vector<vector<double>> test_input;
        double trnTrgtAvg, tstTrgtAvg;

        data_store();
        void clean();
        void assign(vector<vector<double>>, vector<double>);
        void assign(vector<vector<double>>, vector<double>, vector<string>);

        void debug();

        void read(const char* filename);
        void minmaxNormalise();
        void zScoreNormalise();
        void decimalScalingNormalise();

        void split(double trnSize);
        data_store filter(data_store &new_data, double pct);

        // Normalization helper function
        int getNormType();
        vector<double> getFeatMin();
        vector<double> getFeatMax();
        vector<double> getFeatMean();
        vector<double> getFeatStd();

        void reverseMinMax();
        void reverseZScore();
        

        // Misc functions from main

        static void printData(vector<vector<double>> data, vector<double> target);
        static void readAndNormData(data_store &ds, char* fileName, bool isTest);
        static void applyOnTestZScore(vector<vector<double>> &tstData, vector<double> &tstTarget, vector<double> fmu, vector<double> fstd, double trgtMu, double trgtStd);
        static void applyOnTestMinMax(vector<vector<double>> &tstData, vector<double> &tstTarget, vector<double> featMax, vector<double> featMin, double trgtMax, double trgtMin);

        void write_csv(string file, int run, bool isTrain);
        void write_csv_model(string file, int run, bool isTrain, MatrixContinuedFraction& frac);

};

#endif // __DATA_STORE_HPP

