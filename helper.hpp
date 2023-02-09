
#ifndef __HELPER_HPP
#define __HELPER_HPP

/*
 * File:   helper.hpp
 * Author: aciezak
 *
 * Description
 * - static helper functions
 * 
 * Improvement List
 * -
 * 
 * 210830, AC, Document, introduce safe math operations
 * 210000, AC, Initial revision
 */

// Std Lib
#include <cmath>
#include <limits>
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <regex>
#include <stdexcept>
// Local
#include "globals.h"

using namespace std;

class Helper {

    public:

        // File operations
        static bool     is_file_exist(const char *fileName);

        // Safe mathematical operations
        static double   times(double, double);
        static double   divide(double, double);
        static double   add(double, double);
        static double   exp(double);
        static double   pow(double, double);
        static bool     almost_equal(double i, double j, int ulp);
        static bool     is_between(double value, double min, double max);

        static long long    gcd(long long a,  long long b);         // Euclidean Algorithm to find the GCD of a and b
        static void         decimalToRational(double coeff, long long &numer, long long &denom); // Function to convert decimal to Rational

        // String helpers 
        static vector<string> split(string s, string delimiter);
        static string  ltrim(string s);
        static string  rtrim(string s);
        static string  trim(string s);
        
};
#endif
