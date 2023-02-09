
/*
 * File:    args.hpp
 * Authors: aciezak
 * 
 * Description
 * - Encapsulate argument processing logic
 * 
 * Improvements
 * -
 *
 * 210000, AC, initial revision
 * 
 */

#ifndef __ARGS_HPP
#define __ARGS_HPP

// Std lib
#include <cmath>
#include <limits>
#include <iostream>
#include <string>
#include <cstdio>
// Local
#include "globals.h"
#include "helper.hpp"

using namespace std;

class Args {

    public:

        static const string USAGE;
        static void readCmdArguments(int argc, char* argv[]);
        static void printProgramParameters();

};

#endif