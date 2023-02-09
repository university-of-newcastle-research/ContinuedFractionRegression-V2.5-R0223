
/*
 * File:    debug.hpp
 * Authors: aciezak
 * 
 * Description
 * - Debug functions for cleaner and simplified code for output
 * 
 * Improvements
 * -
 *
 * 210000, AC, initial revision
 */

#ifndef DEBUG_HPP
#define DEBUG_HPP

// Std Lib
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

extern ofstream fileSummary;
extern ofstream fileModels;
extern ofstream fileSolution;

extern string   now();
extern void     dbg(ostream &out, string s, bool newline = true, bool time = true);
extern void     dbgc(ostream &out, string s);
extern void     dbgnl(ostream &out, string s = "");
extern void     in();
extern void     out();
extern string   ind();

#endif