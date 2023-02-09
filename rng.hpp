/*
 * File:   rng.hpp
 * Author: Haoyuan Sun
 *
 * Description
 * - Wrapper fo random number generator 
 * 
 * Improvement List
 * -
 * 
 * 210713, AC, Cleanup
 * 190405, HS, Initial version 
 */

using namespace std;

#ifndef __RNG_HPP
#define __RNG_HPP

#include <random>
#include <iostream>
#include <stdexcept>

class randint {

    private:
        mt19937 gen;

    public:

        randint() {
            random_device rd;
            gen.seed(rd());
        }

        randint(int seed) {
            gen.seed(seed);
        }

        int operator()(int a, int b) {
            // return uniform value inclusively of a and b
            if(a > b) {
                throw invalid_argument( "Invalid Range a:"+to_string(a)+" b:"+to_string(b) );
            }
             
            return uniform_int_distribution<> (a, b) (gen);
        }
};

class randreal {

    private:

        mt19937 gen;
        uniform_real_distribution<> dis;

    public:

        randreal() {
            random_device rd;
            gen.seed(rd());
        }

        randreal(int seed) {
            gen.seed(seed);
        }

        double operator()() {
            return uniform_real_distribution<> (0.0, 1.0) (gen);
        }

        double operator()(double a, double b) {
            // return uniform value inclusively of a and b
            if(a > b) 
                   throw invalid_argument( "Invalid Range a:"+to_string(a)+" b:"+to_string(b) );
                
            return uniform_real_distribution<> (a, b) (gen);
        }
};

#endif // __RNG_HPP
