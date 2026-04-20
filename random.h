#pragma once

#include <random>
#include "global.h"

/*************************************************/
/*********** Random number generation ************/
/*************************************************/

const int seed = 500;

typedef std::mt19937 MyRng; // Mersenne twister random number generator, good quality and fast.

inline MyRng& globalRng() {
    static MyRng rng;
    static bool seeded = false;
    if (!seeded) {
        rng.seed(seed);
        seeded = true;
    }
    return rng;
}

// Uniform [0,1)
inline double uniform01(MyRng& rng) {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(rng);
}

// Uniform [a,b)
inline double uniformRange(double a, double b, MyRng& rng) {
    return a + (b - a) * uniform01(rng);
}

// Uniform [-1,1)
inline double uniformMinus1to1(MyRng& rng) {
    return 2.0 * uniform01(rng) - 1.0;
}