#pragma once
#include <random>


//// Header file for random number generation functions ////


// TLDR: Use the following line of code to produce a random (ish) number between a and b:
// float randomNumber = rfRange(a, b);

typedef std::mt19937 MyRng; // certain generator for random numbers

inline MyRng& globalRng() {
    static MyRng rng;
    static bool seeded = false;
    if (!seeded) {
        rng.seed(6);
        seeded = true;
    }
    return rng;
}

static std::uniform_real_distribution<float> rf(0.0f, 1.0f);
double rfRange(double left, double right) {
    auto& rng = globalRng();
    double x = rf(rng);
    if (left < right) {
        return x * (right - left) + left; // Produces a random number between left and right boundaries with uniform dist
    }
    else if (left == right) {
        std::cout << "Warning! Left boundary was equal to right boundary! Using left as 0 and right as 1." << std::endl;
        return x;
    }
    else {
        std::cout << "Warning! Left boundary was larger than right boundary! Flipping the boundaries as easy fix." << std::endl;
        return x * (left - right) + right;
    }
}