#pragma once
#include "main.h"

double (*potential)(double);
double (*potentialDifferential)(double);

/////// Free Particle ///////


class FP
{
public:
    static double potential(double x) {
        return 0;
    }
    static double potentialDifferential(double x) {
        return 0;
    }
};


/////// Quantum Harmonic Oscillator ///////


class QHO
{
public:
    static double potential(double x) {
        return 0.5 * m * pow(omega, 2) * pow(x, 2);
    }
    static double potentialDifferential(double x) {
        return m * pow(omega, 2) * x;
    }
};


/////// Double-Well Potential ///////


class DWP
{
public:

    static double potential(double x) {
        return (lambda / 24) * pow(pow(x, 2) - pow(wellCentres, 2), 2);
    }
    static double potentialDifferential(double x) {
        return (lambda / 6) * ((pow(x, 3) - x * pow(wellCentres, 2)));
    }
};