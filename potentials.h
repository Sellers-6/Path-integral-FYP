#pragma once

#include "global.h"

double (*potential)(double);
double (*potentialDifferential)(double);


/////// Quantum Harmonic Oscillator ///////


class QHO
{
public:
    inline static double potential(double x) {
        return 0.5 * m * omega * omega * x * x;
    }
    inline static double potentialDifferential(double x) {
        return m * omega * omega * x;
    }
};


/////// Double-Well Potential ///////


class DWP
{
public:
    inline static double potential(double x) {
        return lambda * (x * x - wellCentres * wellCentres) * (x * x - wellCentres * wellCentres);
    }
    inline static double potentialDifferential(double x) {
        return lambda * 4 * (((x * x * x) - x * (wellCentres * wellCentres)));
    }
};