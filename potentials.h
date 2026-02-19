#pragma once

#include "global.h"

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
        return 0.5 * m * omega * omega * x * x;
    }
    static double potentialDifferential(double x) {
        return m * omega * omega * x;
    }
};


/////// Anharmonic Oscillator ///////


class AHO
{
public:
    static double potential(double x) {
        return 0.5 * m * omega * omega * x * x + quarticFactor * (x * x * x * x);
    }
    static double potentialDifferential(double x) {
        return m * omega * omega * x;
    }
};


/////// Double-Well Potential ///////


class DWP
{
public:
    
    static double potential(double x) {
        return (lambda / 24) * (x * x - wellCentres * wellCentres) * (x * x - wellCentres * wellCentres);
    }
    static double potentialDifferential(double x) {
        return (lambda / 6) * (((x * x * x) - x * (wellCentres * wellCentres)));
    }
};