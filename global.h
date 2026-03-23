#pragma once

#include <vector>
#include <complex>      // Used for complex numbers in path integral calculations

///// Vectors for data /////

extern std::vector<double> E0Therm;
extern std::vector<double> E0;
extern std::vector<double> accRateTherm;
extern std::vector<double> accRate;
extern std::vector<double> histogram;
extern std::vector<double> GTwo;
extern std::vector<double> GFour;
extern std::vector<double> instantons;
extern std::vector<double> antiInstantons;
extern std::vector<double> headerInfo;   // Stores the parameters of the run, written to the h5 file as an attribute for reference



///// Lattice parameters /////

extern const int N;
extern std::vector<double> positions;
extern const double a;

///// QHO specific parameters /////

extern const int m;
extern const int omega;

///// AHO specific parameters /////

extern const double quarticFactor;

///// DWP specific parameters /////

extern const double lambda;
extern const double wellCentres;

///// Mathematical constants /////

const double pif = 3.14159265358979323846;      // Pi as a double 
const std::complex <double> i(0.0, 1.0);	    // Imaginary unit