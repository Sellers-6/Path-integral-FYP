#pragma once

#include <vector>
#include <complex>

/*****************************************************************/
/*********** External declaration of global variables ************/
/*****************************************************************/

///// Vectors for data /////

extern std::vector<double> E0Therm;
extern std::vector<double> E0;
extern std::vector<double> accRateTherm;
extern std::vector<double> accRate;
extern std::vector<double> histogram;
extern std::vector<double> Gx1x1;
extern std::vector<double> Gx2x2;
extern std::vector<double> instantons;
extern std::vector<double> antiInstantons;
extern std::vector<double> headerInfo;

///// Lattice parameters /////

extern int N;
extern std::vector<double> positions;
extern double a;

///// QHO specific parameters /////

extern const int m;
extern const int omega;

///// DWP specific parameters /////

extern const double lambda;
extern double wellCentres;

///// Mathematical constants /////

const double pif = 3.14159265358979323846;      // Pi
const std::complex <double> i(0.0, 1.0);	    // Imaginary unit