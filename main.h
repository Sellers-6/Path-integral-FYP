#pragma once

#include <iostream>     // Used for standard input and output streams
#include <vector>       // Used for dynamic arrays such as paths
#include <math.h>       // Used for mathematical functions such as pow()
#include <string>       // Generally useful for dealing with strings
#include <thread>       // Used to run the window in a separate thread

#include "random.h"
#include "potentials.h"
//#include "csv.h" // Legacy csv writing functions, replaced by h5 files
#include "h5.h"
#include "window.h"


///// Simulation settings /////

bool takeMeasuresFlag = true;    // Flag to determine whether to take measures after thermalisation (used to test thermalisation)

///// Acceptance rate settings /////

const double epsilon = 0.2;				        // Maximum random displacement for Metropolis algorithm, decreasing epsilon increases acceptance rate
const int accRateInterval = 1000;               // Number of sweeps between recording the acceptance rate of the Metropolis algorithm

///// Decorrelation settings /////

const int decorrelation = 250;			        // Number of sweeps between taking measures of the path to reduce correlation between successive measures
const int measures = 100;                       // Number of measures taken after thermalisation

///// Thermalisation settings /////

const double thermalisationConstant = 0.0004;   // Constant for checking thermalisation, decreasing the constant makes the check more strict, increasing it makes it less strict
const int thermalisationCheckLimit = 15;        // How many consecutive times the thermalisation check must be passed before we consider the system thermalised
const int thermalisationDecrement = 5;          // How much the thermalisation check counter is decremented by if the check fails
const int thermalisationMaximum = 100000;       // Maximum number of iterations for thermalisation, system is assumed to be thermalised after this many sweeps 
const int thermalisationMinimum = 1000;       // Minimum number of iterations for thermalisation
const int thermalisationMeasureInterval = 10;    // Number of MC sweeps performed between measuring parameters during thermalisation

///// Repeats /////

const int repeats = 5;                          // Number of repeats for finding standard error

///// Lattice parameters /////

const int N = 5000;												// Number of lattice points
std::vector<double> positions = std::vector<double>(N, 0.0);	// Lattice points (represents the "path" of the particle)
const double a = 0.1;											// Lattice spacing (I don't understand why this must be set to ~ 0.01 - 0.1)

///// QHO specific parameters /////

const int m = 1;                    // Unit mass
const int omega = 1;                // Unit frequency

///// DWP specific parameters /////

const double lambda = 2.0;          // Coupling constant, increasing this deepens the wells and increases the barrier between them
const double wellCentres = 2.0;     // Well centre positions, increasing this moves the wells further apart

///// Vectors to store data /////

std::vector<double> E0Therm;     // Vector to store the evolution of the ground state energy during thermalisation
std::vector<double> E0Decorr;        // Vector displaying the evolution of the ground state energy estimates over iterations
std::vector<double> accRateTherm;     // Vector to store information about the acceptance rate over measures
std::vector<double> accRateDecorr;	    // Vector to store the two-point correlation function 
std::vector<double> GTherm;  // Vector to store the positions per path, used to produce the wavefunction
std::vector<double> GDecorr;              // Vector to store ground state energies
std::vector<double> psiTherm;              // Vector to store first excited energies
std::vector<double> psiDecorr;              // Vector to store first excited energies
std::vector<double> thermSweeps;              

///// Boundary conditions ///// 

int start;                          // Determines the left most point in the lattice which the Metropolis algorithm will loop through, useful for switching between periodic and Dirichlet boundary conditions (0 for periodic, 1 for Dirichlet)
int end;                            // Like the start variable, determines which point in the lattice the Metropolis algorithm will loop to (N for periodic, N-1 for Dirichlet)

///// Counters /////

int sweep;                          // Counts how many metropolis sweeps have been performed, useful for checking when to take measures
int measureCount;				    // Counter for the number of measures taken within one test
int acceptedMoves;				    // Counts how many of the metropolis updates were accepted 

///// Function declarations /////

void chooseSystem();                                      // Allows user to choose which system to perform metropolis algorithm on

void metropolisRepeat(bool winOn, std::string boundary, std::string system);      // Repeats the metropolis function "repeats" times for histogram production

void metropolis(bool winOn, std::string boundary, std::string system, int test);            // Parent function for metropolisUpdate, runs loops to execute sufficient updates to the path

void metropolisUpdate(bool winOn, double (*potential)(double));  // Mechanics of the Metropolis algorithm such as checking whether a proposed update is accepted

void initialise(std::string boundary, std::string system);

const int thermalise(bool winOn, double (*potentialDifferential)(double), double (*potential)(double));

void takeThermMeasures(std::vector<double> positions, double (*potentialDifferential)(double), double (*potential)(double));

void takeMeasures(std::vector<double> positions, double (*potentialDifferential)(double), double (*potential)(double));

///// Helper functions /////

void setBoundary(std::string boundary);

double(*findPotential(const std::string& system))(double);

double(*findPotentialDifferential(const std::string& system))(double);

static double euclideanActionElement(double xin, double xfin, double (*potential)(double));  // Calculates the Euclidean action for two points on a path

static double E0Calc(const std::vector<double>& path, double (*potentialDifferential)(double), double (*potential)(double)); // Calculates the ground state energy estimate for a given path using the virial theorem

static std::vector<double> twoPointCorrelator(const std::vector<double>& path); // Calculates the two point correlator 

void setBoundary(std::string boundary) {
    if (boundary == "Periodic") {         // Periodic boundary conditions
        start = 0, end = N;
    }
    else if (boundary == "Dirichlet") {    // Dirichlet boundary conditions
        start = 1, end = N - 1;        // Effectively fixes the endpoints to 0, as they are never updated. This way we can apply Dirichlet boundary conditions to the double well potential. We effectively set one of the wells to be the "zero point" of the potential, so the endpoints are fixed at the minimum of one of the wells
    }
}

double(*findPotential(const std::string& system))(double) {
    if (system == "FP") {       // Free particle
        return FP::potential;
    }
    else if (system == "QHO") { // Quantum harmonic oscillator
        return QHO::potential;
    }
    else if (system == "DWP") { // Double-well potential
        return DWP::potential;
    }
}

double(*findPotentialDifferential(const std::string& system))(double) {
    if (system == "FP") {       // Free particle
        return FP::potentialDifferential;
    }
    else if (system == "QHO") { // Quantum harmonic oscillator
        return QHO::potentialDifferential;
    }
    else if (system == "DWP") { // Double-well potential
        return DWP::potentialDifferential;
    }
}

static double euclideanActionElement(double xin, double xfin, double (*potential)(double)) {
    return ((m / (2 * a)) * pow((xfin - xin), 2)) + (a * potential(xin));
}	    				

static double E0Calc(const std::vector<double>& path, double (*potentialDifferential)(double), double (*potential)(double)) {
    double E0Count = 0;
    for (int i = 0; i < path.size(); i++) {
        E0Count += (0.5 * path[i] * potentialDifferential(path[i]) + potential(path[i]));
	}
	return E0Count / N;
}

static std::vector<double> twoPointCorrelator(const std::vector<double>& path) {
    std::vector<double> correlationTemp(N, 0.0);

    for (int t = 0; t < N; t++) {
        for (int n = 0; n < N; n++) {
            correlationTemp[n] += (path[(t + n) % N] * path[t]) / N;
        }
    }
    return correlationTemp;
}