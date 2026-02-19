#pragma once

#include <iostream>     // Used for standard input and output streams
#include <vector>       // Used for dynamic arrays such as paths
#include <string>       // Generally useful for dealing with strings
#include <thread>       // Used to run the window in a separate thread
#include <numeric>  	// Used for std::accumulate 
//#include <omp.h>        // Used for parallelisation of the metropolis function, massively reduces the code execution time

#include "random.h"
#include "potentials.h"
//#include "csv.h" // Legacy csv writing functions, replaced by h5 files
#include "h5.h"
#include "window.h"


///// Simulation settings /////

bool takeMeasuresFlag = true;    // Flag to determine whether to take measures after thermalisation (used to test thermalisation)

///// Acceptance rate settings /////

const double epsilon = 0.3;				        // Maximum random displacement for Metropolis algorithm, decreasing epsilon increases acceptance rate. Found to be best to set epsilon ~ 0.3 for both DWP and QHO, both of which result in an acceptance rate of ~ 74%
const int accRateInterval = 1000;               // Number of sweeps between recording the acceptance rate of the Metropolis algorithm

///// Decorrelation settings /////

const int decorrelation = 250;			        // Number of sweeps between taking measures of the path to reduce correlation between successive measures
const int measures = 100;                       // Number of measures taken after thermalisation

///// Thermalisation settings /////

const double acceptableError = 0.01;              // Ratio of the standard error to the mean for the ground state energy, used as a criterion for thermalisation
const int thermalisationMaximum = 200000;       // Maximum number of iterations for thermalisation, system is assumed to be thermalised after this many sweeps 
const int thermalisationMinimum = 1000;       // Minimum number of iterations for thermalisation
const int thermalisationInterval = 100;    // Number of MC sweeps performed between measuring parameters during thermalisation
// Be careful when changing the thermalisationInterval to be too small; this can massively increase file size
std::vector<double> E0ThermTemp;            // Used for creating batches in one iteration of the thermalisation process

///// Repeats /////

//const int threads = 6;                          // Number of threads to run in parallel, set to the number of cores on my computer (not yet implemented)
const int repeats = 5;                          // Number of repeats for finding standard error (threads * repeats measures are taken in total)
//const bool multThreads = false;                      // Flag to determine whether to run the metropolis function in multiple threads (not yet working)

///// Lattice parameters /////

const int N = 5000;												// Number of lattice points
std::vector<double> positions = std::vector<double>(N, 0.0);	// Lattice points (represents the "path" of the particle)
const double a = 0.1;											// Lattice spacing (I don't understand why this must be set to ~ 0.01 - 0.1)

///// QHO specific parameters /////

const int m = 1;                    // Unit mass
const int omega = 1;                // Unit frequency

///// AHO specific parameters /////

const double quarticFactor = 1;     // Quartic factor for the anharmonic oscillator

///// DWP specific parameters /////

const double lambda = 2.0;          // Coupling constant, increasing this deepens the wells and increases the barrier between them
const double wellCentres = 2.0;     // Well centre positions, increasing this moves the wells further apart

///// Vectors to store data /////

std::vector<double> E0Therm;        //
std::vector<double> E0Decorr;       //
std::vector<double> accRateTherm;   //
std::vector<double> accRateDecorr;	//
std::vector<double> xTherm;         //
std::vector<double> xDecorr;        //
std::vector<double> GTwoDecorr;     //
std::vector<double> GFourDecorr;    //   
std::vector<double> thermSweeps;    //           

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

void takeThermMeasures(std::vector<double>& positions, double (*potentialDifferential)(double), double (*potential)(double));

bool checkThermalised();

void takeMeasures(std::vector<double>& positions, double (*potentialDifferential)(double), double (*potential)(double));

///// Helper functions /////

static void setBoundary(std::string boundary) {
    if (boundary == "Periodic") {         // Periodic boundary conditions
        start = 0, end = N;
    }
    else if (boundary == "Dirichlet") {    // Dirichlet boundary conditions
        start = 1, end = N - 1;        // Effectively fixes the endpoints to 0, as they are never updated. This way we can apply Dirichlet boundary conditions to the double well potential. We effectively set one of the wells to be the "zero point" of the potential, so the endpoints are fixed at the minimum of one of the wells
    }
}

static double(*findPotential(const std::string& system))(double) {
    if (system == "FP") {       // Free particle
        return FP::potential;
    }
    else if (system == "QHO") { // Quantum harmonic oscillator
        return QHO::potential;
    }
    else if (system == "AHO") { // Anharmonic oscillator
        return AHO::potential;
    }
    else if (system == "DWP") { // Double-well potential
        return DWP::potential;
    }
}

static double(*findPotentialDifferential(const std::string& system))(double) {
    if (system == "FP") {       // Free particle
        return FP::potentialDifferential;
    }
    else if (system == "QHO") { // Quantum harmonic oscillator
        return QHO::potentialDifferential;
    }
    else if (system == "AHO") { // Anharmonic oscillator
        return AHO::potentialDifferential;
    }
    else if (system == "DWP") { // Double-well potential
        return DWP::potentialDifferential;
    }
}

static double vectorMean(const std::vector<double>& vector) {
    return std::accumulate(vector.begin(), vector.end(), 0.0) / vector.size();
}

static double vectorVariance(const std::vector<double>& vector, double vectorMean) {
    double variance = 0.0;
    for (auto x : vector) variance += (x - vectorMean) * (x - vectorMean);
	return variance / (vector.size() - 1);  // Bessel's correction for sample variance
}

static double MCSE(const std::vector<double>& samples) {
    size_t N = samples.size();
    if (N < 2) return NAN;

    size_t batchSize = std::max<size_t>(1, static_cast<size_t>(std::sqrt(N)));
    size_t M = N / batchSize;
    if (M < 2) return NAN; // Need at least 2 batches

    std::vector<double> batchMeans(M, 0.0);
    for (size_t k = 0; k < M; ++k) {
        size_t start = k * batchSize;
        size_t end = start + batchSize;
        double sum = std::accumulate(samples.begin() + start, samples.begin() + end, 0.0);
        batchMeans[k] = sum / batchSize;
    }

    double bmMean = vectorMean(batchMeans);
    double bmVar = vectorVariance(batchMeans, bmMean);
    return std::sqrt(bmVar / M);
}

static double euclideanActionElement(double xin, double xfin, double (*potential)(double)) {
    return ((m / (2 * a)) * (xfin - xin) * (xfin - xin)) + (a * potential(xin));
}	    				

static double E0Calc(const std::vector<double>& positions, double (*potentialDifferential)(double), double (*potential)(double)) {
    double E0Count = 0;
    for (int i = 0; i < positions.size(); i++) {
        E0Count += (0.5 * positions[i] * potentialDifferential(positions[i]) + potential(positions[i]));
	}
	return E0Count / N;
}

static std::vector<double> twoPointCorrelator(const std::vector<double>& positions) {
    std::vector<double> correlationTemp(N, 0.0);

    for (int t = 0; t < N; t++) {
        for (int n = 0; n < N; n++) {
            correlationTemp[n] += (positions[(t + n) % N] * positions[t]) / N;
        }
    }
    return correlationTemp;
}

static double x2Mean(const double* x, size_t N)
{
    double sum = 0.0;

    for (size_t i = 0; i < N; ++i)
        sum += x[i] * x[i];

    return sum / N;
}

static std::vector<double> fourPointCorrelator(const std::vector<double>& positions) {
    std::vector<double> correlationTemp(N, 0.0);
    
    double x2 = x2Mean(positions.data(), positions.size());

    for (int t = 0; t < N; t++) {
        for (int n = 0; n < N; n++) {

            correlationTemp[n] += (positions[(t + n) % N] * positions[(t + n) % N]) * (positions[t] * positions[t]) / N;
        }
    }

    for (int n = 0; n < N; n++) {   // Remove the vacuum piece 

        correlationTemp[n] -= x2 * x2;
    }

    return correlationTemp;
}