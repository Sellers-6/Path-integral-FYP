#pragma once

#include <iostream>     // Used for standard input and output streams
#include <vector>       // Used for dynamic arrays such as paths
#include <string>       // Generally useful for dealing with strings
#include <thread>       // Used to run the window in a separate thread
#include <numeric>  	// Used for std::accumulate 
#include <omp.h>        // Used for parallelisation of the metropolis function, massively reduces the code execution time

#include "random.h"
#include "potentials.h"
//#include "csv.h" // Legacy csv writing functions, replaced by h5 files
#include "h5.h"
#include "window.h"


///// Simulation settings /////

const bool takeMeasuresFlag = true;    // Flag to determine whether to take measures (If this is turned off, thermalisation will only complete after thermalisationMaximum!)
const bool takeThermMeasuresFlag = true;
const int numBins = 100;              // Number of bins for the histogram of positions

///// Acceptance rate settings /////

const double epsilon = 0.4;				        // Maximum random displacement for Metropolis algorithm, decreasing epsilon increases acceptance rate. Want an acceptance rate between 50% and 80%. Lower rate is better for DWP so it can autocorrelate faster.
const int accRateInterval = 1000;               // Number of sweeps between recording the acceptance rate of the Metropolis algorithm

///// Decorrelation settings /////

const int decorrelation = 2500;			        // Number of sweeps between taking measures of the path to reduce correlation between successive measures. Decorrelation takes far longer in the DWP system!
const int measures = 50;                       // Number of measures taken after thermalisation

///// Initialisation settings /////

const bool hot_start = false;
const double max_distance = 4;

///// Thermalisation settings /////

// thermalisationMaximum is the default value for thermalisation sweeps if takeThermMeasuresFlag is set to false
const int thermalisationMaximum = 100000;       // Maximum number of iterations for thermalisation, system is assumed to be thermalised after this many sweeps 
const int thermalisationMinimum = 10000;       // Minimum number of iterations for thermalisation
const int thermalisationInterval = 10;    // Number of MC sweeps performed between measuring parameters during thermalisation
// Be careful when changing the thermalisationInterval to be too small; this can massively increase file size
const double acceptableError = 0.01;              // Ratio of the standard error to the mean for the ground state energy, used as a criterion for thermalisation
std::vector<double> E0ThermTemp;            // Used for creating batches in one iteration of the thermalisation process

///// Repeats /////

int repeats = 32;                          // Number of repeats for finding standard error 
bool multThreads = false;                      // Flag to determine whether to run the metropolis function in multiple threads 

///// Lattice parameters /////

const int N = 10000;												// Number of lattice points. This discretises the imaginary time, so increasing N increases the accuracy of the simulation
std::vector<double> positions = std::vector<double>(N, 0.0);	// Lattice points (represents the "path" of the particle)
const double a = 0.05;											// Lattice spacing. Through the lattice spacing we define beta = N * a, the inverse temperature of the system. Making beta larger allows us to project out the ground state more effectively.
const double aInverse = 1.0 / a;											

///// QHO specific parameters /////

const int m = 1;                    // Unit mass
const int omega = 1;                // Unit frequency

///// AHO specific parameters /////

const double quarticFactor = 1;     // Quartic factor for the anharmonic oscillator

///// DWP specific parameters /////

//const double wellCentres = 2.0;     // Well centre positions, increasing this moves the wells further apart
//const double lambda = 3 / (wellCentres * wellCentres);          // Coupling constant, increasing this deepens the wells and increases the barrier between them

const double wellCentres = 1.8;     // Well centre positions, increasing this moves the wells further apart
const double lambda = 12;          // Coupling constant, increasing this deepens the wells and increases the barrier between them

const double omegaDWP = std::sqrt(8 * (lambda / 24) * wellCentres * wellCentres);  // Frequency of the wells in the double well potential is equal to the square root of the second derivative of the potential at the minima, which is 8 * lambda * wellCentres^2.
// To use that the ground and first excited states are centred around 0.5, we require that omegaDWP = 1, which gives the relation lambda = 1 / (8 * wellCentres^2). 
const int tunnellingThreshold = 0.2 * wellCentres;     // Threshold for determining whether the particle is in the left or right well

///// Vectors to store data /////

std::vector<double> E0Therm;        //
std::vector<double> E0;       //
std::vector<double> accRateTherm;   //
std::vector<double> accRate;	//
std::vector<double> positionsTemp;        //
std::vector<double> GTwo;     //
std::vector<double> GFour;    //   
std::vector<double> thermSweeps;    //   
std::vector<double> histogram;
std::vector<double> instantons;
std::vector<double> antiInstantons;

///// Boundary conditions ///// 

int start;                          // Determines the left most point in the lattice which the Metropolis algorithm will loop through, useful for switching between periodic and Dirichlet boundary conditions (0 for periodic, 1 for Dirichlet)
int end;                            // Like the start variable, determines which point in the lattice the Metropolis algorithm will loop to (N for periodic, N-1 for Dirichlet)

///// Boudaries of the histogram /////

double xMax;                       // Maximum x value for the histogram of positions, set based on the maximum position reached during thermalisation
double xMin;                       // Minimum x value for the histogram of positions, set based on the maximum position reached during thermalisation
double binWidth;                   // Bin width for the histogram, set based on the histogram range and the number of bins

///// Counters /////

int sweep;                          // Counts how many metropolis sweeps have been performed, useful for checking when to take measures
int measureCount;				    // Counter for the number of measures taken within one test
int acceptedMoves;				    // Counts how many of the metropolis updates were accepted 
double vacuumPiece;                 // Counts the vacuum piece of the four point correlator
std::vector<double> GTwoTemp;
std::vector<double> GFourTemp;
std::vector<double> histogramTemp;
std::vector<double> instantonsTemp;	
std::vector<double> antiInstantonsTemp;

///// Shared data between threads /////

struct RepeatData {
    // Path positions
    std::vector<double> positions;        // Current path of the particle
    std::vector<double> positionsTemp;          // All particle positions recorded for decorrelated measurements

    // Measurement vectors (and vacuum piece)
    std::vector<double> E0Temp;     // Ground-state energy measurements during this repeat
    std::vector<double> GTwoTemp;   // Two-point correlator accumulator for this repeat
    std::vector<double> GFourTemp;  // Four-point correlator accumulator for this repeat
    double vacuumPiece;             // For four-point correlator subtraction (per-thread)
    std::vector<double> accRate;    // Acceptance rate per repeat (decorrelation steps)
    std::vector<double> histogramTemp;
    std::vector<double> instantonsTemp;
    std::vector<double> antiInstantonsTemp;

    // Thermalisation data
    std::vector<double> E0Therm;          // Ground-state energy during thermalisation
    std::vector<double> E0ThermTemp;      // Temporary accumulator to check thermalisation
    std::vector<double> accRateTherm;     // Acceptance rate during thermalisation

    // Metropolis counters
    int sweep;           // Current sweep number
    int measureCount;    // Number of measurements taken so far
    int acceptedMoves;   // Count of accepted moves since last measurement

	// Tunnelling data



    // Constructor to initialise vectors
    RepeatData(int N = 0, int numBins = 0) {
        positions = std::vector<double>(N, 0.0);
        positionsTemp.clear();

        E0Temp.clear();
        GTwoTemp = std::vector<double>(N, 0.0);
        GFourTemp = std::vector<double>(N, 0.0);
        vacuumPiece = 0.0;
        accRate.clear();
        histogramTemp = std::vector<double>(numBins, 0.0);
        instantonsTemp.clear();
        antiInstantonsTemp.clear();

        E0Therm.clear();
        E0ThermTemp.clear();
        accRateTherm.clear();

        sweep = 0;
        measureCount = 0;
        acceptedMoves = 0;
    }
};

///// Function declarations /////

void chooseSystem();                                      // Allows user to choose which system to perform metropolis algorithm on

void metropolisRepeat(bool winOn, std::string boundary, std::string system);      // Repeats the metropolis function "repeats" times

void metropolis(bool winOn, std::string boundary, std::string system, int repeat, std::mt19937& rng, RepeatData& data,
    double (*potential)(double), double (*potentialDifferential)(double));            // Parent function for metropolisUpdate, runs loops to execute sufficient updates to the path

void metropolisUpdate(bool winOn, double (*potential)(double), std::mt19937& rng, RepeatData& data);  // Mechanics of the Metropolis algorithm such as checking whether a proposed update is accepted

void initialise(std::string boundary, std::string system, std::mt19937& rng, RepeatData& data);

const int thermalise(bool winOn, double (*potentialDifferential)(double), double (*potential)(double), std::mt19937& rng, RepeatData& data);

void takeThermMeasures(std::vector<double>& positions, double (*potentialDifferential)(double), double (*potential)(double), RepeatData& data);

bool checkThermalised(const RepeatData& data);

void takeMeasures(std::vector<double>& positions, double (*potentialDifferential)(double), double (*potential)(double), RepeatData& data);

///// Helper functions /////

static void setBoundary(std::string boundary) {
    if (boundary == "Periodic") {         // Periodic boundary conditions
        start = 0, end = N;
    }
    else if (boundary == "Dirichlet") {    // Dirichlet boundary conditions
        start = 1, end = N - 1;        // Effectively fixes the endpoints to 0, as they are never updated. This way we can apply Dirichlet boundary conditions to the double well potential. We effectively set one of the wells to be the "zero point" of the potential, so the endpoints are fixed at the minimum of one of the wells
    }
}

inline double(*findPotential(const std::string& system))(double) {
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

inline double(*findPotentialDifferential(const std::string& system))(double) {
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

static double E0Calc(const std::vector<double>& positions, double (*potentialDifferential)(double), double (*potential)(double)) {
    double E0Count = 0;
    for (int i = 0; i < positions.size(); i++) {
        E0Count += (0.5 * positions[i] * potentialDifferential(positions[i]) + potential(positions[i]));
	}
	return E0Count / N;
}

static std::vector<double> twoPointCorrelator(const std::vector<double>& positions) {

    std::vector<double> correlationTemp(N, 0.0);

    int halfTime = N / 2;   // Use symmetry: only compute up to N/2

    for (int t = 0; t < N; t++) {

        double position_t = positions[t];

        // Since the modulo function is expensive, we wrap around by taking away N
        // To avoid branching in the inner loop, we split it into two parts: one for the non-wrapping case and one for the wrapping case.

        // Case 1: no wrap
        int maxNoWrap = std::min(halfTime, N - t - 1);
        for (int n = 0; n <= maxNoWrap; n++) {
            correlationTemp[n] += positions[t + n] * position_t;
        }

        // Case 2: wrap
        for (int n = maxNoWrap + 1; n <= halfTime; n++) {
            correlationTemp[n] += positions[t + n - N] * position_t;
        }
    }

    // Mirror symmetric half
    for (int n = 1; n < halfTime; n++) {
        correlationTemp[N - n] = correlationTemp[n];
    }

    // Normalise
    for (int n = 0; n < N; n++) {
        correlationTemp[n] /= N;
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

static std::vector<double> fourPointCorrelator(const std::vector<double>& positions, double& vacuumPiece) {
    std::vector<double> correlationTemp(N, 0.0);

    vacuumPiece += x2Mean(positions.data(), positions.size());

    int halfTime = N / 2;

	// Precompute x(t)^2 for all t to avoid redundant calculations in the inner loop
    std::vector<double> x2(N);
    for (int i = 0; i < N; i++)
        x2[i] = positions[i] * positions[i];

    for (int t = 0; t < N; t++) {

        double positionSquared = x2[t];

        // Case 1: no wrap
        int maxNoWrap = std::min(halfTime, N - t - 1);
        for (int n = 0; n <= maxNoWrap; n++) {
            correlationTemp[n] += x2[t + n] * positionSquared;
        }

        // Case 2: wrap
        for (int n = maxNoWrap + 1; n <= halfTime; n++) {
            correlationTemp[n] += x2[t + n - N] * positionSquared;
        }
    }

    // Mirror the symmetric half
    for (int n = 1; n < halfTime; n++) {
        correlationTemp[N - n] = correlationTemp[n];
    }

	// Normalise
    for (int n = 0; n < N; n++) {
        correlationTemp[n] /= N;
	}

    return correlationTemp;
}

static int whichWell(double x, double threshold)
{
    if (x > threshold)  return 1;   // right well
	if (x < -threshold) return -1;  // left well
    return 0;                       // barrier
}