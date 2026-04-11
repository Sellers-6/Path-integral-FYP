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
const bool modifyGTwo = false;

///// Acceptance rate settings /////

const double epsilon = 0.1;				        // Maximum random displacement for Metropolis algorithm, decreasing epsilon increases acceptance rate. Want an acceptance rate between 50% and 80%. Lower rate is better for DWP so it can autocorrelate faster.
const int accRateInterval = 1000;               // Number of sweeps between recording the acceptance rate of the Metropolis algorithm

///// Decorrelation settings /////

int decorrSweeps;                               // Set by user input based on the system being simulated
const int decorrSweepsQHO = 50000;			        // Number of sweeps between taking measures of the path to reduce correlation between successive measures
const int decorrSweepsDWP = 2000;			        // Decorrelation takes longer in the DWP system
const int measures = 50;                       // Number of measures taken after thermalisation

///// Thermalisation settings /////

int thermSweeps;                                // Set by user input based on the system being simulated
const int thermSweepsQHO = 250000;       // Number of iterations for thermalisation, system is assumed to be thermalised after this many sweeps 
const int thermSweepsDWP = 20000;     // Thermalisation also takes longer in the DWP system
const int thermInterval = 10;    // Number of MC sweeps performed between measuring parameters during thermalisation

///// Initialisation settings /////

const bool hot_start = false;
const bool split_wells = false;
const double max_distance = 4;
int side = 1;    // For the split wells initialisation, determines which well the particle starts in (1 for right, -1 for left)

///// Repeats /////

const int repeats = 32;                          // Number of repeats for finding standard error 
bool multThreads = false;                      // Flag to determine whether to run the metropolis function in multiple threads, changed by user input

///// Lattice parameters /////

const int N = 50000;												// Number of lattice points. This discretises the imaginary time, so increasing N increases the accuracy of the simulation
std::vector<double> positions = std::vector<double>(N, 0.0);	// Lattice points (represents the "path" of the particle)
const double a = 0.01;											// Lattice spacing. Through the lattice spacing we define beta = N * a, the inverse temperature of the system. Making beta larger allows us to project out the ground state more effectively.
const double aInverse = 1.0 / a;											

///// QHO specific parameters /////

const int m = 1;                    // Unit mass
const int omega = 1;                // Unit frequency

///// DWP specific parameters /////

const double wellCentres = 1.5;     // Well centre positions, increasing this moves the wells further apart
const double lambda = 1;          // Coupling constant, increasing this deepens the wells and increases the barrier between them

const double omegaDWP = std::sqrt(8 * lambda * wellCentres * wellCentres);  // Frequency of the wells in the double well potential is equal to the square root of the second derivative of the potential at the minima, which is 8 * lambda * wellCentres^2.
// To use that the ground and first excited states are centred around 0.5, we require that omegaDWP = 1, which gives the relation lambda = 1 / (8 * wellCentres^2). 
const double tunnellingThreshold = 0.2 * wellCentres;     // Threshold for determining whether the particle is in the left or right well

///// Vectors to store data /////

std::vector<double> E0Therm;        //
std::vector<double> E0;       //
std::vector<double> accRateTherm;   //
std::vector<double> accRate;	//
std::vector<double> GTwo;     //
std::vector<double> GFour;    //
std::vector<double> histogram;
std::vector<double> instantons;
std::vector<double> antiInstantons;
std::vector<double> headerInfo;

///// Boundary conditions ///// 

int start;                          // Determines the left most point in the lattice which the Metropolis algorithm will loop through, useful for switching between periodic and Dirichlet boundary conditions (0 for periodic, 1 for Dirichlet)
int end;                            // Like the start variable, determines which point in the lattice the Metropolis algorithm will loop to (N for periodic, N-1 for Dirichlet)

///// Boudaries of the histogram /////

double xMax;                       // Maximum x value for the histogram of positions, set based on the maximum position reached during thermalisation
double xMin;                       // Minimum x value for the histogram of positions, set based on the maximum position reached during thermalisation
double binWidth;                   // Bin width for the histogram, set based on the histogram range and the number of bins

///// Shared data between threads /////

struct RepeatData {
    // Path positions
    std::vector<double> positions;          // All particle positions recorded for decorrelated measurements

    // Measurement vectors (and vacuum piece)
    std::vector<double> E0Temp;     // Ground-state energy measurements during this repeat
    std::vector<double> GTwoTemp;   // Two-point correlator accumulator for this repeat
    std::vector<double> GFourTemp;  // Four-point correlator accumulator for this repeat
    double vacuumPiece;             // For four-point correlator subtraction (per-thread)
    double xMean;                   // For two-point connected correlator subtraction (per-thread)
    std::vector<double> accRateTemp;    // Acceptance rate per repeat (decorrelation steps)
    std::vector<double> histogramTemp;
    std::vector<double> instantonsTemp;
    std::vector<double> antiInstantonsTemp;

    // Thermalisation data
    std::vector<double> E0ThermTemp;      // Temporary accumulator to check thermalisation
    std::vector<double> accRateThermTemp;     // Acceptance rate during thermalisation

    // Metropolis counters
    int sweep;           // Current sweep number
    int measureCount;    // Number of measurements taken so far
    int acceptedMoves;   // Count of accepted moves since last measurement

    // Constructor to initialise vectors
    RepeatData(int N = 0, int numBins = 0) {
        positions = std::vector<double>(N, 0.0);

        E0Temp.clear();
        GTwoTemp = std::vector<double>(N, 0.0);
        GFourTemp = std::vector<double>(N, 0.0);
        vacuumPiece = 0.0;
        xMean = 0.0;
        accRateTemp.clear();
        histogramTemp = std::vector<double>(numBins, 0.0);
        instantonsTemp.clear();
        antiInstantonsTemp.clear();

        E0ThermTemp.clear();
        accRateThermTemp.clear();

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

void thermalise(bool winOn, double (*potentialDifferential)(double), double (*potential)(double), std::mt19937& rng, RepeatData& data);

void takeThermMeasures(std::vector<double>& positions, double (*potentialDifferential)(double), double (*potential)(double), RepeatData& data);

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
    if (system == "QHO") { // Quantum harmonic oscillator
        return QHO::potential;
    }
    else if (system == "DWP") { // Double-well potential
        return DWP::potential;
    }
}

inline double(*findPotentialDifferential(const std::string& system))(double) {
    if (system == "QHO") { // Quantum harmonic oscillator
        return QHO::potentialDifferential;
    }
    else if (system == "DWP") { // Double-well potential
        return DWP::potentialDifferential;
    }
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

	int halfTime = N / 2;

    for (int t = 0; t < N; t++) {

        double position_t = positions[t];

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

void constructHeaderInfo(const std::string& system) {
    headerInfo.clear();
	headerInfo.push_back(N);    // Lattice points
	headerInfo.push_back(a);    // Lattice spacing
	headerInfo.push_back(epsilon);  // Maximum random displacement for Metropolis algorithm
	headerInfo.push_back(accRateInterval);  // Number of sweeps between recording the acceptance rate of the Metropolis algorithm
	headerInfo.push_back(decorrSweeps); // Number of sweeps between taking measures to reduce correlation
	headerInfo.push_back(thermSweeps);  // Number of sweeps for thermalisation, after which we assume the system is thermalised and take measures
	headerInfo.push_back(thermInterval);  // Number of sweeps between measuring parameters during thermalisation
	headerInfo.push_back(measures); // Number of measures taken after thermalisation
    headerInfo.push_back(repeats);   // Number of repeats for finding standard error
    headerInfo.push_back(numBins);  // Number of bins for the histogram of positions
    if (system == "QHO") {
        headerInfo.push_back(m);        // Mass for the quantum harmonic oscillator potential
        headerInfo.push_back(omega);    // Frequency for the quantum harmonic oscillator potential
    }
    else if (system == "DWP") {
        headerInfo.push_back(wellCentres);  // Well centre positions for the double well potential
        headerInfo.push_back(lambda);       // Coupling constant for the double well potential
	}
}