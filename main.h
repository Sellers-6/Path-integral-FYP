#pragma once

#include <iostream>
#include <vector>
#include <numeric>
#include <string>
#include <thread>                               // Used to run the visualisation window in a separate thread.
#include <omp.h>                                // Used for parallelisation, massively reduces the code execution time.

#include "random.h"
#include "potentials.h"
#include "h5.h"
#include "window.h"

/***********************************************************************/
/*********** Parameters, data vectors, and helper functions ************/
/***********************************************************************/

///// Acceptance rate settings /////

double epsilon = 0.4;				            // Maximum random displacement for Metropolis algorithm, decreasing epsilon increases acceptance rate. Want an acceptance rate between 50% and 80%.
const int accRateInterval = 1000;               // Number of sweeps between recording acceptance rate.

///// Decorrelation and Thermalisation /////

int thermSweeps;                                // Sweeps to thermalise the system.
const int thermSweepsQHO = 20000;               // About twice the decorrelation sweeps.
const int thermSweepsDWP = 2500;
const int thermInterval = 10;                   // Number of MC sweeps performed between measuring parameters during thermalisation.

int decorrSweeps;                               // Number of sweeps between taking measures of the path to reduce correlation between successive measures.
const int decorrSweepsQHO = 10000;			     
const int decorrSweepsDWP = 1250;			    // 1250 sweeps for option 7, 2500 sweeps for option 6.

///// Repeats /////

const int repeats = 32;                         // Number of repeats for finding standard error.
const int measures = 50;                        // Number of measures per repeat.
bool multThreads = false;                       // Flag to determine whether to run the metropolis function in multiple threads, changed by user input.

///// Lattice parameters /////

int N = 5000;									// Number of lattice points.
double a = 0.1;								// Lattice spacing. Beta = N * a, the inverse temperature of the system. Making beta larger allows us to project out the ground state more effectively.
double aInverse = 1.0 / a;										
std::vector<double> positions = std::vector<double>(N, 0.0);	// Represents the Euclidean "path" of the particle.


///// QHO specific parameters /////

const int m = 1;                                // Unit mass.
const int omega = 1;                            // Unit harmonic frequency.

///// DWP specific parameters /////

double wellCentres = 1.4;                       // Well centre positions.
const double lambda = 1.0;                      // Coupling constant. Kept as one, changing lambda just changes the scale of the system.

double omegaDWP = std::sqrt(8 * lambda * wellCentres * wellCentres);  // Frequency of the wells in the DWP.
double tunnellingThreshold = 0.2 * wellCentres; // Threshold for determining whether the particle is in the left or right well

///// Histogram for wave function /////

const int numBins = 100;                        // Number of bins for the histogram of positions.
double xMax;                                    // Maximum x value for the histogram of positions, set based on the maximum position reached during thermalisation.
double xMin;                                    // Minimum x value for the histogram of positions, set based on the minimum position reached during thermalisation.
double binWidth;                                // Bin width for the histogram. Set based on xMax, xMin and numBins.

///// Initialisation settings /////

const bool hotStart = false;                   // Start the path with random configuration.
const double maxDistance = 4;                  // Maximum distance from the origin for the initial positions of the path, used for the hot start initialisation.

int side = 1;                                   // Determines which well the particle starts in (1 for right, -1 for left)

///// Flags /////

bool sixFlag = false;                           // Flag to determine whether user selected option six or not.
bool sevenFlag = false;                         // Flag to determine whether user selected option seven or not.
bool eightFlag = false;                         // Flag to determine whether user selected option eight or not.

///// Vectors to store data /////

std::vector<double> accRateTherm;               // Acceptance rate measurements during thermalisation.
std::vector<double> accRate;	                // Acceptance rate measurements after thermalisation.
std::vector<double> E0Therm;                    // Ground state energy measurements during thermalisaiton.
std::vector<double> E0;                         // Decorrelated ground state energy measurements.
std::vector<double> histogram;                  // Histogram of positions for finding the wave function.
std::vector<double> Gx1x1;                      // Two-point connected correlator G(n) = G(x(t), x(t+n)).
std::vector<double> Gx1x2;                      // Two-point connected correlator G(n) = G(x(t), x^2(t+n)) = G(x^2(t), x(t+n)).
std::vector<double> Gx2x2;                      // Two-point connected correlator G(n) = G(x^2(t), x^2(t+n)).
std::vector<double> instantons;                 // Number of instantons.
std::vector<double> antiInstantons;             // Number of anti-instantons.
std::vector<double> headerInfo;                 // Vector to store the parameters of a simulation, useful for data analysis.

///// Shared data between threads /////

struct RepeatData {
    // Path positions
    std::vector<double> positions;              // Particle positions on the lattice.
    std::vector<std::vector<double>> positionsVector; // Vector of all positions measured, stored for computing observables at the end of the repeat.

    // Temporary measurement vectors 
    std::vector<double> accRateThermTemp;
    std::vector<double> accRateTemp;
    std::vector<double> E0ThermTemp;
    std::vector<double> E0Temp;
    std::vector<double> histogramTemp;
	std::vector<double> Gx1x1Temp;
    std::vector<double> Gx1x2Temp;  
    std::vector<double> Gx2x2Temp;
    std::vector<double> instantonsTemp;
    std::vector<double> antiInstantonsTemp;

    // MCMC counters
    int sweep;           // MCMC sweep counter.
    int measureCount;    // Number of measurements taken.
    int acceptedMoves;   // Accepted moves since last measurement of acceptance rate.

    // Constructor to initialise vectors
    RepeatData(int N = 0, int numBins = 0, int measures = 0) {
        positions = std::vector<double>(N, 0.0);
        positionsVector = std::vector<std::vector<double>>(measures, std::vector<double>(N, 0.0));

        accRateThermTemp.clear();
        accRateTemp.clear();
        E0ThermTemp.clear();
        E0Temp.clear();
        histogramTemp = std::vector<double>(numBins, 0.0);
        Gx1x1Temp = std::vector<double>(N, 0.0);
        Gx1x2Temp = std::vector<double>(N, 0.0);
        Gx2x2Temp = std::vector<double>(N, 0.0);
        instantonsTemp.clear();
        antiInstantonsTemp.clear();

        sweep = 0;
        measureCount = 0;
        acceptedMoves = 0;
    }
};

///// Function declarations /////

void chooseSystem();

void metropolisRepeat(std::string system);

void metropolis(std::string system, int repeat, std::mt19937& rng, RepeatData& data, double (*potential)(double), double (*potentialDifferential)(double));

void metropolisSweep(double (*potential)(double), std::mt19937& rng, RepeatData& data);

void initialise(std::string system, std::mt19937& rng, RepeatData& data);

void thermalise(double (*potentialDifferential)(double), double (*potential)(double), std::mt19937& rng, RepeatData& data);

void takeThermMeasures(std::vector<double>& positions, double (*potentialDifferential)(double), double (*potential)(double), RepeatData& data);

void takeMeasures(std::vector<double>& positions, RepeatData& data);

void computeObservables(std::vector<std::vector<double>>& positionsVector, double (*potentialDifferential)(double), double (*potential)(double), RepeatData& data);

///// Helper functions /////

double(*findPotential(const std::string& system))(double) {
    if (system == "QHO") {                      // Quantum harmonic oscillator.
        return QHO::potential;
    }
    else if (system == "DWP") {                 // Double-well potential.
        return DWP::potential;
    }
    else {
		return QHO::potential;                  // Default to QHO potential.
    }
}

double(*findPotentialDifferential(const std::string& system))(double) {
    if (system == "QHO") {
        return QHO::potentialDifferential;
    }
    else if (system == "DWP") {
        return DWP::potentialDifferential;
    }
    else {
        return QHO::potentialDifferential;
    }
}

static double E0Calc(const std::vector<double>& positions, double (*potentialDifferential)(double), double (*potential)(double)) {
    double E0Count = 0;
    for (int i = 0; i < positions.size(); i++) {
        E0Count += (0.5 * positions[i] * potentialDifferential(positions[i]) + potential(positions[i]));
	}
	return E0Count / N;
}

static std::vector<double> correlator(const std::vector<double>& positionsLeft, const std::vector<double>& positionsRight) {

    std::vector<double> correlationTemp(N, 0.0);

	int halfTime = N / 2;

    for (int t = 0; t < N; t++) {
        double position_t = positionsRight[t];

        // Wrap around the correlator instead of using expensive modulo arithmetic
        int maxNoWrap = std::min(halfTime, N - t - 1);
        for (int n = 0; n <= maxNoWrap; n++) {
            correlationTemp[n] += positionsLeft[(t + n)] * position_t;
        }

        for (int n = maxNoWrap + 1; n <= halfTime; n++) {
            correlationTemp[n] += positionsLeft[t + n - N] * position_t;
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

static int whichWell(double x, double threshold) {
    if (x > threshold) { return 1; }            // Right well.
	if (x < -threshold) { return -1; }          // Left well.
    return 0;                                   // Barrier.
}

void constructHeaderInfo(const std::string& system) {
    headerInfo.clear();
	headerInfo.push_back(N);                    // Lattice size.
	headerInfo.push_back(a);                    // Lattice spacing.
	headerInfo.push_back(epsilon);              // Maximum random displacement for MCMC.
	headerInfo.push_back(accRateInterval);      // Number of sweeps between recording the acceptance rate.
	headerInfo.push_back(decorrSweeps);         // Number of sweeps between taking measures.
	headerInfo.push_back(thermSweeps);          // Number of sweeps for thermalisation.
	headerInfo.push_back(thermInterval);        // Number of sweeps between measuring parameters during thermalisation.
	headerInfo.push_back(measures);             // Number of measures taken.
    headerInfo.push_back(repeats);              // Number of repeats/independent threads.
    headerInfo.push_back(numBins);              // Number of bins for the histogram of positions.
    if (system == "QHO") {
        headerInfo.push_back(m);                // Mass for QHO.
        headerInfo.push_back(omega);            // Frequency for QHO.
    }
    else if (system == "DWP") {
        headerInfo.push_back(wellCentres);      // Well centre positions for DWP.
        headerInfo.push_back(lambda);           // Coupling constant for the DWP.
	}
}

void setParameters(int beta, double a, double wellCentres) {
    // Use thermalisation/decorrelation sweeps for the system with smallest a, as this will be the slowest to thermalise and decorrelatese.
	// This is inefficient for systems that are fast to thermalise/decorrelate, but ensures that all systems are sufficiently thermalised and decorrelated.
    
    // Set lattice paramters
    N = beta / a;
	aInverse = 1.0 / a;
    
    // Set parameters for the DWP
    omegaDWP = std::sqrt(8 * lambda * wellCentres * wellCentres);
    tunnellingThreshold = 0.2 * wellCentres;
}