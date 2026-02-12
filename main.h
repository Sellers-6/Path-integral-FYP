#pragma once

#define SDL_MAIN_HANDLED 

#include <iostream>     // Used for standard input and output streams
#include <fstream>      // Used for file input and output
#include <map>          // Used for mapping file keys to file streams for output
#include <stdio.h>      // Used for standard input and output
#include <vector>       // Used for dynamic arrays such as paths
#include <math.h>       // Used for mathematical functions such as pow()
#include <string>       // Used for file name strings
#include <complex>      // Used for complex numbers in path integral calculations
#include <SDL2/SDL.h>   // Used for creating a window to visualise the path updates
#include <thread>       // Used to run the window in a separate thread
#include "lattice.h"    // Custom lattice header file
#include "random.h"     // Custom random number generator header file


///// Variables that can be adjusted to influence the accuracy of results and runtime of the program /////

const int hbar = 1;		                        // Unit (reduced) planck's constant

const int m = 1, omega = 1;                     // Unit mass and frequency (harmonic oscillator); 

const double lambda = 2.0, wellCentres = 2.0;   // Coupling constant and well centre positions (double well potential), adjust to influence the shape of the potential 
                                                // Increasing lambda deepens the wells, increasing wellCentres moves the wells further apart

const double epsilon = 0.2;				        // Maximum random displacement for Metropolis algorithm, decreasing epsilon increases acceptance rate

const int decorrelation = 250;			        // Number of sweeps between taking measures of the path to reduce correlation between successive measures

const int accRateInterval = 1000;               // Number of sweeps between recording the acceptance rate of the Metropolis algorithm

const double thermalisationConstant = 0.0001;   // Constant for checking thermalisation, decreasing the constant makes the check more strict, increasing it makes it less strict
const int thermalisationCheckLimit = 10;        // How many consecutive times the thermalisation check must be passed before we consider the system thermalised
const int thermalisationMaximum = 100000;       // Maximum number of iterations for thermalisation, system is assumed to be thermalised after this many sweeps 
const int thermalisationMinimum = 1000;         // Minimum number of iterations for thermalisation

const int measures = 10;                       // Number of measures taken after thermalisation, adjust to influence accuracy of results 

const int tests = 1000;                            // Number of repeats for finding standard error

//// File naming and creation ////
std::ofstream logfile;

const std::string baseDir = "csv/";

enum class Boundary { Periodic, Dirichlet };
enum class System { QHO, DWP };

std::string boundaryName(Boundary b) {
    return (b == Boundary::Periodic) ? "Periodic" : "Dirichlet";
}

std::string systemName(System s) {
    return (s == System::QHO) ? "QHO" : "DWP";
}

const std::vector<std::string> quantities = {"E0Thermalisation", "E0Evolution", "acceptanceRate", "correlation", "waveFunction", "E0", "E1", "accRate"};

using FileKey = std::tuple<std::string, Boundary, System>;
std::map<FileKey, std::ofstream> files;

void openFiles(std::string boundary, std::string system) {
    for (const auto& q : quantities) {

        std::string filename =
            baseDir + q +
            boundary +
            system +
            ".csv";

        files.emplace(
            std::make_tuple(q, Boundary(boundary == "Periodic" ? Boundary::Periodic : Boundary::Dirichlet),
                System(system == "QHO" ? System::QHO : System::DWP)),
            std::ofstream(filename)
        );

    }
}


///// Mathematical constants /////

const double pif = 3.14159265358979323846;      // Pi as a double 
const std::complex <double> i(0.0, 1.0);	    // Imaginary unit


///// Vectors and arrays containing file input /////

// temp stuff for debugging //
int mestes = measures * tests;
int ntes = N * tests;

///

std::vector<double> E0Thermalising; // Vector to store the evolution of the ground state energy during thermalisation
std::vector<double> E0Evolution;    // Vector displaying the evolution of the ground state energy estimates over iterations
std::vector<double> acceptanceRate; // Vector to store information about the acceptance rate over measures
std::vector<double> G(ntes, 0.0);	    // Vector to store the two-point correlation function 
std::vector<std::vector<double>> psi(mestes, std::vector<double>(N, 0.0));    // Array containing the all the positions per path, used to reproduce the wavefunction
std::vector<double> E0Vec;          // Vector to store ground state energies
std::vector<double> E1Vec;          // Vector to store first excited energies
std::vector<double> accRateVec;     // Vector to store acceptance rates


///// Variable and counter declarations (should not need to be changed) /////

// Boundary conditions // 

int start;                          // Determines the left most point in the lattice which the Metropolis algorithm will loop through, useful for switching between periodic and Dirichlet boundary conditions (0 for periodic, 1 for Dirichlet)
int end;                            // Like the start variable, determines which point in the lattice the Metropolis algorithm will loop to (N for periodic, N-1 for Dirichlet)

// Thermalisation variables and counters //

bool thermalised;                   // Bool value to check if system is thermalised yet
int thermalisationCheck;      // Counter for how many consecutive times the thermalisation check has been passed, this is used to ensure the system is properly thermalised and not just fluctuating around a value for a few iterations
int thermalisationRepeats;      // Counter for how many sweeps have been performed during thermalisation
int thermalisationCount = 0;    // Counter for how many times the system has been considered thermalised, must reach thermalisationCheckLimit to be considered fully thermalised
const int thermalisationMeasures = 10000; // Number of measures to take during thermalisation to check the evolution of the ground state energy
const int thermalisationMeasureInterval = thermalisationMaximum / thermalisationMeasures; // Number of iterations between recording the evolution of the ground state energy during thermalisation, adjust to influence the detail of the thermalisation evolution output (decreasing this increases detail but also runtime and file size)
double currentE0, oldE0 = 0.0, E0Delta;   // Variables for tracking the evolution of the ground state energy estimate, useful for checking thermalisation 

// Measurement variables and counters //

int repeat;                         // Counts how many metropolis sweeps have been performed, useful for finding acceptance rate
int measureCount;				    // Counter for the number of measures taken
double E0Avg;				        // Counter for the ground state energy

// Acceptance rate variables and counters //

int acceptedMoves;				    // Counts how many of the metropolis updates were accepted 
int acceptedMovesPrevious;	        // Stores a previous value of acceptedMoves to find how many moves were accepted in some interval

// Window variables //

bool winRunning = false;			// Default status of the window
bool pathUpdated = false;			// Bool value to indicate an update for the window
const int delay = 0;                // Delay after updating window (milliseconds)
                                    // Do not change from zero unless creating screenshots or analysing some specific path

///// Function declarations /////

void chooseSystem();                                      // Allows user to choose which system to perform metropolis algorithm on

void metropolisRepeat(bool winOn, int BCs, int sys);      // Repeats the metropolis function "repeats" times for histogram production

void metropolis(bool winOn, int BCs, int sys, int test);            // Parent function for metropolisUpdate, runs loops to execute sufficient updates to the path

void metropolisUpdate(bool winOn, double (*potential)(double));  // Mechanics of the Metropolis algorithm such as checking whether a proposed update is accepted

void initialise(int sys);

void thermalise(bool winOn, double (*potentialDifferential)(double), double (*potential)(double));

void takeMeasures(std::vector<double> positions, double (*potentialDifferential)(double), double (*potential)(double), int test);

void createFiles(int BCs, int sys);

void writeVector(std::ofstream& file, std::ostringstream& buf, const std::vector<double>& data, const std::string& header);

void closeAllFiles();

void window(const std::vector<double>& path, bool& runningFlag); // Function for the visualisation of the path in time space


///// Short functions /////

static double euclideanActionElement(double xin, double xfin, double (*potential)(double));  // Calculates the Euclidean action for two points on a path
static double E0Calc(const std::vector<double>& path, double (*potentialDifferential)(double), double (*potential)(double)); // Calculates the ground state energy estimate for a given path using the virial theorem
static std::vector<double> twoPointCorrelator(const std::vector<double>& path); // Calculates the two point correlator 

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
            correlationTemp[n] += path[(t + n) % N] * path[t];
        }
    }

    for (int i = 0; i < N; i++) {
        correlationTemp[i] /= (N * measures);
    }

    return correlationTemp;
}