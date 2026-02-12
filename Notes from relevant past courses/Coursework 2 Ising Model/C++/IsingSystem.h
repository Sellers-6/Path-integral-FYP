//
// IsingSystem.h
//

#pragma once

#define _USE_MATH_DEFINES
#include <gl/freeglut.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <random>
#include <string>
#include <sstream>
#include "Window.h"
#include "rnd.h"

using namespace std;

class IsingSystem {
private:
	Window* win;  // pointer to the window in which the system is running
	// (this is not actually used but we keep it, just in case)

	rnd rgen;

	ofstream logfile;
	char q1[9] = "IMQ1.csv";
	char q2[9] = "IMQ2.csv";
	char q3[9] = "IMQ3.csv";
	char q4[9] = "IMQ4.csv";
	char q5[9] = "IMQ5.csv";
	char q6[9] = "IMQ6.csv";

	static const int gridSize = 250;
	static const int gridSize3D = 40;

	int** grid;  // this will be a 2d array that stores whether each site is occupied
	int*** grid3D; // this will be a 3d array that stores whether each site is occupied for the 3D IM
	int isActive;	// this variable is 1 if the system is running, 0 for paused
	int slowNotFast;	// this variable is 1 if we are in "slow" mode, 0 in fast mode
	double inverseTemperatureBeta = 0.25;	// this is the inverse temperature parameter beta

public:
	// constructor
	IsingSystem();
	// destructor
	~IsingSystem();

	void readCsv(); // finds the specific heat capacity to modulate n and n_0
	const int maxsize = 202;
	int readcount = 0;
	int location = 0;
	double temperatures[202] = { 0 };
	double energies[202] = { 0 };
	double heatcapacities[202] = { 0 };
	double magnetisms[202] = { 0 };
	string line;

	int returnGridSize() { return gridSize; }

	void createWindow(Window* set_win); // function for displaying the window/turning on openGL

	bool glutOn = false; // flag for whether the glutMainLoop is active

	// increase or decrease the temperature (arbitrary increments)
	void Hotter() { inverseTemperatureBeta -= 0.05; }
	void Colder() { inverseTemperatureBeta += 0.05; }

	void setTemperature(double TT) { inverseTemperatureBeta = 1.0 / TT; } // Useful for taking steady inrememnts in temperature rather than beta

	// NOTE! the grid is not accessed directly, only via these functions
	void setGrid(int pos[], int val);
	int readGrid(int pos[]);
	void flipSpin(int pos[]);

	// "user interface" for running and slow variables
	int  isRunning() { return isActive; }
	void setRunning() { isActive = 1; }
	void pauseRunning() { isActive = 0; }
	int isSlow() { return slowNotFast; }
	void setSlow() { slowNotFast = 1; }
	void setFast() { slowNotFast = 0; }

	// these functions are used to update the system
	double computeLocalField(int pos[]);
	void attemptSpinFlip();
	void MCsweep();

	void setSeed(int s) { rgen.setSeed(s); }; // sets the seed 
	void DrawSquares();	// draws the system as squares
	void Reset();	// reset temperature to default value and re-initialise the grid
	void resetGrid(); // re-initialise the grid only

	// this returns "setpos" as the neighbour of "pos" in a direction dependant on "val", val integer between 0 and 3
	void setPosNeighbour(int setpos[], int pos[], int val, int distance);

	// These are the new functions created to analyse the Ising model;
	double computeMagnetisation();
	double computeLocalInteraction(int pos[], int distance);
	double computeEnergy();
	int computeCorrelation(int distance);
	double computeCorrelationLarge(int distance);
	int findn(int temperature);
	int findn_0(int temperature);

	// These are the new functions created to analyse the 3D Ising model;
	int return3DGridSize() { return gridSize3D; }
	void setGrid3D(int pos[], int val);
	int readGrid3D(int pos[]);
	void flipSpin3D(int pos[]);
	double computeLocalField3D(int pos[]);
	void attemptSpinFlip3D();
	void MCsweep3D();
	void Reset3D(); 
	void resetGrid3D(); 
	void setPosNeighbour3D(int setpos[], int pos[], int val, int distance);
	double computeMagnetisation3D();
	double computeLocalInteraction3D(int pos[], int distance);
	double computeEnergy3D();
	double computeCorrelationLarge3D(int distance);
	const int distanceRepeats3D = gridSize3D - 1;
	double correlation3D[gridSize3D];

	int sweepNumber;	// counts the number of sweeps performed

	int measurements;	// Amount of times we take a measure (m)
	int repeats;		// Repeats for questions 1 and 5 which replaces measurements
	int waitingSweeps;	// Sweeps before starting measurements (n_0)
	int spacingSweeps;	// Sweeps between computing properties of the IM (n)
	int sweeps;		// Number of sweeps performed in question 1 and 5

	void csvCreation(int qNum);
	int progress = 0;	// For the progress bar

	int betaRepeats;
	double initialBeta;
	double finalBeta;

	double temperature;
	int tempRepeats;
	double initialTemp;
	double finalTemp;

	const int distanceRepeats = gridSize - 1;
	double correlation[gridSize];
	int leftSpin;
	int rightSpin;

	void question1(); // Where you study how quickly equilibrium is reached	

	void question2(); // Where you calculate the results to compare against theory

	void question3(); // Where you calculate the correlation function

	void question4(); // Where you study how quickly equilibrium is reached	(3D)

	void question5(); // Where you calculate the results for 3D

	void question6(); // Where you calculate the correlation function (3D)
};