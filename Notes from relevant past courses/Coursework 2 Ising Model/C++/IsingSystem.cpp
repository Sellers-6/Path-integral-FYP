//
//  IsingSystem.cpp
//

#include "IsingSystem.h"

// colours
namespace colours {
	// note the f's here avoid warnings by telling C++ to read the numbers as floats and not doubles 
	GLfloat blue[] = { 0.1f, 0.3f, 0.9f, 1.0f };   // blue
	GLfloat red[] = { 1.0f, 0.2f, 0.1f, 0.2f };   // red
	GLfloat green[] = { 0.3f, 0.6f, 0.3f, 1.0f };   // green
}

// constructor
IsingSystem::IsingSystem() {
	// Set the initial conditions here!
	slowNotFast = 1;
	isActive = 0;

	// Allocate memory for the grid, remember to free the memory in destructor
	//   the point here is that each row of the grid is an array
	//   the grid itself is a an array of pointers, one for each row

	grid = new int* [gridSize];	// Here we allocate the array of pointers

	for (int i = 0; i < gridSize; i++) {
		grid[i] = new int[gridSize];		// Now allocate the indidual rows
	}

	grid3D = new int** [gridSize3D];

	for (int i = 0; i < gridSize3D; i++) {
		grid3D[i] = new int* [gridSize3D];  // Allocate the 2D layer
		for (int j = 0; j < gridSize3D; j++) {
			grid3D[i][j] = new int[gridSize3D]; // Allocate the 3D layer
		}
	}

	Reset();	// this sets the temperatre and initialises the spins grid
}

void IsingSystem::readCsv() {
	ifstream fin;
	fin.open("exactisingdatacsv.csv");

	getline(fin, line);

	while (getline(fin, line) && readcount < maxsize) {
		location = line.find(",");
		temperatures[readcount] = stod(line.substr(0, location));
		line = line.substr(location + 1, line.length());

		location = line.find(",");
		energies[readcount] = stod(line.substr(0, location));
		line = line.substr(location + 1, line.length());

		location = line.find(",");
		heatcapacities[readcount] = stod(line.substr(0, location));
		magnetisms[readcount] = stod(line.substr(location + 1, line.length()));

		readcount++;
	}
	
	//for (int i = 0; i < maxsize - 1; i++) {
	//	cout << heatcapacities[i] << endl;
	//} // Testing that the heat capacities were read correctly
	fin.close();
}

void IsingSystem::createWindow(Window* set_win) {
	std::cout << "Creating window, gridSize " << gridSize << endl;
	win = set_win;
}

void IsingSystem::Reset() {
	double initialTemp = 2.3;
	sweepNumber = 0;

	setTemperature(initialTemp);

	// set the grid to -1
	for (int i = 0; i < gridSize; i++) {
		for (int j = 0; j < gridSize; j++) {
			// position is (i,j)
			int pos[2] = { i, j };
			// set this spin to state -1
			setGrid(pos, -1);
		}
	}
}

void IsingSystem::resetGrid() {
	// set the grid to -1
	for (int i = 0; i < gridSize; i++) {
		for (int j = 0; j < gridSize; j++) {
			// position is (i,j)
			int pos[2] = { i, j };
			// set this spin to state -1
			setGrid(pos, -1);
		}
	}
	sweepNumber = 0;
}

// destructor
IsingSystem::~IsingSystem() {
	// Close the file (if open)
	if (logfile.is_open())
		logfile.close();

	// Delete the window
	if (glutOn == true) {
		delete win;
	}

	// Delete the grid
	// First we delete the individual rows
	for (int i = 0; i < gridSize; i++)
		delete[] grid[i];
	// Finally delete the array of pointers
	delete[] grid;

	// Delete the 3D grid
	// First we delete the individual rows and columns
	for (int i = 0; i < gridSize3D; i++) {
		for (int j = 0; j < gridSize3D; j++) {
			delete[] grid3D[i][j];
		}
		delete[] grid3D[i];
	}
	// Finally delete the array of pointers
	delete[] grid3D;
}

// this draws the 2D system 
void IsingSystem::DrawSquares() {

	double drawScale = 2.0 / (gridSize * 1.1);

	// draw the particles
	double halfSize = 0.5;
	int halfGrid = gridSize / 2;
	for (int x = 0; x < gridSize; x++) {
		for (int y = 0; y < gridSize; y++) {

			double vec[2];
			vec[0] = x - halfGrid;
			vec[1] = y - halfGrid;

			// openGL magic
			glPushMatrix();
			// choose a color
			if (grid[x][y] == -1)
				glColor4fv(colours::green);
			else
				glColor4fv(colours::blue);
			// draw a rectangle for the particle
			glRectd(drawScale * (vec[0] - halfSize),
				drawScale * (vec[1] - halfSize),
				drawScale * (vec[0] + halfSize),
				drawScale * (vec[1] + halfSize));
			// openGL magic
			glPopMatrix();
		}
	}
	// print some information (at top left)
	// this ostringstream is a way to make a string with numbers and words (similar to cout << ... )
	ostringstream str;
	str << "beta " << inverseTemperatureBeta << " size " << gridSize;
	win->displayString(str, -0.9, 0.94, colours::red);
}

// attempt N spin flips, where N is the number of spins
void IsingSystem::MCsweep() {
	for (int i = 0; i < gridSize * gridSize; i++) {
		attemptSpinFlip();
	}
	sweepNumber++;
}

// here we attempt to flip a spin and accept/reject with Metropolis rule
void IsingSystem::attemptSpinFlip() {
	int pos[2];

	// random site
	pos[0] = rgen.randomInt(gridSize);
	pos[1] = rgen.randomInt(gridSize);

	double hloc = computeLocalField(pos);

	double dE = 2.0 * hloc * readGrid(pos);
	if (dE < 0) {
		flipSpin(pos);
		//cout << "flipped by default at: " << pos[0] << pos[1] << endl;
	}
	else if (rgen.random01() < exp(-dE)) {
		//cout << "flipped by chance at: " << pos[0] << pos[1] << endl;
		flipSpin(pos);
	}
	/*else {
		cout << "didn't flip at: " << pos[0] << pos[1] << endl;
	}*/
}

// NOTE: this returns the local field *divided by the temperature* (dimensionless quantity)
double IsingSystem::computeLocalField(int pos[]) {
	double result = 0.0;
	for (int i = 0; i < 4; i++) {
		int nborPos[2];
		setPosNeighbour(nborPos, pos, i, 1);
		result += readGrid(nborPos);
	}
	result *= inverseTemperatureBeta;
	return result;
}

double IsingSystem::computeLocalInteraction(int pos[], int distance) {
	double result = 0.0;
	int nborPos[2];
	setPosNeighbour(nborPos, pos, 0, distance);
	result += readGrid(nborPos) * readGrid(pos);
	setPosNeighbour(nborPos, pos, 2, distance);
	result += readGrid(nborPos) * readGrid(pos);
	return result;
}

// set the value of a grid cell for a particular position
void IsingSystem::setGrid(int pos[], int val) {
	grid[pos[0]][pos[1]] = val;
}

// read the grid cell for a given position
int IsingSystem::readGrid(int pos[]) {
	return grid[pos[0]][pos[1]];
}

// flip the grid cell for a given position
void IsingSystem::flipSpin(int pos[]) {
	grid[pos[0]][pos[1]] = -grid[pos[0]][pos[1]];
}

double IsingSystem::computeMagnetisation() {
	double totalSpin = 0;
	int pos[2];
	for (int i = 0; i < gridSize; i++) {
		for (int j = 0; j < gridSize; j++) {
			pos[0] = i;
			pos[1] = j;
			totalSpin += readGrid(pos);
		}
	}
	return totalSpin / (gridSize * gridSize);
}

double IsingSystem::computeEnergy() {
	double totalEnergy = 0;
	int pos[2];
	for (int i = 0; i < gridSize; i++) {
		for (int j = 0; j < gridSize; j++) {
			pos[0] = i;
			pos[1] = j;
			totalEnergy += computeLocalInteraction(pos, 1);
		}
	}
	return -totalEnergy / (gridSize * gridSize);
}

int IsingSystem::computeCorrelation(int distance) {
	int posLeft[2];
	posLeft[0] = 6;
	posLeft[1] = 6;
	int posRight[2];
	posRight[0] = 6 + distance;
	posRight[1] = 6;
	leftSpin = readGrid(posLeft); // These can be changed to any spin pair in the grid
	rightSpin = readGrid(posRight);
	return leftSpin * rightSpin;
}

double IsingSystem::computeCorrelationLarge(int distance) {
	double totalCorrelation = 0;
	int pos[2];
	for (int i = 0; i < gridSize; i++) {
		for (int j = 0; j < gridSize; j++) {
			pos[0] = i;
			pos[1] = j;
			totalCorrelation += computeLocalInteraction(pos, distance);
		}
	}
	return totalCorrelation / (gridSize * gridSize * 2); // dividing by 2 becuase both up and right neighbours were considered
}

int IsingSystem::findn(int temperature) {
	return 3000 * heatcapacities[temperature] + 400;
}

int IsingSystem::findn_0(int temperature) {
	return 300 * heatcapacities[temperature] + 20;
}

// send back the position of a neighbour of a given grid cell
// NOTE: we take care of periodic boundary conditions, also positions are integers now not doubles
void IsingSystem::setPosNeighbour(int setpos[], int pos[], int val, int distance) {
	switch (val) {
	case 0:
		setpos[0] = (pos[0] + distance) % gridSize;
		setpos[1] = pos[1];
		break;
	case 1:
		setpos[0] = (pos[0] - distance + distance * gridSize) % gridSize;
		setpos[1] = pos[1];
		break;
	case 2:
		setpos[0] = pos[0];
		setpos[1] = (pos[1] + distance) % gridSize;
		break;
	case 3:
		setpos[0] = pos[0];
		setpos[1] = (pos[1] - distance + distance * gridSize) % gridSize;
		break;
	}
}

//
// 3D Ising model section:
//

void IsingSystem::Reset3D() {
	double initialTemp = 4.0;
	sweepNumber = 0;

	setTemperature(initialTemp);

	// set the grid to -1
	for (int i = 0; i < gridSize3D; i++) {
		for (int j = 0; j < gridSize3D; j++) {
			for (int k = 0; k < gridSize3D; k++) {
				// position is (i,j)
				int pos[3] = { i, j, k };
				// set this spin to state -1
				setGrid3D(pos, -1);
			}
		}
	}
}

void IsingSystem::resetGrid3D() {
	// set the grid to -1
	for (int i = 0; i < gridSize3D; i++) {
		for (int j = 0; j < gridSize3D; j++) {
			for (int k = 0; k < gridSize3D; k++) {
				// position is (i,j)
				int pos[3] = { i, j, k };
				// set this spin to state -1
				setGrid3D(pos, -1);
			}
		}
	}
	sweepNumber = 0;
}

void IsingSystem::MCsweep3D() {
	for (int i = 0; i < gridSize3D * gridSize3D * gridSize3D; i++) {
		attemptSpinFlip3D();
	}
	sweepNumber++;
}

void IsingSystem::attemptSpinFlip3D() {
	int pos[3];

	// random site
	pos[0] = rgen.randomInt(gridSize3D);
	pos[1] = rgen.randomInt(gridSize3D);
	pos[2] = rgen.randomInt(gridSize3D);

	double hloc = computeLocalField3D(pos);

	double dE = 2.0 * hloc * readGrid3D(pos);
	if (dE < 0)
		flipSpin3D(pos);
	else if (rgen.random01() < exp(-dE))
		flipSpin3D(pos);
}

double IsingSystem::computeLocalField3D(int pos[]) {
	double result = 0.0;
	for (int i = 0; i < 6; i++) {
		int nborPos[3];
		setPosNeighbour3D(nborPos, pos, i, 1);
		result += readGrid3D(nborPos);
	}
	result *= inverseTemperatureBeta;
	return result;
}

double IsingSystem::computeLocalInteraction3D(int pos[], int distance) {
	double result = 0.0;
	int nborPos[3];
	setPosNeighbour3D(nborPos, pos, 0, distance);
	result += readGrid3D(nborPos) * readGrid3D(pos);
	setPosNeighbour3D(nborPos, pos, 2, distance);
	result += readGrid3D(nborPos) * readGrid3D(pos);
	setPosNeighbour3D(nborPos, pos, 4, distance);
	result += readGrid3D(nborPos) * readGrid3D(pos);
	return result;
}

void IsingSystem::setGrid3D(int pos[], int val) {
	grid3D[pos[0]][pos[1]][pos[2]] = val;
}

int IsingSystem::readGrid3D(int pos[]) {
	return grid3D[pos[0]][pos[1]][pos[2]];
}

void IsingSystem::flipSpin3D(int pos[]) {
	grid3D[pos[0]][pos[1]][pos[2]] = -grid3D[pos[0]][pos[1]][pos[2]];
}

double IsingSystem::computeMagnetisation3D() {
	double totalSpin = 0;
	int pos[3];
	for (int i = 0; i < gridSize3D; i++) {
		for (int j = 0; j < gridSize3D; j++) {
			for (int k = 0; k < gridSize3D; k++) {
				pos[0] = i;
				pos[1] = j;
				pos[2] = k;
				totalSpin += readGrid3D(pos);
			}
		}
	}
	return totalSpin / (gridSize3D * gridSize3D * gridSize3D);
}

double IsingSystem::computeEnergy3D() {
	double totalEnergy = 0;
	int pos[3];
	for (int i = 0; i < gridSize3D; i++) {
		for (int j = 0; j < gridSize3D; j++) {
			for (int k = 0; k < gridSize3D; k++) {
				pos[0] = i;
				pos[1] = j;
				pos[2] = k;
				totalEnergy += computeLocalInteraction3D(pos, 1);
			}
		}
	}
	return -totalEnergy / (gridSize3D * gridSize3D * gridSize3D);
}

void IsingSystem::setPosNeighbour3D(int setpos[], int pos[], int val, int distance) {
	switch (val) {
	case 0:
		setpos[0] = (pos[0] + distance) % gridSize3D;
		setpos[1] = pos[1];
		setpos[2] = pos[2];
		break;
	case 1:
		setpos[0] = (pos[0] - distance + distance * gridSize3D) % gridSize3D;
		setpos[1] = pos[1];
		setpos[2] = pos[2];
		break;
	case 2:
		setpos[0] = pos[0];
		setpos[1] = (pos[1] + distance) % gridSize3D;
		setpos[2] = pos[2];
		break;
	case 3:
		setpos[0] = pos[0];
		setpos[1] = (pos[1] - distance + distance * gridSize3D) % gridSize3D;
		setpos[2] = pos[2];
		break;
	case 4:
		setpos[0] = pos[0];
		setpos[1] = pos[1];
		setpos[2] = (pos[2] + distance) % gridSize3D;
		break;
	case 5:
		setpos[0] = pos[0];
		setpos[1] = pos[1];
		setpos[2] = (pos[2] - distance + distance * gridSize3D) % gridSize3D;
		break;
	}
}

double IsingSystem::computeCorrelationLarge3D(int distance) {
	double totalCorrelation = 0;
	int pos[3];
	for (int i = 0; i < gridSize3D; i++) {
		for (int j = 0; j < gridSize3D; j++) {
			for (int k = 0; k < gridSize3D; k++) {
				pos[0] = i;
				pos[1] = j;
				pos[2] = k;
				totalCorrelation += computeLocalInteraction3D(pos, distance);
			}
		}
	}
	return totalCorrelation / (gridSize3D * gridSize3D * gridSize3D * 3);
}

//
// End of 3D section
//

void IsingSystem::csvCreation(int qNum) {
	Reset();
	Reset3D();
	std::cout << "Creating csv for question " << qNum << ", progress(%) :" << endl;
	std::cout << "|0        10        20        30        40        50        60        70        80        90      100|" << endl << "|";
	switch (qNum) {
	case 1:
		logfile.open(q1);
		logfile << "Sweep" << "," << "Repeat" << "," << "Temperature" << "," << "Magnetisation" << "," << "Energy" << endl;
		initialTemp = 1.0;
		finalTemp = 3.0;
		sweeps = 1000;
		repeats = 30;
		tempRepeats = 8;
		question1();
		break;
	case 2:
		logfile.open(q2);
		logfile << "Measurement" << "," << "Temperature" << "," << "Magnetisation" << "," << "Energy" << "," << "Gridsize" << endl;
		initialTemp = 1.0;
		finalTemp = 3.0;
		tempRepeats = 200;
		measurements = 30;
		question2();
		break;
	case 3:
		logfile.open(q3);
		logfile << "Temperature" << "," << "Distance" << "," << "Correlation" << endl;
		initialTemp = 1.0;
		finalTemp = 3.0;
		tempRepeats = 8;
		measurements = 30;
		question3();
		break;
	case 4:
		logfile.open(q4);
		logfile << "Sweep" << "," << "Repeat" << "," << "Temperature" << "," << "Magnetisation" << "," << "Energy" << endl;
		initialTemp = 2.5;
		finalTemp = 5.5;
		sweeps = 1000;
		repeats = 30;
		tempRepeats = 6;
		question4();
	case 5:
		logfile.open(q5);
		logfile << "Measurement" << "," << "Temperature" << "," << "Magnetisation" << "," << "Energy" << endl;
		initialTemp = 2.5;
		finalTemp = 5.5;
		tempRepeats = 300;
		measurements = 30;
		waitingSweeps = 2000;
		spacingSweeps = 50;
		question5();
		break;
	case 6:
		logfile.open(q6);
		logfile << "Temperature" << "," << "Distance" << "," << "Correlation" << endl;
		initialTemp = 2.5;
		finalTemp = 5.5;
		tempRepeats = 6;
		measurements = 30;
		waitingSweeps = 2000;
		spacingSweeps = 50;
		question6();
		break;
	}
	std::cout << "|" << endl << "Completed csv creation." << endl;
	progress = 0;
	logfile.close();
}

void IsingSystem::question1() {
	double energy;
	double magnetisation;
	for (int k = 0; k <= tempRepeats; k++) {
		setTemperature(initialTemp + ((finalTemp - initialTemp) / tempRepeats) * k);
		temperature = 1 / inverseTemperatureBeta;
		for (int j = 0; j < repeats; j++) {
			for (int i = 0; i < sweeps; i++) {
				magnetisation = computeMagnetisation();
				energy = computeEnergy();
				logfile << i << "," << j << "," << temperature << "," << magnetisation << "," << energy << endl;
				if ((i + j * sweeps + k * sweeps * repeats) * 100 > progress * sweeps * repeats * (tempRepeats + 1)) {
					std::cout << ".";
					progress++;
				}
				MCsweep();
			}
			resetGrid();
		}
	}
}

void IsingSystem::question2() {
	double energy;
	double magnetisation;
	for (int k = 0; k <= tempRepeats; k++) {
		setTemperature(initialTemp + ((finalTemp - initialTemp) / tempRepeats) * k);
		temperature = 1 / inverseTemperatureBeta;
		waitingSweeps = findn_0(k);
		spacingSweeps = findn(k);
		for (int i = 0; i < waitingSweeps; i++) {
			MCsweep();
		}
		for (int i = 0; i < measurements; i++) {
			magnetisation = computeMagnetisation();
			energy = computeEnergy();
			logfile << i << "," << temperature << "," << magnetisation << "," << energy << "," << gridSize << endl;
			if ((i + k * measurements) * 100 > progress * (tempRepeats + 1) * measurements) {
				std::cout << ".";
				progress++;
			}
			for (int j = 0; j < spacingSweeps; j++) {
				MCsweep();
			}
		}
		resetGrid();
	}
}

void IsingSystem::question3() {
	for (int t = 0; t <= tempRepeats; t++) {
		setTemperature(initialTemp + ((finalTemp - initialTemp) / tempRepeats) * t);
		temperature = 1 / inverseTemperatureBeta;
		waitingSweeps = findn_0((temperature - initialTemp) * 100);
		spacingSweeps = findn((temperature - initialTemp) * 100);
		for (int i = 0; i < waitingSweeps; i++) {
			MCsweep();
		}
		for (int k = 1; k <= distanceRepeats; k++) {
			correlation[k] = 0;
		}
		for (int i = 0; i < measurements; i++) {
			for (int k = 1; k <= distanceRepeats; k++) {
				correlation[k] += computeCorrelationLarge(k);
				if (((k - 1) + i * distanceRepeats + t * measurements * distanceRepeats) * 100 > progress * (tempRepeats + 1) * distanceRepeats * measurements) {
					std::cout << ".";
					progress++;
				}
			}
			for (int j = 0; j < spacingSweeps; j++) {
				MCsweep();
			}
		}
		for (int k = 1; k <= distanceRepeats; k++) {
			correlation[k] = correlation[k] / measurements;
			logfile << temperature << "," << k << "," << correlation[k] << endl;
		}
		resetGrid();
	}
}

void IsingSystem::question4() {
	double energy;
	double magnetisation;
	for (int k = 0; k <= tempRepeats; k++) {
		setTemperature(initialTemp + ((finalTemp - initialTemp) / tempRepeats) * k);
		temperature = 1 / inverseTemperatureBeta;
		for (int j = 0; j < repeats; j++) {
			for (int i = 0; i < sweeps; i++) {
				magnetisation = computeMagnetisation3D();
				energy = computeEnergy3D();
				logfile << i << "," << j << "," << temperature << "," << magnetisation << "," << energy << endl;
				if ((i + j * sweeps + k * sweeps * repeats) * 100 > progress * sweeps * repeats * (tempRepeats + 1)) {
					std::cout << ".";
					progress++;
				}
				MCsweep3D();
			}
			resetGrid3D();
		}
	}
}

void IsingSystem::question5() {
	double energy;
	double magnetisation;
	for (int k = 0; k <= tempRepeats; k++) {
		setTemperature(initialTemp + ((finalTemp - initialTemp) / tempRepeats) * k);
		temperature = 1 / inverseTemperatureBeta;
		for (int i = 0; i < waitingSweeps; i++) {
			MCsweep3D();
		}
		for (int i = 0; i < measurements; i++) {
			magnetisation = computeMagnetisation3D();
			energy = computeEnergy3D();
			logfile << i << "," << temperature << "," << magnetisation << "," << energy << endl;
			if ((i + k * measurements) * 100 > progress * (tempRepeats + 1) * measurements) {
				std::cout << ".";
				progress++;
			}
			for (int j = 0; j < spacingSweeps; j++) {
				MCsweep3D();
			}
		}
		resetGrid3D();
	}
}

void IsingSystem::question6() {
	for (int t = 0; t <= tempRepeats; t++) {
		setTemperature(initialTemp + ((finalTemp - initialTemp) / tempRepeats) * t);
		temperature = 1 / inverseTemperatureBeta;
		for (int i = 0; i < waitingSweeps; i++) {
			MCsweep3D();
		}
		for (int k = 1; k <= distanceRepeats3D; k++) {
			correlation3D[k] = 0;
		}
		for (int i = 0; i < measurements; i++) {
			for (int k = 1; k <= distanceRepeats3D; k++) {
				correlation3D[k] += computeCorrelationLarge3D(k);
				if (((k - 1) + i * measurements + t * measurements * distanceRepeats3D) * 100 > progress * (tempRepeats + 1) * distanceRepeats3D * measurements) {
					std::cout << ".";
					progress++;
				}
			}
			for (int j = 0; j < spacingSweeps; j++) {
				MCsweep3D();
			}
		}
		for (int k = 1; k <= distanceRepeats3D; k++) {
			correlation3D[k] = correlation3D[k] / measurements;
			logfile << temperature << "," << k << "," << correlation3D[k] << endl;
		}
		resetGrid();
	}
}