///// Simulating quantum systems with different potentials using Feynman path integral and Metropolis algorithm /////

#include "main.h"
#include "potentials.h"
#include "h5.h"
#include "random.h"
#include "window.h"
#include "csv.h" // Legacy csv writing functions, replaced by h5 files

int main() {    // Main function to run the Metropolis algorithm with user choice of boundary conditions and visualisation
    int __forceCRTManifestCUR = 0;
    std::cout << "Simulating quantum systems with different potentials using Feynman path integral and Metropolis algorithm" << std::endl;
    std::string choice;
	while (true) {
        chooseSystem();
		std::cout << "Enter choice: ";
        std::cin >> choice;
        if (choice == "1") {
            metropolisRepeat(false, "Periodic", "QHO");
		}
        else if (choice == "2") {
            metropolisRepeat(true, "Periodic", "QHO");
        }
        else if (choice == "3") {
            metropolisRepeat(false, "Dirichlet", "QHO");
        }
        else if (choice == "4") {
            metropolisRepeat(true, "Dirichlet", "QHO");
		}
        else if (choice == "5") {
            metropolisRepeat(false, "Periodic", "DWP");
        }
        else if (choice == "6") {
            metropolisRepeat(true, "Periodic", "DWP");
        }
        else if (choice == "7") {
            metropolisRepeat(false, "Dirichlet", "DWP");
        }
        else if (choice == "8") {
            metropolisRepeat(true, "Dirichlet", "DWP");
        }
        else if (choice == "0") {
            std::cout << "Exiting." << std::endl;
            break;
        }
        else {
            std::cout << "Invalid choice." << std::endl;
        }
    }
    return 0;
}

void chooseSystem() {  // Function to display user choices
    std::cout << "1: Perform Metropolis algorithm with periodic boundary conditions on the QHO system (without path visualisation)" << std::endl;
    std::cout << "2: Perform Metropolis algorithm with periodic boundary conditions on the QHO system (with path visualisation)" << std::endl;
	std::cout << "3: Perform Metropolis algorithm with Dirichlet boundary conditions on the QHO system (without path visualisation)" << std::endl;
	std::cout << "4: Perform Metropolis algorithm with Dirichlet boundary conditions on the QHO system (with path visualisation)" << std::endl;
    std::cout << "5: Perform Metropolis algorithm with periodic boundary conditions on the DWP system (without path visualisation)" << std::endl;
    std::cout << "6: Perform Metropolis algorithm with periodic boundary conditions on the DWP system (with path visualisation)" << std::endl;
    std::cout << "7: Perform Metropolis algorithm with Dirichlet boundary conditions on the DWP system (without path visualisation)" << std::endl;
    std::cout << "8: Perform Metropolis algorithm with Dirichlet boundary conditions on the DWP system (with path visualisation)" << std::endl;
    std::cout << "0: Exit" << std::endl;
    return;
}

void metropolisRepeat(bool winOn, std::string boundary, std::string system) {
    if (boundary == "Periodic") {         // Periodic boundary conditions
        start = 0;
        end = N;
    }
    else if (boundary == "Dirichlet") {    // Dirichlet boundary conditions
        start = 1;
        end = N - 1;        // Effectively fixes the endpoints to 0, as they are never updated
    }
    else {
        std::cerr << "Invalid boundary condition choice." << std::endl;
        return;
    }

    if (system == "FP") {     // Free particle
        potential = FP::potential;
        potentialDifferential = FP::potentialDifferential;
    }
    else if (system == "QHO") { // Quantum harmonic oscillator
        potential = QHO::potential;
        potentialDifferential = QHO::potentialDifferential;
    }
    else if (system == "DWP") { // Double-well potential
        potential = DWP::potential;
        potentialDifferential = DWP::potentialDifferential;
    }
    else {
        std::cerr << "Invalid system choice." << std::endl;
        return;
    }

    winRunning = winOn; // Sets winRunning flag to true if user wanted visualisation
    std::thread windowThread(window, std::ref(positions), std::ref(winRunning)); // Instantiates the window thread, passing the path and winRunning flag by reference

    for (int repeat = 0; repeat < repeats; repeat++) {
        metropolis(winOn, boundary, system, repeat);
    }
    winRunning = false; // Sets winRunning flag to false to close the window thread
    windowThread.join();
    
    openFiles(boundary, system);
	writeFiles(boundary, system);         // Legacy csv writing functions, replaced by h5 files
    closeAllFiles();

    writeData(boundary, system);            // Writes all data to a single h5 file, separated into groups
    return;
}

void metropolis(bool winOn, std::string boundary, std::string system, int repeat) { // Metropolis function which gets called initially, then calls other functions to perform the algorithm
	initialise(boundary, system);   // Sets initial parameters and path

    thermalise(winOn, potentialDifferential, potential);

    takeMeasures(positions, potentialDifferential, potential, repeat);
    std::cout << "Thermalisation " << repeat << " complete after " << thermalisationSweeps << " sweeps" << std::endl;
	double percentMeasured = 0;
    bool percentChanged = true;
    while (measureCount < measures) {
        metropolisUpdate(winOn, potential);
        sweep++;
        if (remainder(sweep, decorrelation) == 0) {
            takeMeasures(positions, potentialDifferential, potential, repeat);
            while ((100 * measureCount) / measures > percentMeasured) {
                percentMeasured++;
                percentChanged = true;
            }
            if (percentChanged == true) {
                //std::cout << percentMeasured << "% done." << std::endl;
            }
            percentChanged = false;
        }
        if (sweep % accRateInterval == 0) {
            acceptanceRate.push_back((double)(acceptedMoves - acceptedMovesPrevious) / (N * (double)accRateInterval));
            acceptedMovesPrevious = acceptedMoves;
        }
    }
    double accRate = (double)acceptedMoves / ((double)sweep * N);
    double E0 = E0Avg / measures;
    double E1 = 0;
    int d = 0;
    while (true) {
        double arg = G[(int)floor(6) + d] / G[(int)floor(6) + 1 + d];
        if (arg > 1.0) {
            E1 = E0 + (log(arg) / a);
            break;
        }
        else {
            d++;    // Arguement was negative, due to noise??
        }
    }
    E0Vec.push_back((double)E0);
    E1Vec.push_back((double)E1);
    accRateVec.push_back((double)accRate);
    return;
}

void metropolisUpdate(bool winOn, double (*potential)(double)) {
    double newPosition;
    for (int i = start; i < end; i++) {
        double y = rfRange(-1, 1);
        newPosition = positions[i] + epsilon * y; // Incrementing one position by float between -epsilon and +epsilon
        double actionDelta = euclideanActionElement(newPosition, positions[(i + 1) % N], potential) -
            euclideanActionElement(positions[i], positions[(i + 1) % N], potential) +
            euclideanActionElement(positions[(i - 1 + N) % N], newPosition, potential) -
            euclideanActionElement(positions[(i - 1 + N) % N], positions[i], potential);
        double positionDelta = newPosition - positions[i];
        if (actionDelta < 0) { // Favorable move
            positions[i] = newPosition; // Updates position stored (no need to store rejected position)
            acceptedMoves++;
        }
        else {  // Unfavorable move
            auto& rng = globalRng();
            float r = static_cast<float>(rfRange(0, 1));
            if (r < exp(-actionDelta)) { // Accept with probability exp(-deltaS)
                positions[i] = newPosition; // Updates position stored (no need to store rejected position)
                acceptedMoves++;
            }
        }
    }
    if (winOn == true) {
        std::this_thread::sleep_for(std::chrono::milliseconds(delay));
    }
    return;
}

void initialise(std::string boundary, std::string system) {   // Initialises variables and path
    positions = std::vector<double>(N, 0.0);
    sweep = 0;
    acceptedMoves = 0;
    measureCount = 0;
    E0Avg = 0;
    thermalised = false;
	oldE0 = 0;
    thermalisationCheck = 0;
    //if (system == "DWP") // DWP requires a different initial (cold) path since the potential wells are not centered around 0
    //{
    //    positions = std::vector<double>(N, wellCentres);
    //}
    return;
}

void thermalise(bool winOn, double (*potentialDifferential)(double), double (*potential)(double)) { // Thermalisation function to reach equilibrium before measurements
    while (thermalised == false) {
        metropolisUpdate(winOn, potential);
        currentE0 = E0Calc(positions, potentialDifferential, potential); // Ground state energy measurement during thermalisation
        E0Delta = currentE0 - oldE0;
        if (sweep % thermalisationMeasureInterval == 0) { // Store E0 during thermalisation for plotting purposes (not necessary for the algorithm itself)
            E0Thermalising.push_back((double)currentE0);
            thermalisationCount++;
        }
        if (sweep % accRateInterval == 0) {
            acceptanceRate.push_back((double)(acceptedMoves - acceptedMovesPrevious) / (double)accRateInterval);
            acceptedMovesPrevious = acceptedMoves;  // For measuring acceptance rate over sweeps
        }
        sweep++;
        if (E0Delta < thermalisationConstant) {
			thermalisationCheck++;  // If E0 is not changing much, increment thermalisation check counter
            if (thermalisationCheck == thermalisationCheckLimit) { // If this is the first time E0 has not changed significantly, store the current sweep count as a potential thermalisation point (if E0 remains stable for long enough, we will assume thermalisation is complete from this point)
                thermalisationCheck = thermalisationCheckLimit * thermalisationMinimum; // Prevents thermalisation from being marked as incomplete after it has completed once
            }
        }
        else {
			thermalisationCheck -= thermalisationDecrement; // If E0 is changing significantly, lower thermalisation check counter
            if (thermalisationCheck < 0) {
                thermalisationCheck = 0; // Prevents thermalisation check counter from going negative
            }
        }
		if (thermalisationCheck >= thermalisationCheckLimit && sweep > thermalisationMinimum) { // If E0 has not changed significantly for a certain number of checks, we can assume thermalisation is complete
            thermalisationSweeps = sweep;
            thermalised = true;
        }
        if (sweep == thermalisationMaximum) {
            std::cout << "Thermalisation failed to converge after " << thermalisationMaximum << " sweeps, proceeding with measurements anyway." << std::endl;
            thermalisationSweeps  = thermalisationMaximum;
            thermalised = true;
        }
		oldE0 = currentE0;
    }
    return;
}

void takeMeasures(std::vector<double> positions, double (*potentialDifferential)(double), double (*potential)(double), int repeat) {    // Takes all measurements at current path state in one function

    double E0temp = E0Calc(positions, potentialDifferential, potential);
    E0Avg += E0temp;
    E0Evolution.push_back((double)E0temp);

    for (int j = 0; j < N; j++) {
        psi[j + N * (measureCount + (repeat * measures))] = positions[j];
    }

    std::vector<double> corr = twoPointCorrelator(positions);
    for (int j = 0; j < N; j++) {
        G[j + (repeat * N)] += corr[j];
    }
    
    measureCount++;
    return;
}