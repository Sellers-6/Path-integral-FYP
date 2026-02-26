///// Simulating quantum systems with different potentials using Feynman path integral and Metropolis algorithm /////
// Throughout the program, hbar is assumed to be 1 (reduced planck's constant) //

#include "main.h"

int main() {    // Accepts user choice of boundary conditions, system type, and window visualisation
    std::cout << "Simulating quantum systems with different potentials using Feynman path integral and Metropolis algorithm" << std::endl;
    std::string choice;
	chooseSystem(); // Print the options to the user
	while (true) {
		std::cout << "Enter choice: ";
        std::cin >> choice;
        if (choice == "1") { metropolisRepeat(false, "Periodic", "QHO"); }
        else if (choice == "2") { metropolisRepeat(true, "Periodic", "QHO"); }
        else if (choice == "3") { metropolisRepeat(false, "Dirichlet", "QHO"); }
        else if (choice == "4") { metropolisRepeat(true, "Dirichlet", "QHO"); }
        else if (choice == "5") { metropolisRepeat(false, "Periodic", "AHO"); }
        else if (choice == "6") { metropolisRepeat(true, "Periodic", "AHO"); }
        else if (choice == "7") { metropolisRepeat(false, "Dirichlet", "AHO"); }
        else if (choice == "8") { metropolisRepeat(true, "Dirichlet", "AHO"); }
        else if (choice == "9") { metropolisRepeat(false, "Periodic", "DWP"); }
        else if (choice == "10") { metropolisRepeat(true, "Periodic", "DWP"); }
        else if (choice == "11") { metropolisRepeat(false, "Dirichlet", "DWP"); }
        else if (choice == "12") { metropolisRepeat(true, "Dirichlet", "DWP"); }
        else if (choice == "13") { // Run multiple systems/BCs in one go, used to produce data for the report
            metropolisRepeat(false, "Periodic", "QHO");
            //metropolisRepeat(false, "Dirichlet", "QHO");
            metropolisRepeat(false, "Periodic", "AHO");
            //metropolisRepeat(false, "Dirichlet", "AHO");
            metropolisRepeat(false, "Periodic", "DWP");
            //metropolisRepeat(false, "Dirichlet", "DWP");
            std::cout << "Exiting..." << std::endl; break;
		}
        else if (choice == "0") { std::cout << "Exiting..." << std::endl; break; }
        else { std::cerr << "Invalid choice." << std::endl; }
    }
    return 0;
}

void chooseSystem() {  // Function to display user choices
    std::string chooseSystemString = "1: Perform Metropolis algorithm with periodic boundary conditions on the QHO system (without path visualisation) \n"
        "2: Perform Metropolis algorithm with periodic boundary conditions on the QHO system (with path visualisation)\n"
        "3: Perform Metropolis algorithm with Dirichlet boundary conditions on the QHO system (without path visualisation)\n"
        "4: Perform Metropolis algorithm with Dirichlet boundary conditions on the QHO system (with path visualisation)\n"
        "5: Perform Metropolis algorithm with periodic boundary conditions on the AHO system (without path visualisation)\n"
        "6: Perform Metropolis algorithm with periodic boundary conditions on the AHO system (with path visualisation)\n"
        "7: Perform Metropolis algorithm with Dirichlet boundary conditions on the AHO system (without path visualisation)\n"
        "8: Perform Metropolis algorithm with Dirichlet boundary conditions on the AHO system (with path visualisation)\n"
        "9: Perform Metropolis algorithm with Periodic boundary conditions on the DWP system (without path visualisation)\n"
        "10: Perform Metropolis algorithm with Periodic boundary conditions on the DWP system (with path visualisation)\n"
        "11: Perform Metropolis algorithm with Dirichlet boundary conditions on the DWP system (without path visualisation)\n"
        "12: Perform Metropolis algorithm with Dirichlet boundary conditions on the DWP system (with path visualisation)\n"
        "13: Run all systems with periodic boundary conditions in one go without path visualisation (useful for producing data)\n"
        "0: Exit\n"
        "Warning: Program has a tendancy to crash when ran for long times with visualisation\n"
        "Running the algorithm with path visualisation will very slightly increase the runtime of the program";

    std::cout << chooseSystemString << std::endl;
}

void metropolisRepeat(bool winOn, std::string boundary, std::string system) { // Loop over repeats
	// Set up the potentials and boundary conditions for the simulation based on user choice
    setBoundary(boundary);
	potential = findPotential(system);
	potentialDifferential = findPotentialDifferential(system);
    
    // Window thread setup
    winRunning = winOn;
    std::thread windowThread(window, std::ref(positions), std::ref(winRunning));

    // Loop the metropolis function "repeats" times 
    for (int repeat = 0; repeat < repeats; repeat++) {
        metropolis(winOn, boundary, system, repeat);
    }

    // Kill the window thread
    winRunning = false;
    windowThread.join();

    /*if (multThreads == true) {
        omp_set_num_threads(threads);
        #pragma omp parallel for
        for (int thread = 0; thread < threads; thread++) {
            for (int repeat = 0; repeat < repeats; repeat++) {
                metropolis(winOn, boundary, system, thread + repeat * threads); // Not working yet
            };
        }
    }*/
    
	// Write data to files
    //csvWriteData(boundary, system);       // Legacy csv writing functions, replaced by h5 files
    writeData(boundary, system);            // Writes all data to a single h5 file, separated into groups

	// Clear all vectors for the next run
	E0Therm.clear();
	accRateTherm.clear();
	E0Decorr.clear();
	accRateDecorr.clear();
    xTherm.clear();
	xDecorr.clear();
    GTwoDecorr.clear();
    GFourDecorr.clear();
	thermSweeps.clear();

    // Reprint the options to the user, wait for further input
    chooseSystem();
}

void metropolis(bool winOn, std::string boundary, std::string system, int repeat) { // Metropolis function which gets called initially, then calls other functions to perform the algorithm
    // Set initial path and counters to 0
    initialise(boundary, system);

    // Thermalise the system
    int thermalisationSweeps = thermalise(winOn, potentialDifferential, potential);
    if (takeMeasuresFlag == true) { takeMeasures(positions, potentialDifferential, potential); } // First measurement after thermalisation
    std::cout << "Iteration " << repeat + 1 << " thermalised after " << thermalisationSweeps << " sweeps.";
    thermSweeps.push_back(thermalisationSweeps);

	// Take measures of the path every "decorrelation" sweeps
    if (takeMeasuresFlag == true) {
        while (measureCount < measures) {
            metropolisUpdate(winOn, potential);
            sweep++;
            if (remainder(sweep, decorrelation) == 0) {
                takeMeasures(positions, potentialDifferential, potential);
            }
        }
        std::cout << " Completed measurements for iteration " << repeat + 1 << "." << std::endl;
    }
}

void metropolisUpdate(bool winOn, double (*potential)(double)) {    // The heart of the simulation, the metropolis algorithm function
    double newPosition;
    for (int i = start; i < end; i++) {
        double y = rfRange(-1, 1);  // See random.h, sets a random number between -1 and 1
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
            if (r < exp(-actionDelta)) { // Accept move with probability exp(-deltaS)
                positions[i] = newPosition; // Updates position stored (no need to store rejected position)
                acceptedMoves++;
            }
        }
    }
    // Update the window if user wanted visualisation
    if (winOn == true) {
        std::this_thread::sleep_for(std::chrono::milliseconds(delay));
    }
}

void initialise(std::string boundary, std::string system) {   // Initialises variables and path
    // Reset the path
    positions = std::vector<double>(N, 0.0); // Reset the path
    if (system == "DWP") // DWP requires a different initial (cold) path since the potential wells are not centered around 0
    {
        positions = std::vector<double>(N, wellCentres);
    }

    if (hot_start == true) {
        for (int i = 0; i < N; i++) {
            double y = rfRange(-1, 1);
            positions[i] = y * max_distance;
        }
    }

    // Reset the counters for sweeps and measures
    sweep = 0;
    measureCount = 0;

    E0ThermTemp.clear();
}

const int thermalise(bool winOn, double (*potentialDifferential)(double), double (*potential)(double)) { // Thermalisation function to reach equilibrium before measurements
	bool thermalised = false;
    while (thermalised == false) {
        metropolisUpdate(winOn, potential);
        sweep++;
        if (takeMeasuresFlag == true) {
            if ((sweep - 1) % thermalisationInterval == 0) { // Store E0 during thermalisation for plotting purposes (not necessary for the algorithm itself)
                takeThermMeasures(positions, potentialDifferential, potential);
                thermalised = checkThermalised();
            }
        }
        if (sweep >= thermalisationMaximum) {
            std::cout << "Failed to thermalise after " << thermalisationMaximum << " sweeps, proceeding with measurements anyway." << std::endl;
            return thermalisationMaximum;
            thermalised = true;
        }
    }
    return sweep;
}

void takeThermMeasures(std::vector<double>& positions, double (*potentialDifferential)(double), double (*potential)(double)) {
    // Record acceptance rate between decorrelations
    accRateTherm.push_back((double)(acceptedMoves) / (N * (double)thermalisationInterval));
    acceptedMoves = 0;

    // Record ground state energy
    E0Therm.push_back(E0Calc(positions, potentialDifferential, potential));
    E0ThermTemp.push_back(E0Calc(positions, potentialDifferential, potential));

    // Record all the positions of the particle (technically all the information we need, other vectors are for convenience)
    xTherm.insert(xTherm.end(), positions.begin(), positions.end());
}

bool checkThermalised() {
    std::vector<double> E0Samples = E0ThermTemp;

    double E0Mean = vectorMean(E0Samples);
    double E0MCSE = MCSE(E0Samples);

    //std::cout << "E0 = " << E0Mean << " +- " << E0MCSE << " (MCSE)" << std::endl;
    
    if (E0MCSE / E0Mean < acceptableError) {
        if (sweep > thermalisationMinimum) {
            return true;
        }
    }
    return false;
}

void takeMeasures(std::vector<double>& positions, double (*potentialDifferential)(double), double (*potential)(double)) {    // Takes all measurements at current path state in one function
    // Record acceptance rate between decorrelations
    accRateDecorr.push_back((double)(acceptedMoves) / (N * (double)decorrelation));
    acceptedMoves = 0;

    // Record ground state energy
    E0Decorr.push_back(E0Calc(positions, potentialDifferential, potential));

    // Record all the positions of the particle (technically all the information we need, other vectors are for convenience)
    xDecorr.insert(xDecorr.end(), positions.begin(), positions.end());

    // Compute the correlators once
    std::vector<double> tempCorrTwo = twoPointCorrelator(positions);
    std::vector<double> tempCorrFour = fourPointCorrelator(positions);

    // Write the correlators
    GTwoDecorr.insert(GTwoDecorr.end(), tempCorrTwo.begin(), tempCorrTwo.end());
    GFourDecorr.insert(GFourDecorr.end(), tempCorrFour.begin(), tempCorrFour.end());
    
    // Increment measure count
    measureCount++;
}