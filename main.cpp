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
        if (choice == "1") { multThreads = true;  metropolisRepeat(false, "Periodic", "QHO"); }
        else if (choice == "2") { metropolisRepeat(true, "Periodic", "QHO"); }
        else if (choice == "3") { multThreads = true;  metropolisRepeat(false, "Dirichlet", "QHO"); }
        else if (choice == "4") { metropolisRepeat(true, "Dirichlet", "QHO"); }
        else if (choice == "5") { multThreads = true;  metropolisRepeat(false, "Periodic", "AHO"); }
        else if (choice == "6") { metropolisRepeat(true, "Periodic", "AHO"); }
        else if (choice == "7") { multThreads = true;  metropolisRepeat(false, "Dirichlet", "AHO"); }
        else if (choice == "8") { metropolisRepeat(true, "Dirichlet", "AHO"); }
        else if (choice == "9") { multThreads = true;  metropolisRepeat(false, "Periodic", "DWP"); }
        else if (choice == "10") { metropolisRepeat(true, "Periodic", "DWP"); }
        else if (choice == "11") { multThreads = true;  metropolisRepeat(false, "Dirichlet", "DWP"); }
        else if (choice == "12") { metropolisRepeat(true, "Dirichlet", "DWP"); }
        else if (choice == "13") { // Run multiple systems/BCs in one go, used to produce data for the report
            multThreads = true; 
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
    std::string chooseSystemString = "1: Perform Metropolis algorithm with periodic boundary conditions on the QHO system (with multi-threading) \n"
        "2: Perform Metropolis algorithm with periodic boundary conditions on the QHO system (with path visualisation)\n"
        "3: Perform Metropolis algorithm with Dirichlet boundary conditions on the QHO system (with multi-threading)\n"
        "4: Perform Metropolis algorithm with Dirichlet boundary conditions on the QHO system (with path visualisation)\n"
        "5: Perform Metropolis algorithm with periodic boundary conditions on the AHO system (with multi-threading)\n"
        "6: Perform Metropolis algorithm with periodic boundary conditions on the AHO system (with path visualisation)\n"
        "7: Perform Metropolis algorithm with Dirichlet boundary conditions on the AHO system (with multi-threading)\n"
        "8: Perform Metropolis algorithm with Dirichlet boundary conditions on the AHO system (with path visualisation)\n"
        "9: Perform Metropolis algorithm with Periodic boundary conditions on the DWP system (with multi-threading)\n"
        "10: Perform Metropolis algorithm with Periodic boundary conditions on the DWP system (with path visualisation)\n"
        "11: Perform Metropolis algorithm with Dirichlet boundary conditions on the DWP system (with multi-threading)\n"
        "12: Perform Metropolis algorithm with Dirichlet boundary conditions on the DWP system (with path visualisation)\n"
        "13: Run all systems with periodic boundary conditions in one go without path visualisation (useful for producing data)\n"
        "0: Exit\n"
        "Warning: Program has a tendancy to crash when ran for long times with visualisation";

    std::cout << chooseSystemString << std::endl;
}

void metropolisRepeat(bool winOn, std::string boundary, std::string system) { // Loop over repeats
	// Set up the potentials and boundary conditions for the simulation based on user choice
    setBoundary(boundary);
	potential = findPotential(system);
	potentialDifferential = findPotentialDifferential(system);

    // Set histogram range
    if (system == "FP") {
        xMax = 4.0; // Free particle isn't really a bound system, so this range exists only to allow production of a histogram
        xMin = -4.0;
    }
	else if (system == "QHO" || system == "AHO") { // Anharmonic oscillator has the same quadratic term as the harmonic oscillator, so we can use that to set the histogram range
        double sigmaQHO = 1.0 / (std::sqrt(2.0 * m * omega)); 
        xMax = ceil(sigmaQHO * 4.0);  // Set the maximum x value for the histogram to be 4 standard deviations of the analytic ground state wavefunction
        xMin = -xMax;
    }
    else if (system == "DWP") {
		double sigmaDWP = 1.0 / (std::sqrt(omegaDWP)); // Use the frequency of the wells to calculate the standard deviation of the wavefunction in each well
        xMax = ceil(wellCentres + 4.0 * sigmaDWP); 
		xMin = -xMax; 
	}	
    // Set the bin width for the histogram based on the range of positions and the number of bins
    binWidth = (xMax - xMin) / numBins;

    if (multThreads == true) {  // Running with multiple threads, much faster.
        std::vector<RepeatData> repeatResults(repeats); // Store results of all repeats

        #pragma omp parallel for
        for (int r = 0; r < repeats; ++r) {
            std::mt19937 rng(seed + r);

            RepeatData data(N, numBins);

            metropolis(winOn, boundary, system, r, rng, data, potential, potentialDifferential);

            repeatResults[r] = data;

            std::cout << "Finished collecting data for iteration " << r + 1 << std::endl;
        }
        // Merge thread-safe results after parallel region
        E0Therm.clear();
        accRateTherm.clear(); 
        E0.clear();
        GTwo.clear();
        GFour.clear();
        positionsTemp.clear();
		histogram.clear();
        accRate.clear();

        for (int r = 0; r < repeats; ++r) {
            const auto& data = repeatResults[r];

			// Ground-state energy during thermalisation
            E0Therm.insert(E0Therm.end(), data.E0Therm.begin(), data.E0Therm.end());

			// Acceptance rates during thermalisation
            accRateTherm.insert(accRateTherm.end(), data.accRateTherm.begin(), data.accRateTherm.end());

            // Ground-state energy
            E0.insert(E0.end(), data.E0Temp.begin(), data.E0Temp.end());

            // Positions for decorrelated measurements
            positionsTemp.insert(positionsTemp.end(), data.positionsTemp.begin(), data.positionsTemp.end());

            // Acceptance rates
            accRate.insert(accRate.end(), data.accRate.begin(), data.accRate.end());

            // Correlators — sum them up for later averaging
            if (GTwo.empty()) {
                GTwo = data.GTwoTemp;
            }
            else {
                for (size_t i = 0; i < GTwo.size(); ++i)
                    GTwo[i] += data.GTwoTemp[i];
            }

            if (GFour.empty()) {
                GFour = data.GFourTemp;
            }
            else {
                for (size_t i = 0; i < GFour.size(); ++i)
                    GFour[i] += data.GFourTemp[i];
            }
        }

        // Average correlators over repeats
        for (size_t i = 0; i < GTwo.size(); ++i)
            GTwo[i] /= repeats;

        for (size_t i = 0; i < GFour.size(); ++i)
            GFour[i] /= repeats;
    }
    else { // Single-threaded version. Slower but allows visualisation.
        RepeatData data(N, numBins);
        
        // Window thread setup
        winRunning = winOn;
        std::thread windowThread(window, std::ref(data.positions), std::ref(winRunning)); 
        for (int repeat = 0; repeat < repeats; repeat++) {
            std::mt19937 rng(seed + repeat);

            metropolis(winOn, boundary, system, repeat, rng, data, potential, potentialDifferential);

            E0.insert(E0.end(), data.E0Temp.begin(), data.E0Temp.end());
            GTwo.insert(GTwo.end(), data.GTwoTemp.begin(), data.GTwoTemp.end());
            GFour.insert(GFour.end(), data.GFourTemp.begin(), data.GFourTemp.end());
            positionsTemp.insert(positionsTemp.end(), data.positionsTemp.begin(), data.positionsTemp.end());
            accRate.insert(accRate.end(), data.accRate.begin(), data.accRate.end());
            std::cout << "Finished collecting data for iteration " << repeat + 1 << std::endl;
        }
        winRunning = false;
        windowThread.join();
    }
    
	// Write data to files
    //csvWriteData(boundary, system);       // Legacy csv writing functions, replaced by h5 files
    writeData(boundary, system);            // Writes all data to a single h5 file, separated into groups

	// Clear all vectors for the next run
    E0Therm.clear();
    accRateTherm.clear();
    E0.clear();
    accRate.clear();
    positionsTemp.clear();
    GTwo.clear();
    GFour.clear();
	histogramTemp.clear();
    thermSweeps.clear();

    // Reprint the options to the user, wait for further input
    chooseSystem();
}

void metropolis(bool winOn, std::string boundary, std::string system, int repeat, std::mt19937& rng, RepeatData& data,
    double (*potential)(double), double (*potentialDifferential)(double)) { // Metropolis function which gets called initially, then calls other functions to perform the algorithm
    // Set initial path and counters to 0
    initialise(boundary, system, rng, data);

    // Thermalise the system
    int thermalisationSweeps = thermalise(winOn, potentialDifferential, potential, rng, data);

    if (takeMeasuresFlag == true) { takeMeasures(data.positions, potentialDifferential, potential, data); } // First measurement after thermalisation
    //std::cout << "Iteration " << repeat + 1 << " thermalised after " << thermalisationSweeps << " sweeps.";
    thermSweeps.push_back(thermalisationSweeps);

	// Take measures of the path every "decorrelation" sweeps
    if (takeMeasuresFlag == true) {
        while (data.measureCount < measures) {
            metropolisUpdate(winOn, potential, rng, data);
            data.sweep++;
            if (remainder(data.sweep, decorrelation) == 0) {
                takeMeasures(data.positions, potentialDifferential, potential, data);
            }
        }
        for (int n = 0; n < N; ++n) {
            data.GTwoTemp[n] /= measures;
        }
        data.vacuumPiece /= measures;  // Average of the vacuum piece
        for (int n = 0; n < N; ++n) {
            data.GFourTemp[n] = data.GFourTemp[n] / measures - data.vacuumPiece * data.vacuumPiece;
        }
        GTwo.insert(GTwo.end(), data.GTwoTemp.begin(), data.GTwoTemp.end());
        GFour.insert(GFour.end(), data.GFourTemp.begin(), data.GFourTemp.end());

		for (int i = 0; i < numBins; i++) { // Average the histogram over measures and normalise it
            data.histogramTemp[i] /= (measures * N);
        }
		histogram.insert(histogram.end(), data.histogramTemp.begin(), data.histogramTemp.end());

        //std::cout << " Completed measurements for iteration " << repeat + 1 << "." << std::endl;
    }
}

void metropolisUpdate(bool winOn, double (*potential)(double), std::mt19937& rng, RepeatData& data) {    // The heart of the simulation, the metropolis algorithm function
    double newPosition;
    for (int i = start; i < end; i++) {
        double y = uniformMinus1to1(rng);  // Sets a random number between -1 and 1
        newPosition = data.positions[i] + epsilon * y; // Incrementing one position by float between -epsilon and +epsilon
		double oldPosition = data.positions[i];
		double leftPosition = data.positions[(i - 1 + N) % N];   // Previous position, with periodic BCs
		double rightPosition = data.positions[(i + 1) % N];      // Next position, with periodic BCs
        double kineticDelta =
            (rightPosition - newPosition) * (rightPosition - newPosition) 
            - (rightPosition - oldPosition) * (rightPosition - oldPosition)
            + (newPosition - leftPosition) * (newPosition - leftPosition)
            - (oldPosition - leftPosition) * (oldPosition - leftPosition);

        double potentialDelta = potential(newPosition) - potential(oldPosition);

        double actionDelta = (m * 0.5 * aInverse) * kineticDelta + a * potentialDelta;
        // Metropolis acceptance
        if (actionDelta < 0 || uniform01(rng) < exp(-actionDelta)) { // Accepted move
            data.positions[i] = newPosition;
            data.acceptedMoves++;
        }
    }
    // Update the window if user wanted visualisation
    if (winOn == true) {
        std::this_thread::sleep_for(std::chrono::milliseconds(delay));
    }
}

void initialise(std::string boundary, std::string system, std::mt19937& rng, RepeatData& data) {   // Initialises variables and path
    // Reset the path
    data.positions = std::vector<double>(N, 0.0); // Reset the path
    if (system == "DWP") // DWP requires a different initial (cold) path since the potential wells are not centered around 0
    {
        data.positions = std::vector<double>(N, wellCentres);
    }
    
    data.GTwoTemp = std::vector<double>(N, 0.0);
    data.GFourTemp = std::vector<double>(N, 0.0);
    data.E0Temp.clear();

    data.sweep = 0;
    data.measureCount = 0;
    data.acceptedMoves = 0;
    data.vacuumPiece = 0.0;

    if (hot_start) {
        for (int i = 0; i < N; i++)
            data.positions[i] = uniformMinus1to1(rng) * max_distance;
    }
}

const int thermalise(bool winOn, double (*potentialDifferential)(double), double (*potential)(double), std::mt19937& rng, RepeatData& data) { // Thermalisation function to reach equilibrium before measurements
	bool thermalised = false;
    if (takeThermMeasuresFlag == true) {
        while (thermalised == false) {
            metropolisUpdate(winOn, potential, rng, data);
            data.sweep++;

            if ((data.sweep - 1) % thermalisationInterval == 0) { // Store E0 during thermalisation
                takeThermMeasures(data.positions, potentialDifferential, potential, data);
                thermalised = checkThermalised(data);   // Only check thermalisation every "thermalisationInterval" sweeps to reduce calls of the expensive MCSE function
            }
            if (data.sweep >= thermalisationMaximum) {
                std::cout << "Failed to thermalise after " << thermalisationMaximum << " sweeps, proceeding with measurements anyway." << std::endl;
                return thermalisationMaximum;
            }
        }
    }
    else { // If we are not taking measures during thermalisation, just run the sweeps without checking for thermalisation
        while (data.sweep < thermalisationMaximum) {
            metropolisUpdate(winOn, potential, rng, data);
            data.sweep++;
        }
	}
    return data.sweep;
}

void takeThermMeasures(std::vector<double>& positions, double (*potentialDifferential)(double), double (*potential)(double), RepeatData& data) {
    // Record acceptance rate between decorrelations
    data.accRateTherm.push_back((double)(data.acceptedMoves) / (N * (double)thermalisationInterval));
    data.acceptedMoves = 0;

    // Record ground state energy
    data.E0Therm.push_back(E0Calc(positions, potentialDifferential, potential));
    data.E0ThermTemp.push_back(E0Calc(positions, potentialDifferential, potential));
}

bool checkThermalised(const RepeatData& data) {
    double E0Mean = vectorMean(data.E0ThermTemp);
    double E0MCSE = MCSE(data.E0ThermTemp);

    if (E0MCSE / E0Mean < acceptableError) {
        if (data.sweep > thermalisationMinimum) {
            return true;
        }
    }
    return false;
}

void takeMeasures(std::vector<double>& positions, double (*potentialDifferential)(double), double (*potential)(double), RepeatData& data) {    // Takes all measurements at current path state in one function
    // Record acceptance rate between decorrelations
    data.accRate.push_back((double)(data.acceptedMoves) / (N * (double)decorrelation));
    data.acceptedMoves = 0;

    // Record ground state energy
    double E0 = E0Calc(positions, potentialDifferential, potential);
    data.E0Temp.push_back(E0);

    // Record all the positions of the particle
    data.positionsTemp.insert(data.positionsTemp.end(), positions.begin(), positions.end());

    // Compute the correlators once
    std::vector<double> tempCorrTwo = twoPointCorrelator(positions);
    std::vector<double> tempCorrFour = fourPointCorrelator(positions, data.vacuumPiece);

    // Write the correlators
    if (data.measureCount == 0) {
        data.GTwoTemp = tempCorrTwo;
    }
    else {
        for (int n = 0; n < N; ++n)
            data.GTwoTemp[n] += tempCorrTwo[n];
    }

    if (data.measureCount == 0) {
        data.GFourTemp = tempCorrFour;
    }
    else {
        for (int n = 0; n < N; ++n)
            data.GFourTemp[n] += tempCorrFour[n];
    }

    for (int t = 0; t < N; t++) {
        double x = positions[t];

        int bin = int((x - xMin) / binWidth);

        if (bin >= 0 && bin < numBins) {
            data.histogramTemp[bin] += 1.0;
        }
    }

    // Increment measure count
    data.measureCount++;
}