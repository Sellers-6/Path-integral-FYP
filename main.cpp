///// Simulating quantum systems with different potentials using Feynman path integral and Metropolis algorithm /////
// Throughout the program, hbar is assumed to be 1 (reduced planck's constant) //

#include "main.h"

int main() {    // Accepts user choice of system type and window visualisation
    std::cout << "Simulating quantum systems with different potentials using the Feynman path integral and Metropolis algorithm";
    std::string choice;
	chooseSystem(); // Print the options to the user
	while (true) {
		std::cout << "Enter choice: ";
        std::cin >> choice;
        if (choice == "1")      { multThreads = true;  metropolisRepeat(false, "QHO"); }
        else if (choice == "2") { metropolisRepeat(true, "QHO"); }
        else if (choice == "3") { multThreads = true;  metropolisRepeat(false, "DWP"); }
        else if (choice == "4") { metropolisRepeat(true, "DWP"); }
        else if (choice == "5") { // Run both systems in one go, used to produce data for the report
            multThreads = true; 
            metropolisRepeat(false, "QHO");
            metropolisRepeat(false, "DWP");
            std::cout << "Exiting..." << std::endl; break;
		}
        else if (choice == "6") {  // Run DWP system with multi-threading with varying lattice spacing and well separation.
            // This is used to find the error in the energy as a function of the lattice spacing, and to find how the tunnelling rate changes with well separation.
            multThreads = true;
            sixFlag = true;
            std::vector<double> wellCentresVec = { 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0 }; // Varying the separation of the wells.
            std::vector<double> latticeSpacingsVec = { 0.1 }; // Varying the lattice spacing to see how it affects the error in the energy.
            std::vector<double> esplinonVec = { 0.4 };
			int epsilonIndex = 0; // Index to keep track of which epsilon value we are using.
            for (double latticeSpacing : latticeSpacingsVec) {
                a = latticeSpacing; // Update the lattice spacing.
				epsilon = esplinonVec[epsilonIndex]; // Update epsilon based on the lattice spacing, to keep the acceptance rate at ~ 60%. This was found through trial and error, and is not an exact relationship.
                for (double wellCentresElement : wellCentresVec) {
					std::cout << "\n Running DWP system with lattice spacing a = " << latticeSpacing << " and well centres = " << wellCentresElement << std::endl;
                    wellCentres = wellCentresElement; // Update the well centres based on the separation.
                    setParameters(10000, a, wellCentres); // Use beta = 10000 so that tunneling is captured in the all regimes of well separation.
                    metropolisRepeat(false, "DWP");
                }
				epsilonIndex++;
            }
        }
        else if (choice == "0") { std::cout << "Exiting..." << std::endl; break; }
        else { std::cerr << "Invalid choice." << std::endl; }

        // Reprint the options to the user, wait for further input
        chooseSystem();
    }
    return 0;
}

void chooseSystem() {  // Function to display user choices
    std::string chooseSystemString = "\n1: Perform Metropolis algorithm on the QHO system (with multi-threading) \n"
        "2: Perform Metropolis algorithm on the QHO system (with path visualisation)\n"
        "3: Perform Metropolis algorithm on the DWP system (with multi-threading)\n"
        "4: Perform Metropolis algorithm on the DWP system (with path visualisation)\n"
        "5: Run both systems with multi-threading, run this to reproduce results from the report\n"
		"6: Run DWP system with multi-threading with varying lattice spacing and well separation,\n"
        "   this is used to produce data for the error analysis and tunnelling rate plots in the report\n"
        "0: Exit\n"
        "Warning: Program has a tendancy to crash when ran for long times with visualisation";

    std::cout << chooseSystemString << std::endl;
}

void metropolisRepeat(bool winOn, std::string system) { // Loop over repeats
    std::cout << "Performing MCMC algorithm for the " << system << "." << std::endl;

	// Set up the potentials for the simulation based on user choice
	potential = findPotential(system);
	potentialDifferential = findPotentialDifferential(system);

	// Set histogram range, deccorelation sweeps, and thermalisation sweeps based on the system being simulated
	if (system == "QHO") { // Anharmonic oscillator has the same quadratic term as the harmonic oscillator, so we can use that to set the histogram range
        double sigmaQHO = 1.0 / (std::sqrt(2.0 * m * omega)); 
        xMax = ceil(sigmaQHO * 4.0) + 1;  // Set the maximum x value for the histogram to be 4 standard deviations of the analytic ground state wavefunction (plus padding)
        xMin = -xMax;
        decorrSweeps = decorrSweepsQHO;
		thermSweeps = thermSweepsQHO;
    }
    else if (system == "DWP") {
		double sigmaDWP = 1.0 / (std::sqrt(omegaDWP)); // Use the frequency of the wells to calculate the standard deviation of the wavefunction in each well
        xMax = ceil(wellCentres + 4.0 * sigmaDWP) + 1; 
		xMin = -xMax;
        decorrSweeps = decorrSweepsDWP;
		thermSweeps = thermSweepsDWP;
	}	
    // Set the bin width for the histogram based on the range of positions and the number of bins
    binWidth = (xMax - xMin) / numBins;

    if (multThreads == true) {  // Running with multiple threads, much faster.
        std::cout << "Using multiple threads to speed up the simulation. Taking " << measures << " measures and running " << repeats << " repeats." << std::endl;
        std::cout << "                                           Percent complete                                           " << std::endl;
        std::cout << "<       10        20        30        40        50        60        70        80        90        100>" << std::endl;
        std::cout << "<";
        std::vector<RepeatData> repeatResults(repeats); // Store results of all repeats
        double iterationNumber = 1.0;
        double percentDone = 0.0;
        #pragma omp parallel for
        for (int r = 0; r < repeats; ++r) {
            std::mt19937 rng(seed + r);

            RepeatData data(N, numBins, measures);

            metropolis(winOn, system, r, rng, data, potential, potentialDifferential);

            repeatResults[r] = data;

            while ((iterationNumber / repeats) * 100.0 > percentDone) {
                std::cout << "-";
                percentDone++;
            }
            
            iterationNumber++;
        }
        std::cout << ">";

        // Merge thread-safe results after parallel region
        E0Therm.clear();
        accRateTherm.clear(); 
        E0.clear();
        Gx1x1.clear();
        Gx1x2.clear();
        Gx1x3.clear();
        Gx2x2.clear();
        Gx2x3.clear();
        Gx3x3.clear();
        positions.clear();
		histogram.clear();
		instantons.clear();
		antiInstantons.clear();
        accRate.clear();

        for (int r = 0; r < repeats; ++r) {
            const auto& data = repeatResults[r];

			// Acceptance rate during thermalisation
			accRateTherm.insert(accRateTherm.end(), data.accRateThermTemp.begin(), data.accRateThermTemp.end());

			// Acceptance rate during measurements
			accRate.insert(accRate.end(), data.accRateTemp.begin(), data.accRateTemp.end());

			// Ground-state energy during thermalisation
            E0Therm.insert(E0Therm.end(), data.E0ThermTemp.begin(), data.E0ThermTemp.end());

            // Ground-state energy
            E0.insert(E0.end(), data.E0Temp.begin(), data.E0Temp.end());


            // Correlators — sum them up for later averaging
            if (Gx1x1.empty()) {
                Gx1x1 = data.Gx1x1Temp;
            }
            else {
                for (size_t i = 0; i < Gx1x1.size(); ++i)
                    Gx1x1[i] += data.Gx1x1Temp[i];
            }

            if (Gx1x2.empty()) {
                Gx1x2 = data.Gx1x2Temp;
            }
            else {
                for (size_t i = 0; i < Gx1x2.size(); ++i)
                    Gx1x2[i] += data.Gx1x2Temp[i];
            }

            if (Gx1x3.empty()) {
                Gx1x3 = data.Gx1x3Temp;
            }
            else {
                for (size_t i = 0; i < Gx1x3.size(); ++i)
                    Gx1x3[i] += data.Gx1x3Temp[i];
            }

            if (Gx2x2.empty()) {
                Gx2x2 = data.Gx2x2Temp;
            }
            else {
                for (size_t i = 0; i < Gx2x2.size(); ++i)
                    Gx2x2[i] += data.Gx2x2Temp[i];
            }

            if (Gx2x3.empty()) {
                Gx2x3 = data.Gx2x3Temp;
            }
            else {
                for (size_t i = 0; i < Gx2x3.size(); ++i)
                    Gx2x3[i] += data.Gx2x3Temp[i];
            }

            if (Gx3x3.empty()) {
                Gx3x3 = data.Gx3x3Temp;
            }
            else {
                for (size_t i = 0; i < Gx3x3.size(); ++i)
                    Gx3x3[i] += data.Gx3x3Temp[i];
            }

            // Histogram for the wavefunction
            if (histogram.empty()) {
                histogram = data.histogramTemp;  
            }
            else {
                for (int i = 0; i < numBins; ++i)
                    histogram[i] += data.histogramTemp[i];
            }

            instantons.insert(instantons.end(), data.instantonsTemp.begin(), data.instantonsTemp.end());
            antiInstantons.insert(antiInstantons.end(), data.antiInstantonsTemp.begin(), data.antiInstantonsTemp.end());
        }

        // Average correlators over repeats
        for (size_t i = 0; i < Gx1x1.size(); ++i)
            Gx1x1[i] /= repeats;

        for (size_t i = 0; i < Gx2x2.size(); ++i)
            Gx2x2[i] /= repeats;
    }
    else { // Visualisation.
        std::cout << "Providing a visualisation of the evolution of the path. Running " << repeats << " repeats and taking " << measures << " measures." << std::endl;

        RepeatData data(N, numBins);    // Still use this for the single threaded version, though it's not necessary it is convenient
        
        // Window thread setup
        std::thread windowThread(window, std::ref(data.positions), std::ref(winOn)); 
        std::mt19937 rng(seed);

        while (true) {
            metropolisUpdate(winOn, potential, rng, data);
        }
    }

    constructHeaderInfo(system);
    
	// Write data to files
    //csvWriteData(system);       // Legacy csv writing functions, replaced by h5 files
    if (!sixFlag) { writeData(system); }            // Writes all data to a single h5 file, separated into groups
	else if (sixFlag) { writeData2(system, std::to_string(wellCentres), std::to_string(a)); } // Writes data for option 10
	
    // Clear all vectors for the next run
    E0Therm.clear();
    accRateTherm.clear();
    E0.clear();
    accRate.clear();
    positions.clear();
    Gx1x1.clear();
    Gx1x2.clear();
    Gx1x3.clear();
    Gx2x2.clear();
    Gx2x3.clear();
    Gx3x3.clear();
	histogram.clear();
	instantons.clear();
	antiInstantons.clear();
	headerInfo.clear();
}

void metropolis(bool winOn, std::string system, int repeat, std::mt19937& rng, RepeatData& data,
    double (*potential)(double), double (*potentialDifferential)(double)) { // Metropolis function which gets called initially, then calls other functions to perform the algorithm
    // Set initial path and counters to 0
    initialise(system, rng, data);

    // Thermalise the system
    thermalise(winOn, potentialDifferential, potential, rng, data);

    if (takeMeasuresFlag == true) { takeMeasures(data.positions, data); } // First measurement after thermalisation

	// Take measures of the path every "decorrelation" sweeps
    if (takeMeasuresFlag == true) {
        while (data.measureCount < measures) {
            metropolisUpdate(winOn, potential, rng, data);
            data.sweep++;
            if (remainder(data.sweep, decorrSweeps) == 0) {
                takeMeasures(data.positions, data);
            }
        }
		computeObservables(data.positionsVector, potentialDifferential, potential, data); // Compute the observables based on the measured paths
    }
}

void metropolisUpdate(bool winOn, double (*potential)(double), std::mt19937& rng, RepeatData& data) {    // The core of the simulation, the metropolis algorithm function
    double newPosition;
    for (int i = 0; i < N; i++) {
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

void initialise(std::string system, std::mt19937& rng, RepeatData& data) {   // Initialises variables and path
    // Reset the path
    data.positions = std::vector<double>(N, 0.0); // Reset the path
    if (system == "DWP") // DWP requires a different initial (cold) path since the potential wells are not centered around 0
    {
        data.positions = std::vector<double>(N, wellCentres * side);
		side *= -1; // Start the next run in the opposite well to encourage tunnelling and faster thermalisation
    }
    
    data.Gx1x1Temp = std::vector<double>(N, 0.0);
    data.Gx2x2Temp = std::vector<double>(N, 0.0);
    data.E0Temp.clear();

    data.sweep = 0;
    data.measureCount = 0;
    data.acceptedMoves = 0;

    if (hot_start) {
        for (int i = 0; i < N; i++)
            data.positions[i] = uniformMinus1to1(rng) * max_distance;
    }

    if (split_wells) {
        for (int i = 0; i < N; i++) {
            if (i < N / 2) {
                data.positions[i] = wellCentres; // Start half the path in one well
            }
            else {
                data.positions[i] = -wellCentres; // Start the other half in the other well
            }
        }
    }
}

void thermalise(bool winOn, double (*potentialDifferential)(double), double (*potential)(double), std::mt19937& rng, RepeatData& data) { // Thermalisation function to reach equilibrium before measurements
    while (data.sweep < thermSweeps) {
        metropolisUpdate(winOn, potential, rng, data);
        data.sweep++;

        if ((data.sweep - 1) % thermInterval == 0) { // Take thermalisation measures
            takeThermMeasures(data.positions, potentialDifferential, potential, data);
        }
    }
}

void takeThermMeasures(std::vector<double>& positions, double (*potentialDifferential)(double), double (*potential)(double), RepeatData& data) {
    // Record acceptance rate between decorrelations
    data.accRateThermTemp.push_back((double)(data.acceptedMoves) / (N * (double)thermInterval));
    data.acceptedMoves = 0;

    // Record ground state energy
    data.E0ThermTemp.push_back(E0Calc(positions, potentialDifferential, potential));
}

void takeMeasures(std::vector<double>& positions, RepeatData& data) {    // Takes all measurements at current path state in one function
    // Record acceptance rate between decorrelations
    data.accRateTemp.push_back((double)(data.acceptedMoves) / (N * (double)decorrSweeps));
    data.acceptedMoves = 0;

	// Record positions for observable calculations later.
    for (int i = 0; i < N; i++) {
        data.positionsVector[data.measureCount][i] = positions[i];
    }

    // Increment measure count
    data.measureCount++;
}

void computeObservables(std::vector<std::vector<double>>& positionsVector, double (*potentialDifferential)(double), double (*potential)(double), RepeatData& data) {

    std::vector<std::vector<double>> O1(measures, std::vector<double>(N));
    std::vector<std::vector<double>> O2(measures, std::vector<double>(N));
    std::vector<std::vector<double>> O3(measures, std::vector<double>(N));

    for (int cfg = 0; cfg < measures; cfg++)
    {
        for (int t = 0; t < N; t++)
        {
            double x = positionsVector[cfg][t];

            O1[cfg][t] = x;
            O2[cfg][t] = x * x;
            O3[cfg][t] = x * x * x;
        }
    }

    double x1Vacuum = 0.0;
    double x2Vacuum = 0.0;
	double x3Vacuum = 0.0;

    for (int cfg = 0; cfg < measures; cfg++)
    {
        for (int t = 0; t < N; t++)
        {
            x1Vacuum += O1[cfg][t];
			x2Vacuum += O2[cfg][t];
            x3Vacuum += O3[cfg][t];
        }
	}

	x1Vacuum /= (measures * N);
	x2Vacuum /= (measures * N);
	x3Vacuum /= (measures * N);

    for (int cfg = 0; cfg < measures; cfg++)
    {
        std::vector<double> correlation11Temp;
        std::vector<double> correlation12Temp;
        std::vector<double> correlation13Temp;
        std::vector<double> correlation22Temp;
        std::vector<double> correlation23Temp;
        std::vector<double> correlation33Temp;
        correlation11Temp = correlator(O1[cfg], O1[cfg]);
        correlation12Temp = correlator(O1[cfg], O2[cfg]);
        correlation13Temp = correlator(O1[cfg], O3[cfg]);
        correlation22Temp = correlator(O2[cfg], O2[cfg]);
        correlation23Temp = correlator(O2[cfg], O3[cfg]);
        correlation33Temp = correlator(O3[cfg], O3[cfg]);
            
        for (int i = 0; i < N; i++) {
            data.Gx1x1Temp[i] += correlation11Temp[i];
            data.Gx1x2Temp[i] += correlation12Temp[i];
            data.Gx1x3Temp[i] += correlation13Temp[i];
            data.Gx2x2Temp[i] += correlation22Temp[i];
            data.Gx2x3Temp[i] += correlation23Temp[i];
            data.Gx3x3Temp[i] += correlation33Temp[i];
        }
    }

    // Average correlators
    for (int n = 0; n < N; ++n) {
        data.Gx1x1Temp[n] /= measures;
		data.Gx1x1Temp[n] -= x1Vacuum * x1Vacuum; 
        data.Gx1x2Temp[n] /= measures;
		data.Gx1x2Temp[n] -= x1Vacuum * x2Vacuum;
        data.Gx1x3Temp[n] /= measures;
		data.Gx1x3Temp[n] -= x1Vacuum * x3Vacuum;
        data.Gx2x2Temp[n] /= measures;
		data.Gx2x2Temp[n] -= x2Vacuum * x2Vacuum;
        data.Gx2x3Temp[n] /= measures;
		data.Gx2x3Temp[n] -= x2Vacuum * x3Vacuum;
        data.Gx3x3Temp[n] /= measures;
		data.Gx3x3Temp[n] -= x3Vacuum * x3Vacuum;
    }

    for (int i = 0; i < measures; i++)
    {
        // Record ground state energy
        double E0 = E0Calc(positionsVector[i], potentialDifferential, potential);
        data.E0Temp.push_back(E0); // Save every measure of the ground state energy to calculate the error later


        // Record histogram data for the wavefunction
        for (int t = 0; t < N; t++) {
            double x = positionsVector[i][t];

            int bin = int((x - xMin) / binWidth);

            if (bin >= 0 && bin < numBins) {
                data.histogramTemp[bin] += 1.0;
            }
        }
        // Find instantons and anti instantons
        int prev = whichWell(positionsVector[i][0], tunnellingThreshold);

        int instantons = 0;
        int antiInstantons = 0;

        for (int j = 1; j < N; j++)
        {
            int curr = whichWell(positionsVector[i][j], tunnellingThreshold);

            if (curr == 0) continue;   // Ignore barrier

            if (prev == -1 && curr == 1)
                instantons++;

            if (prev == 1 && curr == -1)
                antiInstantons++;

            prev = curr;
        }

        // Check for tunnelling in the first position to account for the periodic boundary conditions
        int curr = whichWell(positionsVector[i][0], tunnellingThreshold);

        if (curr != 0)  // Ignore barrier
        {
            if (prev == -1 && curr == 1)
                instantons++;

            if (prev == 1 && curr == -1)
                antiInstantons++;
        }

        data.instantonsTemp.push_back(instantons);
        data.antiInstantonsTemp.push_back(antiInstantons);
    }

    for (int i = 0; i < numBins; i++) { // Average the histogram over measures and normalise it
        data.histogramTemp[i] /= (measures * N);
    }
}