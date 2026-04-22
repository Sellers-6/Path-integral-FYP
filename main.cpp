#include "main.h"

/***********************************************************************/
/*********** Euclidean path propagation with MCMC algorithm ************/
/***********************************************************************/

int main() {
    std::cout << "Euclidean path propagation with MCMC algorithm.";
    std::string choice;
    // Print the simulation options to the user.
    chooseSystem();
	while (true) {
        // Accept user choice of simulation to perform.
		std::cout << "Enter choice: ";
        std::cin >> choice;
        if (choice == "1")      { // Run QHO with multi-threading.
            multThreads = true;  
            std::cout << "Using multiple threads to simulate the QHO. Taking " << measures << " measures and running " << repeats << " repeats." << std::endl;
            metropolisRepeat("QHO"); 
        }
        else if (choice == "2") { // Run QHO with visualisation.
            winRunning = true;
            std::cout << "Providing a visualisation of the Euclidean path for the QHO." << std::endl;
            metropolisRepeat("QHO"); 
        }
        else if (choice == "3") { // Run DWP with multi-threading.
            multThreads = true;  
            std::cout << "Using multiple threads to simulate the DWP. Taking " << measures << " measures and running " << repeats << " repeats." << std::endl;
            metropolisRepeat("DWP"); 
        }
        else if (choice == "4") { // Run DWP with visualisation.
            winRunning = true;
            std::cout << "Providing a visualisation of the Euclidean path for the DWP." << std::endl;
            metropolisRepeat("DWP"); 
        }
        else if (choice == "5") { // Run both systems with multi-threading.
            std::cout << "Using multiple threads to simulate the QHO and DWP. Taking " << measures << " measures and running " << repeats << " repeats." << std::endl;
            multThreads = true;
            metropolisRepeat("QHO");
            metropolisRepeat("DWP");
            std::cout << "Exiting..." << std::endl; break;
		}
        else if (choice == "6") { // Run DWP system with multi-threading with varying beta and well separation.
            multThreads = true;
            sixFlag = true;
            std::vector<double> wellCentresVec =        { 1.000, 1.100, 1.200, 1.300, 1.4000, 1.50000 };
            std::vector<double> thermalisationVec =     { 10000, 10000, 10000, 10000, 100000, 1000000 };
            std::vector<double> betaVec =               { 500, 400, 250, 100, 75, 50, 25, 10 };
            int thermalisationIndex = 0;
			measures = 1; repeats = 32;
            for (double betaElement : betaVec) {
                a = 0.05; epsilon = 0.32;                           // Fix lattice spacing and epsilon.
                for (double wellCentresElement : wellCentresVec) {
                    thermSweepsDWP = thermalisationVec[thermalisationIndex % thermalisationVec.size()];
                    decorrSweepsDWP = thermSweepsDWP / 2;
                    wellCentres = wellCentresElement;               // Update the well centres based on the separation.
                    setParameters(betaElement, a, wellCentres);
                    std::cout << "\n Running DWP system with beta = " << betaElement << " and well centres = " << wellCentresElement << std::endl;
                    metropolisRepeat("DWP");
                    thermalisationIndex++;
                }
            }
        }
        else if (choice == "7") { // Run DWP system with multi-threading with varying lattice spacing and well separation.
            multThreads = true;
            sevenFlag = true;
            std::vector<double> wellCentresVec =        { 1.000, 1.100, 1.200, 1.300, 1.4000, 1.50000 };
            std::vector<double> thermalisationVec =     { 10000, 10000, 10000, 10000, 100000, 1000000 };
            std::vector<double> latticeSpacingsVec =    { 0.5, 0.40, 0.3, 0.25, 0.20, 0.175, 0.15, 0.125, 0.10, 0.075, 0.05 };
            std::vector<double> esplinonVec =           { 0.7, 0.65, 0.6, 0.57, 0.55, 0.520, 0.50, 0.480, 0.45, 0.400, 0.30 };
            int thermalisationIndex = 0;
            int epsilonIndex = 0;
            measures = 1; repeats = 32;
            for (double latticeSpacingsElement : latticeSpacingsVec) {
                a = latticeSpacingsElement;                         // Update the lattice spacing.
                epsilon = esplinonVec[epsilonIndex];                // Update epsilon based on the lattice spacing, smaller lattice spacings require smaller epsilon for good acceptance rates.
                for (double wellCentresElement : wellCentresVec) {
                    thermSweepsDWP = thermalisationVec[thermalisationIndex % thermalisationVec.size()];
                    decorrSweepsDWP = thermSweepsDWP / 2;
                    wellCentres = wellCentresElement;               // Update the well centres based on the separation.
					setParameters(500, a, wellCentres);             // Use the same beta for all runs, large enough to capture the ground state accurately.
                    std::cout << "\n Running DWP system with lattice spacing = " << latticeSpacingsElement << " and well centres = " << wellCentresElement << std::endl;
                    metropolisRepeat("DWP");
                    thermalisationIndex++;
                }
                epsilonIndex++;
            }
        }
        else if (choice == "8") { // Run DWP system with multi-threading with varying well separation.
            multThreads = true;
            eightFlag = true;
            std::vector<double> wellCentresVec =        {   0.025, 0.050, 0.075, 0.100, 0.125, 0.150, 0.175, 0.200, 0.225, 0.250, 0.275, 0.300, 0.325, 0.350, 0.375, 0.400,  0.425,  0.450,  0.475,  0.500,
                                                            0.525, 0.550, 0.575, 0.600, 0.625, 0.650, 0.675, 0.700, 0.725, 0.750, 0.775, 0.800, 0.825, 0.850, 0.875, 0.900,  0.925,  0.950,  0.975,  1.000, 
                                                            1.025, 1.050, 1.075, 1.100, 1.125, 1.150, 1.175, 1.200, 1.225, 1.250, 1.275, 1.300, 1.325, 1.350, 1.375, 1.400,  1.425,  1.450,  1.475,  1.500 };
            std::vector<double> thermalisationVec =     {   10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000,  10000,  10000,  10000,  10000, 
                                                            10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000,  10000,  10000,  10000,  10000, 
                                                            10000, 10000, 10000, 10000, 10000, 10000, 10000, 25000, 25000, 25000, 50000, 50000, 50000, 50000, 75000, 100000, 250000, 500000, 750000, 1000000 };
            int thermalisationIndex = 0;
            measures = 10; repeats = 32;
            a = 0.075; epsilon = 0.35;                              // Found to give 60% acceptance rate with lattice spacing = 0.075.
            for (double wellCentresElement : wellCentresVec) {
                thermSweepsDWP = thermalisationVec[thermalisationIndex];
                decorrSweepsDWP = thermSweepsDWP / 2;
                wellCentres = wellCentresElement;                   // Update the well centres based on the separation.
                setParameters(250, 0.075, wellCentres);             // Run for fixed lattice spacing and inverse temperature beta, found to be big enough in analysis.
                std::cout << "\n Running DWP system with well centres = " << wellCentresElement << std::endl;
                metropolisRepeat("DWP");
                thermalisationIndex++;
            }
        }
        else if (choice == "0") { // Exit the program.
            std::cout << "Exiting..." << std::endl; 
            break;
        }
		else {                    // Invalid choice, reprompt the user.
            std::cerr << "Invalid choice." << std::endl; 
        }
        // Reprint the simulation options to the user.
        chooseSystem();
    }
    return 0;
}

void chooseSystem() {
    // User choices.
    std::string chooseSystemString = "\n1: Perform Metropolis algorithm on the QHO system (with multi-threading) \n"
        "2: Perform Metropolis algorithm on the QHO system (with path visualisation)\n"
        "3: Perform Metropolis algorithm on the DWP system (with multi-threading)\n"
        "4: Perform Metropolis algorithm on the DWP system (with path visualisation)\n"
        "5: Run the QHO and DWP systems with multi-threading\n"
		"6: Run DWP system with multi-threading with varying beta and well separation,\n"
        "   this is used to produce data for the error analysis of finite beta effects on energy levels.\n"
        "7: Run DWP system with multi-threading with varying lattice spacing and well separation,\n"
        "   this is used to produce data for the error analysis of discretisation on energy levels.\n"
        "8: Run DWP system with multi-threading with varying well separation.\n"
        "0: Exit\n"
        "Warning: Program has a tendancy to crash when ran for long times with visualisation";
    std::cout << chooseSystemString << std::endl;
}

void metropolisRepeat(std::string system) { 
	// Set up the potentials for the simulation based on user choice.
	potential = findPotential(system);
	potentialDifferential = findPotentialDifferential(system);

	// Set histogram range, deccorelation sweeps, and thermalisation sweeps based on the system being simulated.
	if (system == "QHO") {
        double sigmaQHO = 1.0 / (std::sqrt(2.0 * m * omega)); 
        xMax = ceil(sigmaQHO * 4.0) + 1;                            // Set the maximum x value for the histogram to be 4 standard deviations of the analytic ground state wavefunction plus padding.
        xMin = -xMax;
        decorrSweeps = decorrSweepsQHO;
		thermSweeps = thermSweepsQHO;
    }
    else if (system == "DWP") {
		double sigmaDWP = 1.0 / (std::sqrt(omegaDWP));              // Use the frequency of the wells to calculate the standard deviation of the wavefunction in each well.
        xMax = ceil(wellCentres + 4.0 * sigmaDWP) + 1; 
		xMin = -xMax;
        decorrSweeps = decorrSweepsDWP;
		thermSweeps = thermSweepsDWP;
	}	
    // Set the bin width for the histogram based on the range of positions and the number of bins.
    binWidth = (xMax - xMin) / numBins;

    if (multThreads == true) { // Run with multiple threads, much faster.
        std::cout << "                                           Percent complete                                           " << std::endl;
        std::cout << "<       10        20        30        40        50        60        70        80        90        100>" << std::endl;
        std::cout << "<";
        std::vector<RepeatData> repeatResults(repeats);             // Store results of all repeats.
        double iterationNumber = 1.0;
        double percentDone = 0.0;
        #pragma omp parallel for
        for (int r = 0; r < repeats; ++r) {
            std::mt19937 rng(seed + r);                             // Use a different seed for each repeat to ensure different random numbers across repeats, but still have reproducibility.

			RepeatData data(N, numBins, measures);                  // Construct a RepeatData object to store the results of this repeat.

            metropolis(system, r, rng, data, potential, potentialDifferential);

            repeatResults[r] = data;

            while ((iterationNumber / repeats) * 100.0 > percentDone) {
                std::cout << "-";
                percentDone++;
            }
            
            iterationNumber++;
        }
        std::cout << ">";

        // Merge thread-safe results after parallel region. Update to include all observables.
        positions.clear();
        
        accRateTherm.clear();
        accRate.clear();
        E0Therm.clear();
        E0.clear();
        histogram.clear();
        Gx1x1.clear();
        Gx2x2.clear();
		instantons.clear();
		antiInstantons.clear();

        for (int r = 0; r < repeats; ++r) {
            const auto& data = repeatResults[r];

			accRateTherm.insert(accRateTherm.end(), data.accRateThermTemp.begin(), data.accRateThermTemp.end());

			accRate.insert(accRate.end(), data.accRateTemp.begin(), data.accRateTemp.end());

            E0Therm.insert(E0Therm.end(), data.E0ThermTemp.begin(), data.E0ThermTemp.end());

            E0.insert(E0.end(), data.E0Temp.begin(), data.E0Temp.end());

            if (histogram.empty()) { histogram = data.histogramTemp; }
            else {
                for (int i = 0; i < numBins; ++i)
                    histogram[i] += data.histogramTemp[i];
            }

            if (Gx1x1.empty()) { Gx1x1 = data.Gx1x1Temp; }
            else { 
                for (size_t i = 0; i < Gx1x1.size(); ++i) { 
                    Gx1x1[i] += data.Gx1x1Temp[i]; 
                } 
            }

            if (Gx2x2.empty()) { Gx2x2 = data.Gx2x2Temp; }
            else {
                for (size_t i = 0; i < Gx2x2.size(); ++i)
                    Gx2x2[i] += data.Gx2x2Temp[i];
            }           

            instantons.insert(instantons.end(), data.instantonsTemp.begin(), data.instantonsTemp.end());
            antiInstantons.insert(antiInstantons.end(), data.antiInstantonsTemp.begin(), data.antiInstantonsTemp.end());
        }

        // Average correlators over repeats.
        for (size_t i = 0; i < Gx1x1.size(); ++i) { Gx1x1[i] /= repeats; }
        for (size_t i = 0; i < Gx2x2.size(); ++i) { Gx2x2[i] /= repeats; }
    }
    else { // Visualisation.
        RepeatData data(N, numBins);                                // Remains convenient to construct the RepeatData struct.
        
        // Window thread setup.
        std::thread windowThread(window, std::ref(data.positions), std::ref(winRunning)); 
        std::mt19937 rng(seed);

        initialise(system, rng, data);

        int currentSweep = 0;
        int checkPointInterval = 10000;
        int checkPoint = checkPointInterval;
		while (true) { // Keep updating the path and the window until the user closes the window.
            metropolisSweep(potential, rng, data);
            currentSweep++;
            if (currentSweep == checkPoint) {
                checkPoint += checkPointInterval;
                std::cout << "Sweep number: " << currentSweep << std::endl; // This can be commented out, but is useful for telling how many sweeps tunneling takes in the DWP system.
            }
        }
    }

    constructHeaderInfo(system);
    
	// Write data to files.
    if (!sixFlag && !sevenFlag && !eightFlag) { writeData(system); }                                // Writes all data to a single h5 file, separated into groups.
	else if (sixFlag) { writeData6(system, std::to_string(wellCentres), std::to_string(N * a)); }   // Writes data for option six.
    else if (sevenFlag) { writeData7(system, std::to_string(wellCentres), std::to_string(a)); }     // Writes data for option seven.
	else if (eightFlag) { writeData8(system, std::to_string(wellCentres)); }                                  // Writes data for option eight.
	
    // Clear all vectors for the next run.
    positions.clear();
    
    accRateTherm.clear();
    accRate.clear();
    E0.clear();
    E0Therm.clear();
    histogram.clear();
    Gx1x1.clear();
    Gx2x2.clear();
	instantons.clear();
	antiInstantons.clear();
	headerInfo.clear();
}

void metropolis(std::string system, int repeat, std::mt19937& rng, RepeatData& data, double (*potential)(double), double (*potentialDifferential)(double)) {
    // Set initial path and counters to 0.
    initialise(system, rng, data);

    // Thermalise the system.
    thermalise(potentialDifferential, potential, rng, data);

    takeMeasures(data.positions, data);                 // First measurement after thermalisation.

	// Take measures of the path every "decorrelation" sweeps.
    while (data.measureCount < measures) {
        metropolisSweep(potential, rng, data);
        data.sweep++;
        if (remainder(data.sweep, decorrSweeps) == 0) { takeMeasures(data.positions, data); }
    }
    // Compute the observables based on the measured paths.
	computeObservables(data.positionsVector, potentialDifferential, potential, data); 
    
}

void metropolisSweep(double (*potential)(double), std::mt19937& rng, RepeatData& data) {
    double newPosition;

    double y;
	double oldPosition;
	double leftPosition;
	double rightPosition;
    double kineticDelta;
    double potentialDelta;
	double actionDelta;
    
    // Left most point, wrap around to the right.
    y = uniformMinus1to1(rng);                                      // Set a random number between -1 and 1.
    newPosition = data.positions[0] + epsilon * y;                  // Increments one position by random number between -epsilon and +epsilon.
    oldPosition = data.positions[0];
    leftPosition = data.positions[N - 1];                           // Previous position.
    rightPosition = data.positions[1];                              // Next position.
    kineticDelta =
        (rightPosition - newPosition) * (rightPosition - newPosition)
        - (rightPosition - oldPosition) * (rightPosition - oldPosition)
        + (newPosition - leftPosition) * (newPosition - leftPosition)
        - (oldPosition - leftPosition) * (oldPosition - leftPosition);

    potentialDelta = potential(newPosition) - potential(oldPosition);

    actionDelta = (m * 0.5 * aInverse) * kineticDelta + a * potentialDelta;
    // Metropolis acceptance.
    if (actionDelta < 0 || uniform01(rng) < exp(-actionDelta)) {    // Accepted move.
        data.positions[0] = newPosition;
        data.acceptedMoves++;
    }

	// Update the middle points.
    for (int i = 1; i < N - 1; i++) {
        y = uniformMinus1to1(rng);                                  // Set a random number between -1 and 1.
        newPosition = data.positions[i] + epsilon * y;              // Increments one position by random number between -epsilon and +epsilon.
		oldPosition = data.positions[i];
		leftPosition = data.positions[i - 1];                       // Previous position.
		rightPosition = data.positions[i + 1];                      // Next position.
        kineticDelta =
            (rightPosition - newPosition) * (rightPosition - newPosition) 
            - (rightPosition - oldPosition) * (rightPosition - oldPosition)
            + (newPosition - leftPosition) * (newPosition - leftPosition)
            - (oldPosition - leftPosition) * (oldPosition - leftPosition);

        potentialDelta = potential(newPosition) - potential(oldPosition);

        actionDelta = (m * 0.5 * aInverse) * kineticDelta + a * potentialDelta;
        // Metropolis acceptance.
        if (actionDelta < 0 || uniform01(rng) < exp(-actionDelta)) { // Accepted move.
            data.positions[i] = newPosition;
            data.acceptedMoves++;
        }
    }

    // Right most point, wrap around to the left.
    y = uniformMinus1to1(rng);                                      // Set a random number between -1 and 1.
    newPosition = data.positions[N - 1] + epsilon * y;              // Increments one position by random number between -epsilon and +epsilon.
    oldPosition = data.positions[N - 1];
    leftPosition = data.positions[N - 2];                           // Previous position.
    rightPosition = data.positions[0];                              // Next position.
    kineticDelta =
        (rightPosition - newPosition) * (rightPosition - newPosition)
        - (rightPosition - oldPosition) * (rightPosition - oldPosition)
        + (newPosition - leftPosition) * (newPosition - leftPosition)
        - (oldPosition - leftPosition) * (oldPosition - leftPosition);

    potentialDelta = potential(newPosition) - potential(oldPosition);

    actionDelta = (m * 0.5 * aInverse) * kineticDelta + a * potentialDelta;
    // Metropolis acceptance.
    if (actionDelta < 0 || uniform01(rng) < exp(-actionDelta)) {    // Accepted move.
        data.positions[N - 1] = newPosition;
        data.acceptedMoves++;
    }

    // Update the window if user wanted visualisation.
    if (winRunning == true) {
        std::this_thread::sleep_for(std::chrono::milliseconds(delay));
    }
}

void initialise(std::string system, std::mt19937& rng, RepeatData& data) {
    // Reset the path.
    data.positions = std::vector<double>(N, 0.0);
    if (system == "DWP")                                            // DWP starts in one of the two wells.
    {
        data.positions = std::vector<double>(N, wellCentres * side);
		side *= -1;                                                 // Start the next run in the opposite well to symmetrise the wave function.
    }

	// Reset counters.
    data.sweep = 0;
    data.measureCount = 0;
    data.acceptedMoves = 0;

	// Initialise the path with random values between -maxDistance and maxDistance.
    if (hotStart) {
        for (int i = 0; i < N; i++) { data.positions[i] = uniformMinus1to1(rng) * maxDistance; }
    }
}

void thermalise(double (*potentialDifferential)(double), double (*potential)(double), std::mt19937& rng, RepeatData& data) {
    // Thermalisation performed before before measurements.
    while (data.sweep < thermSweeps) {
        metropolisSweep(potential, rng, data);
        data.sweep++;

        if ((data.sweep - 1) % thermInterval == 0) { takeThermMeasures(data.positions, potentialDifferential, potential, data); }
    }
}

void takeThermMeasures(std::vector<double>& positions, double (*potentialDifferential)(double), double (*potential)(double), RepeatData& data) {
    // Record acceptance rate.
    data.accRateThermTemp.push_back((double)(data.acceptedMoves) / (N * (double)thermInterval));
    data.acceptedMoves = 0;

    // Record ground state energy.
    data.E0ThermTemp.push_back(E0Calc(positions, potentialDifferential, potential));
}

void takeMeasures(std::vector<double>& positions, RepeatData& data) {
    // Record acceptance rate.
    data.accRateTemp.push_back((double)(data.acceptedMoves) / (N * (double)decorrSweeps));
    data.acceptedMoves = 0;

	// Record positions for observable calculations later.
    for (int i = 0; i < N; i++) {
        data.positionsVector[data.measureCount][i] = positions[i];
    }

    // Increment measure count.
    data.measureCount++;
}

void computeObservables(std::vector<std::vector<double>>& positionsVector, double (*potentialDifferential)(double), double (*potential)(double), RepeatData& data) {

    for (int i = 0; i < measures; i++) {
        // Record ground state energy.
        double E0 = E0Calc(positionsVector[i], potentialDifferential, potential);
        data.E0Temp.push_back(E0);


        // Record histogram data for the wavefunction.
        for (int t = 0; t < N; t++) {
            double x = positionsVector[i][t];

            int bin = int((x - xMin) / binWidth);

            if (bin >= 0 && bin < numBins) { data.histogramTemp[bin] += 1.0; }
        }
    }

    // Average the histogram over measures and normalise it.
    for (int i = 0; i < numBins; i++) { data.histogramTemp[i] /= (measures * N); }

    // Compute correlators.
    std::vector<std::vector<double>> x1(measures, std::vector<double>(N));
    std::vector<std::vector<double>> x2(measures, std::vector<double>(N));

    for (int measure = 0; measure < measures; measure++)
    {
        for (int t = 0; t < N; t++)
        {
            double x = positionsVector[measure][t];

            x1[measure][t] = x;
            x2[measure][t] = x * x;
        }
    }

    // Find vacuum point for x and x squared.
    double x1Vacuum = 0.0;
    double x2Vacuum = 0.0;

    for (int measure = 0; measure < measures; measure++)
    {
        for (int t = 0; t < N; t++)
        {
            x1Vacuum += x1[measure][t];
            x2Vacuum += x2[measure][t];
        }
    }

    x1Vacuum /= (measures * N);
    x2Vacuum /= (measures * N);

    for (int measure = 0; measure < measures; measure++)
    {
        std::vector<double> correlation11Temp;
        std::vector<double> correlation22Temp;

        correlation11Temp = correlator(x1[measure], x1[measure]);
        correlation22Temp = correlator(x2[measure], x2[measure]);

        for (int i = 0; i < N; i++) {
            data.Gx1x1Temp[i] += correlation11Temp[i];
            data.Gx2x2Temp[i] += correlation22Temp[i];
        }
    }

    // Average correlators and remove vacuum point.
    for (int n = 0; n < N; ++n) {
        data.Gx1x1Temp[n] /= measures;
        data.Gx1x1Temp[n] -= x1Vacuum * x1Vacuum;
        data.Gx2x2Temp[n] /= measures;
        data.Gx2x2Temp[n] -= x2Vacuum * x2Vacuum;
    }

    // Find instantons and anti instantons.
    for (int i = 0; i < measures; i++) {
        int prev = whichWell(positionsVector[i][0], tunnellingThreshold);

        int instantons = 0;
        int antiInstantons = 0;

        for (int j = 1; j < N; j++) {
            int curr = whichWell(positionsVector[i][j], tunnellingThreshold);

            if (curr == 0) { continue; }                                // Ignore barrier.

            if (prev == -1 && curr == 1) { instantons++; }

            if (prev == 1 && curr == -1) { antiInstantons++; }

            prev = curr;
        }

        // Check for tunnelling in the first position to account for the periodic boundary conditions.
        int curr = whichWell(positionsVector[i][0], tunnellingThreshold);

        if (curr != 0) { // Ignore barrier
            if (prev == -1 && curr == 1) { instantons++; }

            if (prev == 1 && curr == -1) { antiInstantons++; }
        }

        data.instantonsTemp.push_back(instantons);
        data.antiInstantonsTemp.push_back(antiInstantons);
    }
}