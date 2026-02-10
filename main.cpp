///// Simulating quantum systems with different potentials using Feynman path integral and Metropolis algorithm /////

#include "main.h"
#include "potentials.h"

int main() {    // Main function to run the Metropolis algorithm with user choice of boundary conditions and visualisation
    openAllFiles();
    std::cout << "Simulating quantum systems with different potentials using Feynman path integral and Metropolis algorithm" << std::endl;
    std::string choice;
	while (true) {
        chooseSystem();
		std::cout << "Enter choice: ";
        std::cin >> choice;
        if (choice == "1") {
            metropolis(false, 1, 1);
		}
        else if (choice == "2") {
            metropolis(true, 1, 1);
        }
        else if (choice == "3") {
            metropolis(false, 2, 1);
        }
        else if (choice == "4") {
            metropolis(true, 2, 1);
		}
        else if (choice == "5") {
            metropolis(false, 1, 2);
        }
        else if (choice == "6") {
            metropolis(true, 1, 2);
        }
        else if (choice == "7") {
            metropolis(false, 2, 2);
        }
        else if (choice == "8") {
            metropolis(true, 2, 2);
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

void metropolis(bool winOn, int BCs, int sys) { // Metropolis function which gets called initially, then calls other functions to perform the algorithm

    if (BCs == 1) {         // Periodic boundary conditions
        start = 0;
        end = N;
    }
    else if (BCs == 2) {    // Dirichlet boundary conditions
        start = 1;
		end = N - 1;        // Effectively fixes the endpoints to 0, as they are never updated
    }
    else {
        std::cout << "Invalid boundary condition choice." << std::endl;
        return;
    }

    if (sys == 0) {     // Free particle
        potential = FP::potential;
        potentialDifferential = FP::potentialDifferential;
    }
    else if (sys == 1) { // Quantum harmonic oscillator
        potential = QHO::potential;
        potentialDifferential = QHO::potentialDifferential;
    }
    else if (sys == 2) { // Double-well potential
        potential = DWP::potential;
        potentialDifferential = DWP::potentialDifferential;
    }
    else {
        std::cout << "Invalid system choice." << std::endl;
        return;
    }

	std::thread windowThread(window, std::ref(positions), std::ref(winRunning)); // Instantiates the window thread, passing the path and winRunning flag by reference
	winRunning = winOn; // Sets winRunning flag to true if user wanted visualisation

	initialise(sys);   // Sets initial parameters and path

    thermalise(winOn, potentialDifferential, potential);

    takeMeasures(positions, potentialDifferential, potential);
    std::cout << "Thermalisation complete after " << thermalisationRepeats << " sweeps" << std::endl;
	double percentMeasured = 0;
    bool percentChanged = true;
    while (measureCount < measures) {
        metropolisUpdate(winOn, potential);
        repeat++;
        if (remainder(repeat, decorrelation) == 0) {
            takeMeasures(positions, potentialDifferential, potential);
            while ((100 * measureCount) / measures > percentMeasured) {
                percentMeasured++;
                percentChanged = true;
            }
            if (percentChanged == true) {
                std::cout << percentMeasured << "% done." << std::endl;
            }
            percentChanged = false;
        }
        if (repeat % accRateInterval == 0) {
            acceptanceRate.push_back((double)(acceptedMoves - acceptedMovesPrevious) / (N * (double)accRateInterval));
            acceptedMovesPrevious = acceptedMoves;
        }
    }
    createFiles(BCs, sys);
    //createFiles2();
    std::cout << "Accepted moves: " << acceptedMoves << " out of " << repeat * N << ", giving a total acceptance rate of " << (double)acceptedMoves / ((double)repeat * N) << std::endl;
    //std::cout << "Action of final path: " << action << std::endl;
    double E0 = E0Avg / measures;
    std::cout << "Average ground state energy: " << E0 << std::endl;
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
    std::cout << "Energy of first exited state: " << E1 << std::endl;
	winRunning = false; // Sets winRunning flag to false to close the window thread
    windowThread.join();
    void closeAllFiles();
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
        pathUpdated = true;
        std::this_thread::sleep_for(std::chrono::milliseconds(delay));
    }
    return;
}

void initialise(int sys) {   // Initialises variables and path
    positions = std::vector<double>(N, 0.0);
    repeat = 0;
    acceptedMoves = 0;
    measureCount = 0;
    E0Avg = 0;
    thermalised = false;
	oldE0 = 0;
    thermalisationCheck = 0;
    //if (sys == 2) // DWP requires a different initial (cold) path since the potential wells are not centered around 0
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
        if (repeat % thermalisationMeasureInterval == 0) { // Store E0 during thermalisation for plotting purposes (not necessary for the algorithm itself)
            E0Thermalising.push_back((double)currentE0);
            thermalisationCount++;
        }
        if (repeat % accRateInterval == 0) {
            acceptanceRate.push_back((double)(acceptedMoves - acceptedMovesPrevious) / (double)accRateInterval);
            acceptedMovesPrevious = acceptedMoves;  // For measuring acceptance rate over sweeps
        }
        repeat++;
        if (E0Delta < thermalisationConstant) {
			thermalisationCheck++;  // If E0 is not changing much, increment thermalisation check counter
        }
        else {
			thermalisationCheck = 0; // If E0 is changing significantly, reset thermalisation check counter
        }
		if (thermalisationCheck >= thermalisationCheckLimit && repeat > thermalisationMinimum) { // If E0 has not changed significantly for a certain number of checks, we can assume thermalisation is complete
            thermalisationRepeats = repeat;
            thermalised = true;
        }
        if (repeat == thermalisationMaximum) {
            std::cout << "Thermalisation failed to converge after " << thermalisationMaximum << " repeats, proceeding with measurements anyway." << std::endl;
            thermalisationRepeats  = thermalisationMaximum;
            thermalised = true;
        }
		oldE0 = currentE0;
    }
    return;
}

void takeMeasures(std::vector<double> positions, double (*potentialDifferential)(double), double (*potential)(double)) {    // Takes all measurements at current path state in one function
	    
    double E0temp = E0Calc(positions, potentialDifferential, potential);             // Ground state measurements
    E0Avg += E0temp;
    E0Evolution.push_back((double)E0temp);

    for (int j = 0; j < N; j++) {
        psi[measureCount][j] = positions[j];        // Wavefunction measure (Data processed in python)
    }

    twoPointCorrelator(positions);

    measureCount++;
    return;
}

void createFiles(int BCs, int sys) {
    Boundary b;
	System s;
    if (BCs == 1) {
        b = Boundary::Periodic;
    }
    else if (BCs == 2) {
        b = Boundary::Dirichlet;
	}
    if (sys == 1) {
        s = System::QHO;
    }
    else if (sys == 2) {
        s = System::DWP;
    }
	std::ostringstream bufferString;        // Stringstream buffer to write to files in one go, more efficient than writing line by line

	auto& thermFile = files[{ "E0Thermalisation", b, s }];      // References to the files we want to write to
    auto& evoFile = files[{ "E0Evolution", b, s }];
    auto& accFile = files[{ "acceptanceRate", b, s }];
    auto& corrFile = files[{ "correlation", b, s }];
    auto& waveFile = files[{ "waveFunction", b, s }];

	writeVector(thermFile, bufferString, E0Thermalising, "E0"); // Write all the vectors to their respective files using the writeVector function
    writeVector(evoFile, bufferString, E0Evolution, "E0");
    writeVector(accFile, bufferString, acceptanceRate, "AcceptanceRate");
    writeVector(corrFile, bufferString, G, "Correlation");
    
	bufferString.str("");       // Seperate file writing for the wavefunction, as it is a 2D array rather than a vector
    bufferString.clear();
    for (int i = 0; i < measures; i++) {
        for (int j = 0; j < N; j++) {
            bufferString << psi[i][j];
            if (j != N - 1) bufferString << "\n";
        }
        bufferString << "\n";
    }
    waveFile << bufferString.str();

    return;
}

void writeVector(std::ofstream& file, std::ostringstream& buf, const std::vector<double>& data, const std::string& header) {
	file << header << "\n";     // Writes the header to the file
	buf.str("");                // Clears the stringstream buffer for reuse
    buf.clear();                // Clears any error flags on the stringstream
    for (double val : data)
		buf << val << "\n";     // Writes each value in the vector to the stringstream
	file << buf.str();          // Writes the entire contents of the stringstream to the file at once
}

void closeAllFiles() {
    for (auto& pair : files) {
        pair.second.close();    // Explicitly flush and close all files
    }
}

void window(const std::vector<double>& positions, bool& runningFlag) { 
    
    SDL_Init(SDL_INIT_VIDEO); // Initialises the window

    const int windowWidth = 800;
    const int windowHeight = 600;

    SDL_Window* window = SDL_CreateWindow(
        "Path Visualization",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, // x and y coordinates for centering the window
        windowWidth, windowHeight,  // Dimensions of created window
        0   // Passing 0 is default settings for the window, could decide to make it resizable, fullscreen, etc. here by changing the arguement
    );

    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED); // Creates a renderer, tied to this window. 
                                                                                       //-1 arguement lets SDL pick the renderer. SDL_RENDERER_ACCELERATED requests hardware acceleration

    SDL_Event evt; // evt will store information about a single event (key press, mouse move, window close) 

    while (runningFlag) { // Ensures window stays open while runningFlag true
        
        while (SDL_PollEvent(&evt)) {
            if (evt.type == SDL_QUIT)
                runningFlag = false;
        } // Checks whether window has been closed, sets runningFlag to false to indicate this

        auto minmax = std::minmax_element(positions.begin(), positions.end()); // Find min and max in the vector
        double minVal = *minmax.first;
        double maxVal = *minmax.second;
        if (maxVal - minVal < 1e-6) maxVal = minVal + 1e-6; // Avoids division by zero

        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255); // Sets the background of the window to be black
        if (pathUpdated) { // Neccessary check to ensure that the path does not update visually when it fails the metropolis probablity check
            SDL_RenderClear(renderer); // Clears the screen to the background colour
            SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255); // Sets the draw colour (for the path) to white

            for (int i = 0; i < (int)positions.size() - 1; i++) {
                int x1 = static_cast<int>(static_cast<unsigned long long>(i) * windowWidth / positions.size());
                int x2 = static_cast<int>((static_cast<size_t>(i) + 1) * windowWidth / positions.size()); // Maps the index of the position vector to the x coordinate

                double yNorm1 = (positions[i] - minVal) / (maxVal - minVal);
                double yNorm2 = (positions[static_cast<std::vector<double, std::allocator<double>>::size_type>(i) + 1] - minVal) / (maxVal - minVal); // Normalises the position of the path elements to (0,1), making the y axis vary dynamically with the range of position

                int y1 = windowHeight - static_cast<int>(yNorm1 * windowHeight);
                int y2 = windowHeight - static_cast<int>(yNorm2 * windowHeight); // Maps the position of the path elements to pixels on the window

                SDL_RenderDrawLine(renderer, x1, y1, x2, y2);   // Draws lines between the points 
            }

            SDL_RenderPresent(renderer); // Presents the drawn objects on the window
            pathUpdated = false;  // One draw function has been performed, the window will wait until it needs to update the path again
        }
        SDL_Delay(delay);  // Small delay to avoid overconsumption of cpu resources
    }

    SDL_DestroyRenderer(renderer);  // Frees up the render resources
    SDL_DestroyWindow(window);      // Quits the window
    SDL_Quit();                     // Quits SDL entirely
    return;
}