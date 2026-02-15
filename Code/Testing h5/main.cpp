///// Simulating quantum systems with different potentials using Feynman path integral and Metropolis algorithm /////

#include "main.h"
#include "potentials.h"
#include "h5io.h"


int main() {    // Main function to run the Metropolis algorithm with user choice of boundary conditions and visualisation
    std::cout << "Simulating quantum systems with different potentials using Feynman path integral and Metropolis algorithm" << std::endl;
    std::string choice;
    while (true) {
        chooseSystem();
        std::cout << "Enter choice: ";
        std::cin >> choice;
        if (choice == "1") {
            metropolisRepeat(false, 1, 1);
        }
        else if (choice == "2") {
            metropolisRepeat(true, 1, 1);
        }
        else if (choice == "3") {
            metropolisRepeat(false, 2, 1);
        }
        else if (choice == "4") {
            metropolisRepeat(true, 2, 1);
        }
        else if (choice == "5") {
            metropolisRepeat(false, 1, 2);
        }
        else if (choice == "6") {
            metropolisRepeat(true, 1, 2);
        }
        else if (choice == "7") {
            metropolisRepeat(false, 2, 2);
        }
        else if (choice == "8") {
            metropolisRepeat(true, 2, 2);
        }
        else if (choice == "A") {
            metropolisRepeat(false, 1, 1);
            std::cout << "QHO with periodic boundary conditions complete." << std::endl;
            metropolisRepeat(false, 2, 1);
            std::cout << "QHO with Dirichlet boundary conditions complete." << std::endl;
            metropolisRepeat(false, 1, 2);
            std::cout << "DWP with periodic boundary conditions complete." << std::endl;
            metropolisRepeat(false, 2, 2);
            std::cout << "DWP with Dirichlet boundary conditions complete." << std::endl;
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
    std::cout << "A: Perform all tests (without path visualisation)" << std::endl;
    std::cout << "0: Exit" << std::endl;
    return;
}

void metropolisRepeat(bool winOn, int BCs, int sys) {
    std::string boundary;
    std::string system;
    if (BCs == 1) {         // Periodic boundary conditions
        start = 0;
        end = N;
        boundary = "Periodic";
    }
    else if (BCs == 2) {    // Dirichlet boundary conditions
        start = 1;
        end = N - 1;        // Effectively fixes the endpoints to 0, as they are never updated
        boundary = "Dirichlet";
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
        system = "QHO";
    }
    else if (sys == 2) { // Double-well potential
        potential = DWP::potential;
        potentialDifferential = DWP::potentialDifferential;
        system = "DWP";
    }
    else {
        std::cout << "Invalid system choice." << std::endl;
        return;
    }

    std::thread windowThread(window, std::ref(positions), std::ref(winRunning)); // Instantiates the window thread, passing the path and winRunning flag by reference
    winRunning = winOn; // Sets winRunning flag to true if user wanted visualisation

    for (int repeat = 0; repeat < repeats; repeat++) {
        metropolis(winOn, BCs, sys, repeat);
    }
    winRunning = false; // Sets winRunning flag to false to close the window thread
    windowThread.join();

    writeSimulationRun(boundary, system);


    
    // Remember to reset the vectors here

    return;
}

void metropolis(bool winOn, int BCs, int sys, int repeat) { // Metropolis function which gets called initially, then calls other functions to perform the algorithm
    initialise(sys);   // Sets initial parameters and path

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
    //std::cout << "Accepted moves: " << acceptedMoves << " out of " << sweep * N << ", giving a total acceptance rate of " << accRate << std::endl;
    double E0 = E0Avg / measures;
    //std::cout << "Average ground state energy: " << E0 << std::endl;
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
    //std::cout << "Energy of first exited state: " << E1 << std::endl;
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
    sweep = 0;
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
            thermalisationSweeps = thermalisationMaximum;
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
        psi[measureCount + (repeat * measures)][j] = positions[j];
    }

    std::vector<double> corr = twoPointCorrelator(positions);
    for (int j = 0; j < N; j++) {
        G[j + (repeat * N)] += corr[j];
    }

    measureCount++;
    return;
}

//void writeFiles(std::string boundary, std::string system) {
//    //std::map<std::string, std::vector<double>> observableVecs;
//    initializeObsVecMap();
//
//
//    hid_t file_id = H5Fopen("simulations.h5", H5F_ACC_RDWR, H5P_DEFAULT);
//    if (file_id < 0) {
//        file_id = H5Fcreate("simulations.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
//    }
//
//    writeRun(system, boundary);
//
//    H5Fclose(file_id);
//}

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