#pragma once

//#include "main.h"
#include <fstream>      // Used for file input and output
#include <map>          // Used for mapping file keys to file streams for output
#include <sstream>      // Used for string stream to efficiently write data to files in one go
#include <vector>
//#include "global.h"


void writeFiles(std::string boundary, std::string system);

void writeVector(std::ofstream& file, std::ostringstream& buf, const std::vector<double>& data, const std::string& header);

void closeAllFiles();

std::ofstream logfile;

const std::string baseDir = "csv/";

enum class Boundary { Periodic, Dirichlet };
enum class System { FP, QHO, DWP };

std::string boundaryName(Boundary b) {
    return (b == Boundary::Periodic) ? "Periodic" : "Dirichlet";
}

std::string systemName(System s) {
    switch (s) {
    case System::QHO: return "QHO";
    case System::DWP: return "DWP";
    case System::FP:  return "FP";
    default:          return "Unknown";
    }
}

const std::vector<std::string> quantities = {   "E0Thermalisation", 
                                                "E0Evolution", 
                                                "acceptanceRate", 
                                                "correlation", 
                                                "waveFunction", 
                                                "E0", 
                                                "E1",   
                                                "accRate" };

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
            std::make_tuple(
                q,
                Boundary(boundary == "Periodic" ? Boundary::Periodic : Boundary::Dirichlet),
                System(
                    system == "QHO" ? System::QHO :
                    (system == "DWP" ? System::DWP : System::FP)
                )
            ),
            std::ofstream(filename)
        );

    }
}


void writeFiles(std::string boundary, std::string system) {
    Boundary b;
    System s;
    if (boundary == "Periodic") {
        b = Boundary::Periodic;
    }
    else if (boundary == "Dirichlet") {
        b = Boundary::Dirichlet;
    }
    if (system == "FP") {
        s = System::FP;
    }
    else if (system == "QHO") {
        s = System::QHO;
    }
    else if (system == "DWP") {
        s = System::DWP;
    }
    std::ostringstream bufferString;        // Stringstream buffer to write to files in one go, more efficient than writing line by line

    auto& thermFile = files[{ "E0Thermalisation", b, s }];      // References to the files we want to write to
    auto& evoFile = files[{ "E0Evolution", b, s }];
    auto& accFile = files[{ "acceptanceRate", b, s }];
    auto& corrFile = files[{ "correlation", b, s }];
    auto& waveFile = files[{ "waveFunction", b, s }];
    auto& E0File = files[{ "E0", b, s }];
    auto& E1File = files[{ "E1", b, s }];
    auto& accRateFile = files[{ "accRate", b, s }];

    writeVector(thermFile, bufferString, E0Thermalising, "E0"); // Write all the vectors to their respective files using the writeVector function
    writeVector(evoFile, bufferString, E0Evolution, "E0");
    writeVector(accFile, bufferString, acceptanceRate, "AcceptanceRate");
    writeVector(corrFile, bufferString, G, "Correlation");
    writeVector(waveFile, bufferString, psi, "Position");
    writeVector(E0File, bufferString, E0Vec, "E0");
    writeVector(E1File, bufferString, E1Vec, "E1");
    writeVector(accRateFile, bufferString, accRateVec, "accRate");


    bufferString.str("");       // Seperate file writing for the wavefunction, as it is a 2D array rather than a vector
    bufferString.clear();

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