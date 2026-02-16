#pragma once

#include <fstream>      // Used for file input and output
#include <map>          // Used for mapping file keys to file streams for output
#include <sstream>      // Used for string stream to efficiently write data to files in one go
#include <vector>
#include "global.h"


void writeFiles(std::string boundary, std::string system);

void writeVector(std::ofstream& file, std::ostringstream& buf, const std::vector<double>& data, const std::string& header);

void closeAllFiles();

void csvWriteData(std::string boundary, std::string system);

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

const std::vector<std::string> quantities = {   "E0Therm", 
                                                "E0Decorr", 
                                                "accRateTherm", 
                                                "accRateDecorr", 
                                                "GTherm", 
                                                "GDecorr", 
                                                "psiTherm",   
                                                "psiDecorr",
                                                "thermSweeps"};

using FileKey = std::tuple<std::string, Boundary, System>;
std::map<FileKey, std::ofstream> files;

void openFiles(std::string boundary, std::string system) {
    for (const auto& q : quantities) {

        std::string filename = baseDir + q + boundary + system + ".csv";

        files.emplace(  // Opens only the csv files we want to open
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

    // Stringstream buffer to write to files in one go, more efficient than writing line by line
    std::ostringstream bufferString;

    // References to the files we want to write to
    auto& E0TFile = files[{ "E0Therm", b, s }];      
    auto& E0DFile = files[{ "E0Decorr", b, s }];
    auto& accRateTFile = files[{ "accRateTherm", b, s }];
    auto& accRateDFile = files[{ "accRateDecorr", b, s }];
    auto& gTFile = files[{ "GTherm", b, s }];
    auto& gDFile = files[{ "GDecorr", b, s }];
    auto& psiTFile = files[{ "psiTherm", b, s }];
    auto& psiDFile = files[{ "psiDecorr", b, s }];

    // Write all the vectors to their respective files using the writeVector function
    writeVector(E0TFile, bufferString, E0Therm, "E0"); 
    writeVector(E0DFile, bufferString, E0Decorr, "E0");
    writeVector(accRateTFile, bufferString, accRateTherm, "AcceptanceRate");
    writeVector(accRateDFile, bufferString, accRateDecorr, "AcceptanceRate");
    writeVector(gTFile, bufferString, GTherm, "Correlation");
    writeVector(gDFile, bufferString, GDecorr, "Correlation");
    writeVector(psiTFile, bufferString, psiTherm, "Position");
    writeVector(psiDFile, bufferString, psiDecorr, "Position");

    // Seperate file writing for the wavefunction, as it is a 2D array rather than a vector
    bufferString.str("");       
    bufferString.clear();
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
        pair.second.close();    // Close all files
    }
}

void csvWriteData(std::string boundary, std::string system) {
    openFiles(boundary, system);
    writeFiles(boundary, system);   // Condense all the csv file writing steps into one clean function call for main.cpp
    closeAllFiles();
}