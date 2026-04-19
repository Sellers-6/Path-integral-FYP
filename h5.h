#pragma once

#include "hdf5.h"
#include <string>
#include <sstream>
#include <iostream>
#include "global.h"

// Writes all simulation observables to simulations.h5 into groups
void writeData(const std::string& system);
void writeData2(const std::string& system, const std::string& wellCentres, const std::string& beta);
void writeData3(const std::string& system, const std::string& latticeSpacing);