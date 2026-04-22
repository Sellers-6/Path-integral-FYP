#pragma once

#include "hdf5.h"
#include <string>
#include <sstream>
#include <iostream>
#include "global.h"

/****************************************/
/*********** h5 data storage ************/
/****************************************/

// Write all simulation observables to h5 files in groups
void writeData(const std::string& system);
void writeData6(const std::string& system, const std::string& wellCentres, const std::string& beta);
void writeData7(const std::string& system, const std::string& wellCentres, const std::string& latticeSpacing);
void writeData8(const std::string& system, const std::string& wellCentres);