#pragma once

#include "hdf5.h"
#include <string>
#include <sstream>
#include <iostream>
#include "global.h"

// Writes all simulation observables to simulations.h5 into /Observable/Boundary/System groups
void writeData(const std::string& boundary, const std::string& system);