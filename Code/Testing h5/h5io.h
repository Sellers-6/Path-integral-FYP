#pragma once
#include <string>

// Writes all simulation observables to simulations.h5 into /Observable/Boundary/System groups
void writeData(const std::string& boundary, const std::string& system);