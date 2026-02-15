#pragma once
#include <vector>


//// Header file for lattice parameters ////


const int N = 5000;												// Number of lattice points
std::vector<double> positions = std::vector<double>(N, 0.0);	// Lattice points
const double a = 0.1;											// Lattice spacing (I don't understand why this must be set to ~ 0.01 - 0.1)