#pragma once
#include <vector>
#include <SDL2/SDL.h>
#include <algorithm>

// The main window drawing function
void window(const std::vector<double>& path, bool& runningFlag);

// Shared settings
extern int delay;       // Declared here, defined in window.cpp
extern bool winRunning; // Declared here, defined in window.cpp
