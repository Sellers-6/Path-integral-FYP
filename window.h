#pragma once

#define SDL_MAIN_HANDLED 

#include "global.h"     // Used for global variables which other headers need to access
#include <SDL2/SDL.h>
#include <algorithm>

// The main window drawing function
void window(const std::vector<double>& path, bool& runningFlag);

// Shared settings
extern int delay;       // Declared here, defined in window.cpp
extern bool winRunning;
extern double winSizeIncrement;
