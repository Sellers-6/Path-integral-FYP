#pragma once

#define SDL_MAIN_HANDLED	// Prevents SDL from overriding the main function, which can cause issues.

#include "global.h"
#include <SDL2/SDL.h>		// Used for visualisation of the Euclidean path.
#include <algorithm>

/************************************************************/
/*********** Visualisation of the Euclidean path ************/
/************************************************************/

// Window drawing function
void window(const std::vector<double>& path, bool& runningFlag);

// Shared settings
extern int delay;
extern bool winRunning;