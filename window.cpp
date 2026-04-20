#include "window.h"

/************************************************************/
/*********** Visualisation of the Euclidean path ************/
/************************************************************/

int delay = 0;                                          // Delay in milliseconds between frames.
bool winRunning = false;

// Window visualisation magic, feel free to change windowWidth/windowHeight and minVal/maxVal, everything else should be kept the same.

void window(const std::vector<double>& positions, bool& runningFlag) {
    SDL_Init(SDL_INIT_VIDEO);

    const int windowWidth = 800;
    const int windowHeight = 600;

    SDL_Window* window = SDL_CreateWindow(
        "Path Visualization",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        windowWidth, windowHeight,
        0
    );

    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

    SDL_Event evt;
    double yMax = 0.0;
    while (runningFlag) {
        while (SDL_PollEvent(&evt)) {
            if (evt.type == SDL_QUIT)
                runningFlag = false;
        }

		// Fixed y-axis range based on the expected values of the path, adjust for different potentials.
		double minVal = -6.0;
		double maxVal = 6.0;

        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);

        for (size_t i = 0; i < positions.size() - 1; i++) {
            int x1 = static_cast<int>(i * windowWidth / positions.size());
            int x2 = static_cast<int>((i + 1) * windowWidth / positions.size());

            double yNorm1 = (positions[i] - minVal) / (maxVal - minVal);
            double yNorm2 = (positions[i + 1] - minVal) / (maxVal - minVal);

            int y1 = windowHeight - static_cast<int>(yNorm1 * windowHeight);
            int y2 = windowHeight - static_cast<int>(yNorm2 * windowHeight);

            SDL_RenderDrawLine(renderer, x1, y1, x2, y2);
        }
        SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
        SDL_RenderDrawLine(renderer, 0, windowHeight / 2, windowWidth, windowHeight / 2);
        SDL_RenderPresent(renderer);
        SDL_Delay(delay);
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
}
