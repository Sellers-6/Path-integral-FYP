//
// MAIN_Source.cpp
//

#include <gl/freeglut.h>
#include <iostream>
#include "IsingSystem.h"

bool running = true; // flag for whether the first loop is active
bool glutOn = false; // flag for whether the glut loop is active

using namespace std;

// functions which are needed for openGL go into a namespace so that we can identify them
namespace drawFuncs {
	void input1(unsigned char key, int x, char** y);
	void handleKeypress(unsigned char key, int x, int y);
	void display(void);
	void update(int val);
	void introMessage();
	void windowMessage();
	void quitMe(int val);
	void withWindow(int argc, char** argv);
}

// this is a global pointer, which is how we access the system itself
IsingSystem* sys;

int main(int argc, char** argv) {
	sys = new IsingSystem();

	// this is the seed for the random numbers
	int seed = 6;
	int GridSize = sys->returnGridSize();
	int GridSize3D = sys->return3DGridSize();
	cout << "Seed: " << seed << "; 2D grid size: " << GridSize << "; 3D grid size: " << GridSize3D << endl;
	sys->setSeed(seed);
	sys->readCsv(); // Finds the heat capacities for varying n and n_0

	drawFuncs::introMessage();


	while (running == true) {
		drawFuncs::input1(cin.get(), argc, argv);
	}

	cout << "Exiting main..." << endl;

	return 0;
}

void drawFuncs::withWindow(int argc, char** argv) {
	// turn on glut
	glutInit(&argc, argv);
	int window_size[] = { 480,480 };
	string window_title("2D square Ising simulation");

	Window* win = new Window(window_size, window_title);

	sys->createWindow(win);

	// tell openGL how to redraw the screen and respond to the keyboard
	glutDisplayFunc(drawFuncs::display);	
	glutKeyboardFunc(drawFuncs::handleKeypress); 

	// tell openGL to do its first update after waiting 10ms
	int wait = 10;
	int val = 0;
	glutTimerFunc(wait, drawFuncs::update, val);

	drawFuncs::windowMessage();

	glutOn = true;
	sys->glutOn = true;

	// start the openGL stuff
	glutMainLoop();
}

void drawFuncs::quitMe(int val) {
	cout << "delete system" << endl;
	delete sys;
	cout << "leave loop" << endl;
	if (glutOn == true) {
		glutLeaveMainLoop();
	}
}

void drawFuncs::introMessage() {
	cout << "To remove a possible bottleneck from openGL, the display is now optional." << endl;
	cout << "It is advised not to display the window when producing csv files, as the " << endl;
	cout << "openGL display does not update when corresponding functions are called." << endl;
	cout << "Note that for large grid sizes, csvs will take a long time to create." << endl;
	cout << "When csvs were produced for the report, six different seeds were ran simultaneously" << endl;
	cout << "for one sixth of the original measures to reduce computational time." << endl;
	cout << "Keys (type these and then hit ENTER):" << endl;
	cout << "  ? to print this message" << endl;
	cout << "  w to display the openGL window for the 2D square grid Ising model" << endl;
	cout << "  1-6 to create csv files for questions 1 - 8 respectively" << endl;
	cout << "  q to quit" << endl;
}

void drawFuncs::windowMessage() {
	cout << "Keys (while in graphics window):" << endl;
	cout << "  ? to print this message" << endl;
	cout << "  w to display the window" << endl;
	cout << "  q to quit" << endl;
	cout << "  f for fast; s for slow" << endl;
	cout << "  g to go; p to pause" << endl;
	cout << "  h for hotter; c for colder" << endl;
	cout << "  u for update once" << endl;
	cout << "  r to reset the system" << endl;
	cout << "  m to print current magentisation; e for current energy" << endl;
	cout << "  1-6 to create csv files for questions 1 - 8 respectively" << endl;
}

void drawFuncs::input1(unsigned char key, int argc, char** argv) {
	switch (key) {
	case 'w':
		withWindow(argc, argv);
		break;
	case 'q':
		drawFuncs::quitMe;
		running = false;
		break;
	case '?':
		drawFuncs::introMessage();
		break;
	case '1':
		sys->csvCreation(1);
		break;
	case '2':
		sys->csvCreation(2);
		break;
	case '3':
		sys->csvCreation(3);
		break;
	case '4':
		sys->csvCreation(4);
		break;
	case '5':
		sys->csvCreation(5);
		break;
	case '6':
		sys->csvCreation(6);
		break;
	}
}

// openGL function that deals with the keyboard
void drawFuncs::handleKeypress(unsigned char key, int x, int y) {
	switch (key) {
	case 'q':
		cout << "Quitting..." << endl;
		sys->pauseRunning();
		glutTimerFunc(500, drawFuncs::quitMe, 17);	// wait 500ms for any redrawing etc to finish, then quit
		break;
	case 's':
		cout << "slow" << endl;
		sys->setSlow();
		break;
	case 'f':
		cout << "fast" << endl;
		sys->setFast();
		break;
	case 'h':
		cout << "hotter" << endl;
		sys->Hotter();
		break;
	case 'c':
		cout << "colder" << endl;
		sys->Colder();
		break;
	case 'g':
		cout << "go" << endl;
		sys->setRunning();
		drawFuncs::update(0);
		break;
	case 'p':
		cout << "pause" << endl;
		sys->pauseRunning();
		break;
	case 'u':
		cout << "upd" << endl;
		sys->MCsweep();
		break;
	case 'r':
		cout << "reset" << endl;
		sys->Reset();
		break;
	case '?':
		drawFuncs::windowMessage();
		break;
	case 'm':
		cout << "Magnetism = " << sys->computeMagnetisation() << endl;
		break;
	case 'e':
		cout << "Energy = " << sys->computeEnergy() << endl;
		break;
	case 't':
		cout << "Number of sweeps performed so far = " << sys->sweepNumber << endl;
		break;
	case '1':
		sys->csvCreation(1);
		break;
	case '2':
		sys->csvCreation(2);
		break;
	case '3':
		sys->csvCreation(3);
		break;
	case '4':
		sys->csvCreation(4);
		break;
	case '5':
		sys->csvCreation(5);
		break;
	case '6':
		sys->csvCreation(6);
		break;
	}
	glutPostRedisplay();
}

// this function gets called whenever the algorithm should do its update
void drawFuncs::update(int val) {
	int wait;  // time to wait between updates (milliseconds)

	if (sys->isRunning()) {
		sys->MCsweep();

		glutPostRedisplay(); // tells GLUT that the display has changed
		// if you comment out this line then it will refresh
		// the screen only when you press a key

		if (sys->isSlow() == 1)
			wait = 200;
		else
			wait = 1;
		// tell openGL to call this funtion again after "wait" milliseconds
		glutTimerFunc(wait, drawFuncs::update, 0);
	}
}

// this function redraws the window when necessary
void drawFuncs::display() {
	//  Clear the window or more specifically the frame buffer...
	//  This happens by replacing all the contents of the frame
	//  buffer by the clear color (black in our case)
	glClear(GL_COLOR_BUFFER_BIT);

	// this puts the camera at the origin (not sure why) with (I think) z axis out of page and y axis up
	// there is also the question of the GL perspective which is not set up in any clear way at the moment
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0.0, 0.0, 1.0,   /* camera position */
		0.0, 0.0, -1.0,        /* point to look at */
		0.0, 1.0, 0.0);		   /* up direction */

	sys->DrawSquares();

	//  Swap contents of backward and forward frame buffers
	glutSwapBuffers();
}