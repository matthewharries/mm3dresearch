//KMC Master Code
//3D atomic modeling program using the Kinetic Monte Carlo model
//Authors: Tim Krumweide, Tim Sitze, Hunter, Hyrum, Tanner Cranor

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <stdio.h>      
#include <stdlib.h>     
#include <time.h>
#include <cmath>
#include <stack>

#define n_max 12

using namespace std;

//Global Constants
const int MAXLEDGER = 10000000;
const int MAXGRID = 100;
const int X = 0;
const int Y = 1;
const int Z = 2;

//Global Variables
int ledger[10000000][5];				//contains all atoms and their information
										//name of atom is its index in the ledger (this helps with space and is accounted for the the generation and deletion of atoms using the index stack)
										//ledger[atom#][0] - x coordinate
										//ledger[atom#][1] - y coordinate
										//ledger[atom#][2] - z coordinate
										//ledger[atom#][3] - neighborhood within census
										//ledger[atom#][4] - index within neighborhood (of census)

stack<int> ghostatoms;

const int NEIGHBORLIST[12][3] = //Hyrum's order to make neighborfinding more efficient 
{ { 1, 0, 0 },
{ 0, 1, 0 },
{ 1, -1, 0 },
{ 0, 0, 1 },
{ 1, -1, 1 },
{ -1, 1, 1 },
{ 1, -1, -1 },
{ -1, 1, -1 },
{ 0, 0, -1 },
{ -1, 1, 0 },
{ 0, -1, 0 },
{ -1, 0, 0 } };

int gridarray[100][100][100];
int census[12][10000000];
int censusinterpreter[12];
double censusrates[12];																				//census sum gone, kaput. //
//censusinterpreter[i][0] - # of atoms currently in each neighborhood of the census
//censusrates[i][1] - rate: (12 - i)*exp((0 - a) * i)
//censusrates[i][2] - sum: (censusinterpreter[i]*censusrates[i])

bool atomfloor = false;
int numatoms = 0;

//end Global Variables


void updateAtomNCount(int atom) {
	int atomx = ledger[atom][0];
	int atomy = ledger[atom][1];
	int atomz = ledger[atom][2];
	int neighborx;
	int neighbory;
	int neighborz;
	int n = 0;

	for (int i = 0; i < 12; i++) {
		neighborx = atomx + NEIGHBORLIST[i][0];
		neighbory = atomy + NEIGHBORLIST[i][1];
		neighborz = atomz + NEIGHBORLIST[i][2];
		if (gridarray[neighborx][neighbory][neighborz] > -1) {
			n++;
		}
	}
	ledger[atom][3] = n;
	return;
}

void initializeAtominCensus(int atom) {
	int neighborhood = ledger[atom][3];							//atomn to neighborhood//
	int atomindex;

	for (int i = 0; i < numatoms; i++) {
		atomindex = censusinterpreter[neighborhood];			
		census[neighborhood][atomindex] = i;
		ledger[atom][4] = atomindex;
		censusinterpreter[neighborhood] += 1;
	}
	return;
}

void initializeAtom(int xin, int yin, int zin) {
	ledger[numatoms][X] = xin;										//Updated the initialize to the correct ledger system//
	ledger[numatoms][Y] = yin;
	ledger[numatoms][Z] = zin;
	gridarray[xin][yin][zin] = numatoms;
	numatoms++;
	return;
}

void generateSphereWithRadius(int r) {

	const int sphereBounds = pow(r, 2); // I believe Tim said we'd use a radius of 10, but I wanted to make the function more universally applicable, so I allow for the input of a radius.
	const int sphereCenter = 49; //should be (GRID_SIZE/2 - 1)...I suggest we use more global constants to avoid so many magic numbers
								 //initialize center atom
	initializeAtom(sphereCenter, sphereCenter, sphereCenter);

	int currentAtom = 0;
	while (currentAtom < numatoms) { //while there are still new/unanalyzed atoms in ledger...
		const int x = ledger[currentAtom][X];
		const int y = ledger[currentAtom][Y];
		const int z = ledger[currentAtom][Z];

		//cycle through each neighboring position of currentAtom
		for (int i = 0; i < 12; i++) {

			//get neighbor atom coordinates
			const int newX = x + NEIGHBORLIST[i][X];
			const int newY = y + NEIGHBORLIST[i][Y];
			const int newZ = z + NEIGHBORLIST[i][Z];

			//if already defined as atom, continue to next neighboring coordinates
			if (gridarray[newX][newY][newZ] >= 0) {
				continue;
			}

			//uses following formula to determine which positions should be atoms in sphere:
			// (x - sphereCenter)^2 + (y - sphereCenter)^2 + (z - sphereCenter)^2 <= r^2
			if ((pow((newX - sphereCenter), 2) + pow((newY - sphereCenter), 2) + pow((newZ - sphereCenter), 2)) <= sphereBounds) {
				initializeAtom(newX, newY, newZ);
			}
		}

		currentAtom++;
	}
}

void initializeArray(bool solidfloor, int startradius, int a) {
	//Initialize ledger to absolutely blank (-1's)
	for (int i = 0; i < 10000000; i++) {
		for (int j = 0; j < 5; j++) {
			ledger[i][j] = -1;
		}
	}
	//Initialize Census to absolutely blank
	for (int i = 0; i < 12; i++) {
		for (int j = 0; j < 10000000; j++) {
			census[i][j] = -1;
		}
	}
	//Initialize Census Interpreter
	for (int i = 0; i < 12; i++) {
		censusinterpreter[i] = 0;												//update the censusrates and censusinterpreter to correct format// 
		censusrates[i] = (12 - i)*exp((0 - a) * i);
		//censusrates[i][2] = 0;
	}
	for (int x = 0; x < 100; x++) {
		for (int y = 0; y < 100; y++) {
			for (int z = 0; z < 100; z++) {
				if ((y < 3)) {
					gridarray[x][y][z] = -3;
				}
				else if (((x < 3) || (y < 3) || (z < 3) || (x > 96) || (y > 96) || (z > 96)) && (solidfloor)) {
					gridarray[x][y][z] = -2;
				}
				else {
					gridarray[x][y][z] = -1;
				}
			}
		}
	}
	generateSphereWithRadius(startradius);
	for (int i = 0; i < numatoms; i++) {
		updateAtomNCount(i);
	}
	//Run update n and fill census
	for (int i = 0; i < numatoms; i++) {
		updateAtomNCount(i);
	}
	for (int i = 0; i < numatoms; i++) {
		initializeAtominCensus(i);
	}

	return;
}

unsigned int pickAtom() {
	unsigned int p_sum = 0;

	double census_sum[n_max];														//calculate census_sum here//
	for (int i = 0; i < n_max; i++) {
		census_sum[i] = (censusinterpreter[i] * censusrates[i]);
	}

	for (unsigned int i = 0; i < n_max; i++) p_sum += census_sum[i];
	// don't foget to seed the random number generator inside main()
	double r = ((double)rand() / (RAND_MAX));
	p_sum = p_sum * r;
	unsigned int i = 0;
	while (p_sum > census_sum[i]) {
		p_sum -= census_sum[i];
		if (i < n_max) ++i;
	}
	p_sum = p_sum / census_sum[i];
	return census[i][p_sum];
}

int selectDirection(int atom, int rand3) {
	int count = 0;
	int selector[12];
	memset(selector, -1, sizeof(int) * 12);
	for (int i = 4; i < 16; i++) {
		if (ledger[atom][i] == -1) {
			selector[count] = i;
			count++;
		}
	}
	if (count == 0) return -1;
	return selector[rand3%count] - 4;
}

void atomChangeList(int atom, int newlistnum) {
	int oldlistnum = ledger[atom][3];
	int oldindex = ledger[atom][4];
	int swappedAtom = census[oldlistnum][censusinterpreter[oldlistnum] - 1];
	census[oldlistnum][oldindex] = swappedAtom;
	ledger[swappedAtom][4] = oldindex;
	censusinterpreter[oldlistnum]--;
	census[newlistnum][censusinterpreter[newlistnum]] = atom;
	ledger[atom][3] = newlistnum;
	ledger[atom][4] = censusinterpreter[newlistnum];
	censusinterpreter[newlistnum]++;
}

void updateOldNeighbors(int atom) {
	int newn = 0;
	int x_ = 0, y_ = 0, z_ = 0;
	for (int i = 0; i < 12; i++) {
		x_ = ledger[atom][0] + NEIGHBORLIST[i][0];
		y_ = ledger[atom][1] + NEIGHBORLIST[i][1];
		z_ = ledger[atom][2] + NEIGHBORLIST[i][2];
		if (gridarray[x_][y_][z_] > -1) {
			newn = ledger[atom][3] - 1;
			atomChangeList(gridarray[x_][y_][z_], newn);
		}
	}
}

void updateNewNeighbors(int atom) {
	int newn = 0;
	int x_ = 0, y_ = 0, z_ = 0;
	for (int i = 0; i < 12; i++) {
		x_ = ledger[atom][X] + NEIGHBORLIST[i][X];
		y_ = ledger[atom][Y] + NEIGHBORLIST[i][Y];
		z_ = ledger[atom][Z] + NEIGHBORLIST[i][Z];
		if (gridarray[x_][y_][z_] > -1) {
			newn = ledger[atom][3] + 1;
			atomChangeList(gridarray[x_][y_][z_], newn);
		}
	}
}

void moveAtom(int atom, int direction) {
	//update old neighbors																			//updated all the 0,1,2 to XYZ format//
	updateOldNeighbors(atom);

	//Move atom
	//When flux boundary is introduced we will need to add a case for when the selected direction 
	//takes the atom off the edge of the domain. When that happens, the remaining code in
	//this function is ignored and the 'killAtom' function is called on 'atom'.
	int newrow = NEIGHBORLIST[direction][X];
	int newcol = NEIGHBORLIST[direction][Y];
	int newlay = NEIGHBORLIST[direction][Z];
	gridarray[ledger[atom][X] + newrow][ledger[atom][Y] + newcol][ledger[atom][Z] + newlay] = atom;
	gridarray[ledger[atom][X]][ledger[atom][Y]][ledger[atom][Z]] = -1;
	ledger[atom][X] = ledger[atom][X] + newrow;
	ledger[atom][Y] = ledger[atom][Y] + newcol;
	ledger[atom][Z] = ledger[atom][Z] + newlay;
	int numneighbors = 0;
	for (int i = 0; i < 12; i++) {
		int x_ = 0, y_ = 0, z_ = 0;
		x_ = ledger[atom][X] + NEIGHBORLIST[i][Y];
		y_ = ledger[atom][Y] + NEIGHBORLIST[i][Y];
		z_ = ledger[atom][Z] + NEIGHBORLIST[i][Z];
		if (gridarray[x_][y_][z_] > -1) {
			numneighbors++;
		}
	}
	atomChangeList(atom, numneighbors);

	updateNewNeighbors(atom);
}

//void generateAtom(){
//return;
//}

void AtomDestruct(int atom) {														//This is my code, feel free to look at it and make comments//
	int newn = 0;
	int x = 0, y = 0, z = 0;
	updateOldNeighbors(atom);

	gridarray[ledger[atom][X]][ledger[atom][Y]][ledger[atom][Z]] = -1;

	atomChangeList(atom, 0);

	ledger[atom][X] = -1;
	ledger[atom][Y] = -1;
	ledger[atom][Z] = -1;
	ledger[atom][3] = -1;
	ledger[atom][4] = -1;

	ghostatoms.push(atom);
}

void printAtoms(clock_t t1, int iterations) {
	ofstream fout;
	clock_t t2 = clock();
	float diff((float)t2 - (float)t1);
	float seconds = diff / CLOCKS_PER_SEC;
	cout << iterations << " iterations" << endl << endl;
	cout << seconds << " seconds" << endl << endl;
	cout << iterations / seconds << " iterations per second" << endl << endl;

	fout.open("Random.xyz");
	fout << numatoms;
	fout << endl << "KMC" << endl;

	for (int atom = 0; atom < numatoms; atom++) {
		fout << "X" << " ";
		fout << ledger[atom][X] << " ";
		fout << ledger[atom][Y] << " ";
		fout << ledger[atom][Z] << " ";
		fout << 0.5 << endl;
	}

	fout.close();
}

//int calculateTime(){
//return 0;
//}

//int trajectory(){
//return 0;
//}

int main() {
	int a;
	int solidground;
	int startradius;
	int iterations = 0;
	double iterrand = 0;
	int atompicked = 0;
	int directpicked = 0;
	cout << "Welcome to KMC.  " << endl;
	cout << "Please enter a value for the bond strength \'a\':";
	cin >> a;
	cout << "\nPlease choose an initialization environment" << endl;
	cout << "0 - space (no floor)" << endl;
	cout << "1 - solid floor" << endl;
	cin >> solidground;
	cout << "\nPlease choose an initial sphere radius" << endl;
	cout << "(integer from 2 - 10" << endl;

	cin >> startradius;

	initializeArray(true, startradius, a);

	cout << "Grid Initialized!" << endl;

	srand(time(NULL));

	while (iterations != -1) {
		cout << "Please enter the deisred number of iterations or -1 to quit: ";
		cin >> iterations;
		if (iterations == -1) {
			break;
		}
		clock_t t1 = clock();
		while (iterations > 0) {
			iterrand = rand();

			atompicked = pickAtom();
			directpicked = selectDirection(atompicked, iterrand);
			moveAtom(atompicked, directpicked);

			iterations--;
		}
		printAtoms(t1, iterations);
	}



	system("pause");
	return 0;
}