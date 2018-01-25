#include<math.h>
#include<stdio.h>
#include<string>
#include<iostream>
#include<random>
#include<sstream>
#include<fstream>
#include<time.h>
#include<string.h>
#include<algorithm>
#include<stack>

using namespace std;

int destroyCount = 0;
int createCount = 0;

const int MAX_NUM_ATOMS = 600;
const int INIT_NUM_ATOMS = 500;
const float BOND = 4.0;//Bond strength
const float RATE_OF_GEN = 5.0;

const int RADIUS = 5;
const int NUM_RATES = 13;
const int ROWS = 100;
const int COLS = 100;
const int LAYS = 100;

const int R = 0;
const int C = 1;
const int L = 2;
const int N = 3;
const int Idx = 4;

const int EMPTY_SPACE = -1;
const int EDGE_SPACE = -2;
const int FLOOR_SPACE = -3;
const int GEN_ATOM = -4;

const float RATES[] = { 12,
11 * exp(-BOND),
10 * exp(-BOND * 2),
9 * exp(-BOND * 3),
8 * exp(-BOND * 4),
7 * exp(-BOND * 5),
6 * exp(-BOND * 6),
5 * exp(-BOND * 7),
4 * exp(-BOND * 8),
3 * exp(-BOND * 9),
2 * exp(-BOND * 10),
exp(-BOND * 11),
0 };

/*
all possible movement directions for FCC {row, col, lay}
To get opposite direction: opp_dir = (11-dir).
*/
const int DIR[12][3] = { { 1, 0, 0 },
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

const float PRIM[3][3] = { { 1.0, 0.0, 0.0 },
{ 0.5, sqrt(3.0) / 2.0, 0.0 },
{ 0.5, 1.0 / sqrt(12.0), sqrt(2.0 / 3.0) } };

int genPoints[8][3] = { { 1, 1, 1 },
{ ROWS - 1, 1, 1 },
{ 1, COLS - 1, 1 },
{ 1, 1, LAYS - 1 },
{ ROWS - 1, COLS - 1, 1 },
{ ROWS - 1, 1, LAYS - 1 },
{ 1, COLS - 1, LAYS - 1 },
{ ROWS - 1, COLS - 1, LAYS - 1 } };

int currentNumAtoms = 0;

/*
domain is the spatial field where the atoms 'physically' dwell
domain[x][y][z] = -2 represents a boarder space
domain[x][y][z] = -1 represents an empty space
domain[x][y][z] = [0 , MAX_NUM_ATOMS-1] represents a space occupied by an atom (identified by its number)
*/
int domain[ROWS][COLS][LAYS];

int picker[NUM_RATES - 1];

/*
inverselist is a list of all atoms in numerical order
information about each atom (row, col, lay, number of neighbors, index in atomlist) is stored in an array
*/
int inverselist[MAX_NUM_ATOMS][5];

/*
census stores the size of each row of atomlist
i.e. census[0] is equivalent to atomlist[0].size()
*/
int census[NUM_RATES];

/*
atomlist separates the atoms into arrays based on the number of neighbors an atom has
array[0] stores all atoms that have 0 neighbors. etc.
fully coordinated atoms (12 neighbors) are not stored in atomlist
*/
int atomlist[NUM_RATES][MAX_NUM_ATOMS];

stack<int> unusedAtoms;

float get_PN(){
	float sum = 0;
	for (int i = 0; i < NUM_RATES; i++){ sum += census[i] * RATES[i]; }
	return sum;
}

int pick_atom(float rand1, int rand2){
	float partial = 0;
	for (int i = 0; i < NUM_RATES - 1; i++){
		partial += census[i] * RATES[i];
		if (partial > rand1){
			return atomlist[i][rand2%census[i]];
		}
	}
	return GEN_ATOM;
}

int count_neighbors(int atom){
	int count = 0;
	for (int i = 0; i < NUM_RATES - 1; i++){
		int neighb = domain[inverselist[atom][R] + DIR[i][R]][inverselist[atom][C] + DIR[i][C]][inverselist[atom][L] + DIR[i][L]];
		if (neighb > EMPTY_SPACE){
			count++;
		}
	}
	return count;
}

int grab_neighbors(int atom){
	int count = 0;
	for (int i = 0; i < NUM_RATES - 1; i++){
		int neighb = domain[inverselist[atom][R] + DIR[i][R]][inverselist[atom][C] + DIR[i][C]][inverselist[atom][L] + DIR[i][L]];
		if (neighb > EMPTY_SPACE){
			picker[count] = i;
			count++;
		}
	}
	return count;
}

int grab_directions(int atom){
	int count = 0;
	for (int i = 0; i < NUM_RATES - 1; i++){
		int neighb = domain[inverselist[atom][R] + DIR[i][R]][inverselist[atom][C] + DIR[i][C]][inverselist[atom][L] + DIR[i][L]];
		if (neighb == EMPTY_SPACE){//if neighb == EMPTY_SPACE or EDGE_SPACE then put it in picker!
			picker[count] = i;
			count++;
		}
		else if (neighb == EDGE_SPACE){
			picker[count] = EDGE_SPACE; 
			count++;
		}
	}
	return count;
}

int pick_direction(int atom, int rand3){
	memset(picker, -1, sizeof(int)*(NUM_RATES - 1));
	int count = grab_directions(atom);
	if (count == 0) { return -1; }
	return picker[rand3%count];
}

void atom_change_list(int atom, int newlistnum){
	int oldlistnum = inverselist[atom][N];
	int oldindex = inverselist[atom][Idx];
	int swapped_atom = atomlist[oldlistnum][census[oldlistnum] - 1];
	atomlist[oldlistnum][oldindex] = swapped_atom;
	inverselist[swapped_atom][Idx] = oldindex;
	census[oldlistnum]--;
	atomlist[newlistnum][census[newlistnum]] = atom;
	inverselist[atom][N] = newlistnum;
	inverselist[atom][Idx] = census[newlistnum];
	census[newlistnum]++;
}

void old_neighbors(int atom){
	memset(picker, -1, sizeof(int)*(NUM_RATES - 1));
	int numneigh = grab_neighbors(atom);
	int r, c, l, neighbor, listn;
	for (int i = 0; i < numneigh; i++){
		r = DIR[picker[i]][R];
		c = DIR[picker[i]][C];
		l = DIR[picker[i]][L];
		neighbor = domain[inverselist[atom][R] + r][inverselist[atom][C] + c][inverselist[atom][L] + l];
		listn = inverselist[neighbor][N];
		atom_change_list(neighbor, listn - 1);
	}
}

void new_neighbors(int atom){
	memset(picker, -1, sizeof(int)*(NUM_RATES - 1));
	int numneigh = grab_neighbors(atom);
	int r, c, l, neighbor, listn;
	for (int i = 0; i < numneigh; i++){
		r = DIR[picker[i]][R];
		c = DIR[picker[i]][C];
		l = DIR[picker[i]][L];
		neighbor = domain[inverselist[atom][R] + r][inverselist[atom][C] + c][inverselist[atom][L] + l];
		listn = inverselist[neighbor][N];
		atom_change_list(neighbor, listn + 1);
	}
}

void destroy_atom(int atom){
	destroyCount++;
	old_neighbors(atom);
	unusedAtoms.push(atom);
	domain[inverselist[atom][R]][inverselist[atom][C]][inverselist[atom][L]] = EMPTY_SPACE;
	int oldlistnum = inverselist[atom][N];
	int oldindex = inverselist[atom][Idx];
	int swapped_atom = atomlist[oldlistnum][census[oldlistnum] - 1];
	atomlist[oldlistnum][oldindex] = swapped_atom;
	inverselist[swapped_atom][Idx] = oldindex;
	census[oldlistnum]--;
	memset(inverselist[atom], 0, sizeof(int) * 5);
	currentNumAtoms--;
}

void create_atom(int rand3){
	createCount++;
	int index = rand3 % 8;
	int pointr = genPoints[index][R];
	int pointc = genPoints[index][C];
	int pointl = genPoints[index][L];
	if (domain[pointr][pointc][pointl] == EMPTY_SPACE){
		int newAtom = unusedAtoms.top();
		unusedAtoms.pop();
		domain[pointr][pointc][pointl] = newAtom;
		inverselist[newAtom][R] = pointr;
		inverselist[newAtom][C] = pointc;
		inverselist[newAtom][L] = pointl;
		int newAtomList = count_neighbors(newAtom);
		inverselist[newAtom][N] = newAtomList;
		inverselist[newAtom][Idx] = census[newAtomList];
		atomlist[newAtomList][census[newAtomList]] = newAtom;
		census[newAtomList]++;
		new_neighbors(newAtom);
		currentNumAtoms++;
	}
}

void move_atom(int atom, int direction){
	if (direction < 0){
		//Comment this line to remove the flux boundary conditions
			if (direction == -2) destroy_atom(atom);
		//////////////////////////////////////////////////////////
		return;
	}
	old_neighbors(atom);
	int newrow = DIR[direction][R];
	int newcol = DIR[direction][C];
	int newlay = DIR[direction][L];
	domain[inverselist[atom][R] + newrow][inverselist[atom][C] + newcol][inverselist[atom][L] + newlay] = atom;
	domain[inverselist[atom][R]][inverselist[atom][C]][inverselist[atom][L]] = EMPTY_SPACE;
	inverselist[atom][R] = inverselist[atom][R] + newrow;
	inverselist[atom][C] = inverselist[atom][C] + newcol;
	inverselist[atom][L] = inverselist[atom][L] + newlay;
	int numneighbors = count_neighbors(atom);
	atom_change_list(atom, numneighbors);
	new_neighbors(atom);
}

void export_system(string filename, int currentNumAtoms){
	ofstream file;
	file.open(filename);
	stringstream ss;
	int count = 0;
	float x, y, z;
	float coor[3];
	for (int atom = 0; atom < currentNumAtoms; atom++){
		//if (inverselist[atom][3] == 0) continue;
		count++;
		x = (float)inverselist[atom][R];
		y = (float)inverselist[atom][C];
		z = (float)inverselist[atom][L];
		coor[0] = (x*PRIM[0][0] + y*PRIM[1][0] + z*PRIM[2][0]);
		coor[1] = (x*PRIM[0][1] + y*PRIM[1][1] + z*PRIM[2][1]);
		coor[2] = (x*PRIM[0][2] + y*PRIM[1][2] + z*PRIM[2][2]);
		ss << "H " << coor[0] << " " << coor[1] << " " << coor[2] << "\n";
	}
	file << to_string(count) << "\n\n" << ss.str();
	file.close();
}

void init_rand(){
	memset(atomlist, -1, sizeof(atomlist[0][0])*(NUM_RATES)*MAX_NUM_ATOMS);
	memset(census, 0, sizeof(census[0])*NUM_RATES);
	memset(picker, -1, sizeof(picker[0])*(NUM_RATES - 1));
	for (int i = MAX_NUM_ATOMS - 1; i > INIT_NUM_ATOMS - 1; i--) unusedAtoms.push(i);
	random_device randy;
	mt19937 rn(randy());
	uniform_int_distribution<int> ur(1, ROWS - 2);
	for (int i = 0; i < ROWS; i++){
		for (int j = 0; j < COLS; j++){
			for (int k = 0; k < LAYS; k++){
				domain[i][j][k] = (i == 0 || i == ROWS - 1 || j == 0 ||
					j == COLS - 1 || k == 0 || k == LAYS - 1) ? EDGE_SPACE : EMPTY_SPACE;
			}
		}
	}
	int count = 0;
	int r, c, l;
	while (count < INIT_NUM_ATOMS){
		r = ur(rn);
		c = ur(rn);
		l = ur(rn);
		if (domain[r][c][l] != -1) continue;
		domain[r][c][l] = count;
		inverselist[count][R] = r;
		inverselist[count][C] = c;
		inverselist[count][L] = l;
		count++;
	}
	int number;
	for (int atom = 0; atom < INIT_NUM_ATOMS; atom++){
		number = count_neighbors(atom);
		inverselist[atom][N] = number;
		if (number != NUM_RATES - 1){
			atomlist[number][census[number]] = atom;
			inverselist[atom][Idx] = census[number];
			census[number]++;
		}
	}
	currentNumAtoms = INIT_NUM_ATOMS;
}

void generateSphere(int atom, int& currentNumAtoms){
	int r, c, l;
	int sphereEdge = pow(RADIUS, 2);
	for (int i = 0; i < NUM_RATES - 1; i++){
		r = inverselist[atom][R] + DIR[i][R];
		c = inverselist[atom][C] + DIR[i][C];
		l = inverselist[atom][L] + DIR[i][L];
		if (domain[r][c][l] != -1) continue;
		if (pow(r-inverselist[0][R], 2) + pow(c-inverselist[0][C], 2) + pow(l-inverselist[0][L], 2) <= sphereEdge && currentNumAtoms < INIT_NUM_ATOMS){
			inverselist[currentNumAtoms][R] = r;
			inverselist[currentNumAtoms][C] = c;
			inverselist[currentNumAtoms][L] = l;
			domain[r][c][l] = currentNumAtoms;
			currentNumAtoms++;
			generateSphere(currentNumAtoms - 1, currentNumAtoms);
		}
	}
}

void init_circ(int& currentNumAtoms){
	memset(atomlist, -1, sizeof(atomlist[0][0])*(NUM_RATES)*MAX_NUM_ATOMS);
	memset(census, 0, sizeof(census[0])*NUM_RATES);
	memset(picker, -1, sizeof(picker[0])*(NUM_RATES - 1));
	for (int i = MAX_NUM_ATOMS-1; i > INIT_NUM_ATOMS-1; i--) unusedAtoms.push(i);
	for (int i = 0; i < ROWS; i++){
		for (int j = 0; j < COLS; j++){
			for (int k = 0; k < LAYS; k++){
				domain[i][j][k] = (i == 0 || i == ROWS - 1 || j == 0 ||
					j == COLS - 1 || k == 0 || k == LAYS - 1) ? EDGE_SPACE : EMPTY_SPACE;
			}
		}
	}
	int centerRow = (ROWS / 2) - 1;
	int centerCol = (COLS / 2) - 1;
	int centerLay = (LAYS / 2) - 1;
	domain[centerRow][centerCol][centerLay] = 0;
	inverselist[0][R] = centerRow;
	inverselist[0][C] = centerCol;
	inverselist[0][L] = centerLay;
	currentNumAtoms = 1;

	generateSphere(0, currentNumAtoms);

	int number;
	for (int atom = 0; atom < currentNumAtoms; atom++){
		number = count_neighbors(atom);
		inverselist[atom][N] = number;
		atomlist[number][census[number]] = atom;
		inverselist[atom][Idx] = census[number];
		census[number]++;
	}
}

string atomToString(int atom){
	return "Row: " + to_string(inverselist[atom][0]) + "\tCol: " + to_string(inverselist[atom][1]) + "\tLay: " +
		to_string(inverselist[atom][2]) + "\tList: " + to_string(inverselist[atom][3]) + "\tIndex: " +
		to_string(inverselist[atom][4]) + "\n";
}

int* get_atom(int atom){
	return inverselist[atom];
}

void printDomain(){
	for (int i = 0; i < ROWS; i++){
		for (int j = 0; j < COLS; j++){
			for (int k = 0; k < LAYS; k++){
				if (domain[i][j][k] == -1) cout << "  ";
				else cout << "* ";
			}
			cout << "\n";
		}
		cout << "\n";
	}
}

int main()
{
	init_circ(currentNumAtoms);
	//init_rand();
	srand(time(0));
	random_device rd;
	mt19937 rng(rd());
	uniform_real_distribution<double> urd(0.0, 1.0);
	clock_t t;
	t = clock();
	long double number;
	float rand1;
	int numIt = 200000000;
	unsigned int i;
	int rand2, rand3, atom, direction;

	for (i = numIt; i != -1; i--){
		if (i % 50000000 == 0) {
			cout << i << "\t" << unusedAtoms.size() << "\t" << destroyCount << "\t" << createCount << "\n";
			export_system("systemat" + to_string(i) + ".xyz", currentNumAtoms);
		}
		number = urd(rng) * 10000000000;
		rand1 = (number - (long long)number) * get_PN();
		rand2 = (long long)number % 100000;
		rand3 = (long long)number / 100000;
		atom = pick_atom(rand1, rand2);
		if (atom < 0){
			fprintf(stderr, "no atom picked: %d\n", i);
			continue;
		}
		/*if (atom == GEN_ATOM){
			if (!unusedAtoms.empty()) create_atom(rand3);
			continue;
		}*/
		direction = pick_direction(atom, rand3);
		/*if (direction < 0) {
			fprintf(stderr, "No Direction picked\n");
			cout << atom << " " << i << " " << atomToString(atom);
			continue;
		}*/
		move_atom(atom, direction);
	}
	t = clock() - t;
	printf("Time to iterate %d times: %d\n", i, t);
	//system("pause");
	return 0;
}