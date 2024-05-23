#ifndef SWAP_H
#define SWAP_H

// Libraries
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <sstream> 
#include <iomanip>
#include <vector> 

// Global variables
//  Run parameters
extern const int tau; //Correlation max-duration
extern const int steps; //Monte Carlo sweeps
extern const double T; //Temperature in units of 1/k_B

// Simulation parameters
const int N = 2000; //Number of particles
const double Size = 44.721359550000003; //Size of the system
const double sigmaMax = 1.613048; //Maximum diameter of particles
const double rSkin = 1.25 * sigmaMax; //Radius of neighbours included in NL (e.g. 1.8)
const double rC = 1.25 * sigmaMax; //Cutoff radius for calculating potential
const double rNL = pow(rC+rSkin,2); //NL radius squared
const double deltaMax = 0.12; //Max particle displacement
const double deltaSMax = 0.2; //Max diameter difference for swap
const double RUpdate = pow(rSkin,2)/4; //When R2Max exceeds this, update NL
const double x_max = 1.3; // maximal value of r/s for the real neighbours

// Some constants
const double c0 = -28/pow(1.25,12);
const double c2 = 48/pow(1.25,14);
const double c4 = -21/pow(1.25,16);
const double pi = 3.14159265358979323846;

// Arrays
extern double X[N], Y[N], S[N], X0[N], Y0[N];
extern double Xfull[N], Yfull[N], Xref[N], Yref[N], Xtw[N], Ytw[N];
// X0 initial position at last neighbour list update
// Xfull real positions (not taking into account periodic boundaries)
// Xref positition at t=0
// Xtw position at last aging update


//  Neighbour Lists
extern std::vector< std::vector<int> > NL;
extern std::vector< std::vector<int> > nn_0, nn_tw;
// nn_0 nearest neighbours at t=0
// nn_tw nearest neighbours at last aging update

//  Function prototypes
double bcs(double a, double b);
int Find(double arr[], int len, double seek);
void UpdateList(), UpdateAge();
double PairPotential(double x1, double y1, double s1, double x2, double y2, double s2);
double V(double xj, double yj, double rj, int j);
std::vector<int> effective_neighbours(int j), nearest_neighbours(int j, double x);
double VTotal(), CBLoc(int j), CB(), MSD(), FS(double theta);
void TryDisp(int j), TrySwap(int j, int k), MC(std::string out, int tw, int ss);

//  Random number between 0 and 1
#define ranf() \
    ((double)rand()/(1.0+RAND_MAX)) //check random numbers
    
#endif