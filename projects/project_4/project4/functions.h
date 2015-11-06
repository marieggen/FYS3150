#ifndef FUNCTIONS
#define FUNCTIONS

#include <fstream>
#include <armadillo>
using namespace arma;
using namespace std;

vec anaMeanValues(double T);

void Metropolis(int nSpin, mat& s, double T, double& E, double& M,long& idum);

void numMeanValues(double E, double M, vec& values);

vec outputValues(double T, int N, int MCc, vec& numValues);

void print(vec anaValues, vec outputVales, double T);

void WriteToFile(ofstream &analyticalFile, ofstream &numericalFile, vec anaValues, vec outputVales, double T);



#endif // FUNCTIONS
