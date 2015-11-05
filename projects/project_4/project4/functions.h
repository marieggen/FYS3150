#ifndef FUNCTIONS
#define FUNCTIONS

#endif // FUNCTIONS

#include <armadillo>
using namespace arma;

vec anaMeanValues(double T);

void Metropolis(int nSpin, mat& s, double T, double& E, double& M);

vec numMeanValues(double E, double M, vec& values);

vec outputValues(double T, int N, int MCc, vec& numValues);



