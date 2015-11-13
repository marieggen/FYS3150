#ifndef FUNCTIONS
#define FUNCTIONS

#include <fstream>
#include <armadillo>
using namespace arma;
using namespace std;

vec anaMeanValues(double T, int N);

void Metropolis(vec w, int nSpin, mat& s, double T, double& E, double& M,
                long& idum,int& count);

void numMeanValues(int N, double E, double M, vec& values,
                   vec& P_E);

vec outputValues(double T, int N, int MCc, vec& numValues,
                 vec& Evariance, vec& Mvariance, int numMCC,
                 vec& VarE, vec& VarM, int& num, vec& EvarianceValues);

void print(vec anaValues, vec outputVales, double T);

void WriteToFileT(ofstream &analyticalFile, ofstream &numericalFile, vec anaValues,
                  vec outputVales, double T, double MCC, int count);

void WriteOutMatrix(mat s, int number);

void initializeSpin(double T, int N, mat& s, long& idum,
                    double& E, double& M);

void SaveValues(int N, int cycle, vec numValues, vec& meanE, vec& Evalues,
                vec& MC_cycles, double E);

void ising(int N, int Tsteps, int MCCmax, int MCClimit,int my_rank,
          int numprocs, vec& numValues);

//void SaveAllValues(int N, int cycle, double T, vec& numValues, vec& meanE, vec& meanCv,
//vec& meanM, vec& meanX, vec& MC_cycles, double E, double M);

//void WriteToFileMCC(int MCC, ofstream &numericalFile, vec& meanE, vec& meanCv,
//vec& meanM, vec& meanX, vec& MC_cycles, vec count);

#endif // FUNCTIONS
