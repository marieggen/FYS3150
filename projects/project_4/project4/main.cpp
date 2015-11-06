#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <functions.h>
#include "time.h"

using namespace std;
using namespace arma;




int main()
{
    long idum = -1;
    int N=2; //num. of spins in each dim.
    mat s(N,N); s.fill(1); //spin matrix
    double E = -8.0; //init value for low T
    double M = 4.0; // init value for low T

    double Tstart = 1.0;
    double Tfinish = 4.0;
    int Tsteps = 1000;
    double Tstep = (Tfinish - Tstart)/(double) Tsteps;
    double T = Tstart;

//    double MCCstart = 1.0;
//    double MCCfinish = 2.0;
//    int MCCsteps = 100;
//    double MCCstep = (MCCfinish - MCCstart)/(double) steps;
//    double MCC = MCCstart;

    int MCc = 1e6; //#MC-cycles

    ofstream analyticalFile("../project4/ana.txt");
    ofstream numericalFile("../project4/num.txt");

    for(int i=0 ; i<Tsteps ; i++){
        vec numValues(5);
        numValues.fill(0);
        for(int cycle=0 ; cycle<MCc ; cycle++){
            Metropolis(N, s, T, E, M, idum);
            numMeanValues(E, M, numValues);
        }
        vec outputVales = outputValues(T, N,MCc,numValues);
        vec anaValues = anaMeanValues(T);
        //print(anaValues, outputVales,T);
        WriteToFile(analyticalFile,numericalFile,anaValues,outputVales, T);

        T += Tstep;
    }

    analyticalFile.close();
    numericalFile.close();

    return 0;
}











