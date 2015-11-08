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
    long idum = -clock();
    int N=2; //num. of spins in each dim.
    mat s(N,N);
    vec w(17); w.fill(0);

    double Tstart = 1.0;
    double Tfinish = 2.0;
    int Tsteps = 1;
    double Tstep = (Tfinish - Tstart)/(double) Tsteps;
    double T = Tstart;
    int MCC = 1e3; //#MC-cycles

    vec meanE(MCC);
    vec meanCv(MCC);
    vec meanM(MCC);
    vec meanX(MCC);
    vec MC_cycles(MCC);

    ofstream analyticalFile("../project4/ana_b.txt");
    ofstream numericalFile("../project4/num_b.txt");

    for(int i=0 ; i<Tsteps ; i++){
        double E = 0.0;
        double M = 0.0;
        initializeSpin(T, N, s, idum, E, M);
        for(int k=-8 ; k<=8 ; k+=4){w(k+8) = exp(-k/T);}

        vec numValues(5);
        numValues.fill(0);

        for(int cycle=0 ; cycle<MCC ; cycle++){

            Metropolis(w, N, s, T, E, M, idum);

            numMeanValues(E, M, numValues);
            //SaveAllValues(N, cycle, T, numValues, meanE, meanCv,
             //             meanM, meanX, MC_cycles, E, M);
        }

        vec outputVales = outputValues(T,N,MCC,numValues);
        vec anaValues = anaMeanValues(T, N);
        print(anaValues, outputVales,T);
        WriteToFileT(analyticalFile,numericalFile,anaValues,outputVales, T, MCC);
        //WriteToFileMCC(MCC, numericalFile, meanE, meanCv, meanM,
          //             meanX, MC_cycles);

        T += Tstep;
    }

    analyticalFile.close();
    numericalFile.close();

    return 0;
}











