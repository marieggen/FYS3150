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
    int N=20; //num. of spins in each dim.
    mat s(N,N);
    vec w(17); w.fill(0);

    double Tstart = 1.0;
    double Tfinish = 2.0;
    int Tsteps = 1;
    double Tstep = (Tfinish - Tstart)/(double) Tsteps;
    double T = Tstart;

    double MCCstart = 0.0;
    double MCCfinish = 5000.0;
    int MCCsteps = 5000;
    double MCCstep = (MCCfinish - MCCstart)/(double) MCCsteps;
    double MCC = MCCstart;
    //int MCC = 1e2; //#MC-cycles

    vec meanE(MCC);
    vec meanCv(MCC);
    vec meanM(MCC);
    vec meanX(MCC);
    vec MC_cycles(MCC);

    ofstream analyticalFile("../project4/ana_mcc.txt");
    ofstream numericalFile("../project4/num_mcc.txt");

    for(int i=0 ; i<Tsteps ; i++){
        double E = 0.0;
        double M = 0.0;
        initializeSpin(T, N, s, idum, E, M);
        vec anaValues = anaMeanValues(T, N);
        for(int k=-8 ; k<=8 ; k+=4){w(k+8) = exp(-k/T);}

        for(int cycles=0 ; cycles<=MCCsteps ; cycles++){

            int count = 0;
            vec numValues(5);
            numValues.fill(0);

            for(int cycle=0 ; cycle<MCC ; cycle++){
                Metropolis(w, N, s, T, E, M, idum, count);

                numMeanValues(E, M, numValues);
                //SaveAllValues(N, cycle, T, numValues, meanE, meanCv,
                  //            meanM, meanX, MC_cycles, E, M);
            }

            vec outputVales = outputValues(T,N,MCC,numValues);
            //print(anaValues, outputVales,T);
            WriteToFileT(analyticalFile,numericalFile,anaValues,outputVales,
                         T, MCC, count);
            //WriteToFileMCC(MCC, numericalFile, meanE, meanCv, meanM,
              //             meanX, MC_cycles, countVec);

            MCC += MCCstep;

        }

        T += Tstep;
    }

    analyticalFile.close();
    numericalFile.close();

    return 0;
}











