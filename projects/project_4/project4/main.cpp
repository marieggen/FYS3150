#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <functions.h>
#include "time.h"
#include "mpi.h"

using namespace std;
using namespace arma;
clock_t start, finish;

<<<<<<< HEAD
int main(int argc, char* argv[])
{
    int numprocs, my_rank;

    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    int N = 100;
    int Tsteps = 16;
    int MCCmax = 1e7;
    int MCClimit = 5000;


    vec numValues(5);
    start = clock();
    ising(N,Tsteps, MCCmax, MCClimit, my_rank, numprocs, numValues);
    finish = clock();
    double time_spent = (finish - start)/double(CLOCKS_PER_SEC);

    cout << "Time spent when N = " << N << ", Tsteps = " <<
            Tsteps << ", and #MCC = " << MCCmax+MCClimit <<
            " is t = " << time_spent << " sec." << endl;




    /*
    for(int j=2 ; j<10 ; j++){
        //Check change in variance
        //for different intarvals of
        //MC cycles.
        cout << "E interval " << j-1 << " and " << j <<
                ": " << VarE(j-1)/VarE(j) << endl;
        cout << "M interval " << j-1 << " and " << j <<
                " :" << VarM(j-1)/VarM(j) << endl;
    }

    double totE = MCC-MCClimit;

    for(int k=0 ; k<length ; k++){
        P_E(k) = P_E(k)/totE;
    }

    ofstream P_Efile("../project4/P_E_T1.txt");
    for(int i=0 ; i<length ; i++)
    {
        P_Efile << setprecision(16) << P_E(i) << endl;
    }
    P_Efile.close();*/

    MPI_Finalize ();

    //cout << "Evar " << EvarianceValues << endl;

=======



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

>>>>>>> parent of 874c528... d) almost done
    return 0;
}











