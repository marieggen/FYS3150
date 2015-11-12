#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <functions.h>
#include "time.h"
#include "mpi.h"

using namespace std;
using namespace arma;

int main(int argc, char* argv[])
{
    int numprocs, my_rank;

    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    long idum = -clock();
    int N=20; //num. of spins in each dim.
    mat s(N,N);
    vec w(17); w.fill(0);

    double Tstart = 1.0;
    double Tfinish = 2.0;
    int Tsteps = 1;
    double Tstep = (Tfinish - Tstart)/(double) Tsteps;
    double T = Tstart;

//    int NoOfTempPerCPU = Tsteps/numprocs;
//    int Tbegin =

    double T_interval = Tsteps/(double) numprocs;
    double Tbegin = Tstart + (T_interval*my_rank*Tstep);
    double Tend = Tstart + (T_interval*(my_rank+1)*Tstep);
    int T_count = 0;

    mat localValues(Tsteps,6);
    localValues.fill(0);
    mat totalValues(Tsteps,6);
    localValues.fill(0);

    int MCClimit = 0;
    double MCCstart = 0.0;
    double MCCfinish = 5000;
    int MCCsteps = 1;
    double MCCstep = (MCCfinish - MCCstart)/(double) MCCsteps;
    double MCC = 1e5+MCClimit;
    int numMCC = 0;
    int num = 0;

    vec Evalues(MCC);
    vec meanE(MCC);
    vec meanCv(MCC);
    vec meanM(MCC);
    vec meanX(MCC);
    vec MC_cycles(MCC);

    vec Evariance(MCCsteps); Evariance.fill(0);
    vec Mvariance(MCCsteps); Mvariance.fill(0);
    vec EvarianceValues(MCCsteps);
    vec VarE(10);
    vec VarM(10);
    int length = 2*N*N+2*N*N;
    vec P_E(length);

    ofstream analyticalFile("../project4/ana_test.txt");
    ofstream numericalFile("../project4/num_test.txt");



    //for(int i=0 ; i<Tsteps ; i++){
    for(double T=Tbegin ; T<Tend ; T+=Tstep){
        double E = 0.0;
        double M = 0.0;
        initializeSpin(T, N, s, idum, E, M);
        vec anaValues = anaMeanValues(T, N);
        for(int k=-8 ; k<=8 ; k+=4){w(k+8) = exp(-k/T);}

        for(int cycles=0 ; cycles<MCCsteps ; cycles++){

            int count = 0;
            vec numValues(5);
            numValues.fill(0);

            for(int cycle=MCClimit ; cycle<MCC ; cycle++){
                Metropolis(w, N, s, T, E, M, idum, count);
                numMeanValues(N, E, M, numValues,P_E);
            }
            cout << "Temperature: " << T << endl;
            cout << numValues << endl;



            vec outputVales = outputValues(T,N,MCC,numValues,
                                           Evariance, Mvariance,numMCC,
                                           VarE, VarM, num, EvarianceValues);
            //print(anaValues, outputVales,T);
            //WriteToFileT(analyticalFile,numericalFile,anaValues,outputVales,
            //           T, MCC, count);

            MCC += MCCstep;
            numMCC += 1;
        }

        // Write everything to file

        T_count += 1;
    }




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
    analyticalFile.close();
    numericalFile.close();
    MPI_Finalize ();

    //cout << "Evar " << EvarianceValues << endl;

    return 0;
}










