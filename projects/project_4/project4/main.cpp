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
clock_t start, finish;

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

    return 0;
}










