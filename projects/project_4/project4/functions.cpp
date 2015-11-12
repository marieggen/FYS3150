#include <iostream>
#include <armadillo>
#include <cmath>
#include <cstdlib>
#include "lib.cpp"
#include <fstream>
#include <iomanip>

using namespace std;
using namespace arma;
/*
void WriteOutMatrix(mat s, ofstream outfile, int number){
    char *filename;
    sprintf(filename,"%d",number);
    s.save("hei.dat",raw_ascii);
    return;
}*/

vec anaMeanValues(double T, int N){
    /*All the values in this function
    is defined to be dimensionless.
    Find expectation values of energy,
    magnetization, specific heat capacity
    and susceptibility (analytically calc.)*/

    double beta = 1.0/T;
    double B = 8.0*beta;

    double Z = 4.0*cosh(B)+12.0;
    double E = -(8.0*sinh(B))/(cosh(B) + 3.0);
    double Cv = -(((64.0)/(cosh(B) + 3.0))*
                  (-cosh(B) + ((sinh(B)*sinh(B))/
                               (cosh(B) + 3.0))))/(T*T);
    double M = 0.0;
    double absM = (16.0 + (8.0*exp(B)))/Z;
    double X = ((8.0*(exp(B) + 1.0))/(cosh(B) + 3.0))/(T);
    double absX = ( ((8.0*(exp(B) + 1.0))/(cosh(B) + 3.0))
                    - (absM*absM) )/(T);

    vec values(7);
    values[0] = Z; values[1] = E/(N*N); values[2] = Cv/(N*N);
    values[3] = M/(N*N); values[4] = X/(N*N); values[5] = absM/(N*N);
    values[6] = absX/(N*N);

    return values;
}


void Metropolis(vec w, int nSpin, mat& s, double T, double& E, double& M,
                long& idum, int& count){
    /*Flip random spins, one at a time, and uses
     Metropolis to decide if we want to keep that
     change*/

    for(int i=0 ; i<(nSpin*nSpin) ; i++){
        //pick random spin
        int x = (int) (ran0(&idum)*(double)nSpin);
        int y = (int) (ran0(&idum)*(double)nSpin);

        double deltaE = 2*s(x,y)*
                (s(x,(y+nSpin+1)%nSpin)+
                 s(x,(y+nSpin-1)%nSpin)+
                 s((x+nSpin+1)%nSpin,y)+
                 s((x+nSpin-1)%nSpin,y));//change in energy

        if( ran0(&idum) <= w(deltaE+8) ){
            count += 1;
            s(x,y) *= -1; //flip spin
            M += 2.0*s(x,y);
            E += deltaE;
        }

    }
    return;
}



void numMeanValues(int N, double E, double M, vec& values,
                   vec& P_E){
    /*Find expectation values of energy,
    magnetization, specific heat capacity
    and susceptibility (numerical calc.)*/
    values(0) += E; values(1) += M;
    values(2) += E*E; values(3) += M*M;
    values(4) += abs(M);

    //cout << E+800 << endl;
    P_E(E + (2*N*N)) += 1;

    return;
}


vec outputValues(double T, int N, int MCc, vec& numValues,
                 vec& Evariance, vec& Mvariance, int numMCC,
                 vec& VarE, vec& VarM, int& num, vec& EvarianceValues){
    double norm = 1.0/((double) MCc); //divide by tot. num. of cycles
    double perSpin = 1.0/((double) N*(double) N);
    double Eavg = numValues(0)*norm;
    double Mavg = numValues(1)*norm;
    double E2avg = numValues(2)*norm;
    double M2avg = numValues(3)*norm;
    double absMavg = numValues(4)*norm;

    vec outputValues(6);

    outputValues(0) = Eavg*perSpin;//Eavg
    outputValues(1) = ((E2avg - Eavg*Eavg)/(T*T))*perSpin;//CVavg
    outputValues(2) = Mavg*perSpin;//Mavg
    outputValues(3) = ((M2avg - Mavg*Mavg)/T)*perSpin;//Xavg
    outputValues(4) = absMavg*perSpin;//absMavg
    outputValues(5) = ((M2avg - absMavg*absMavg)/T)*perSpin;//absXavg


    Evariance(numMCC) = outputValues(1)*(T*T);
    Mvariance(numMCC) = outputValues(3)*T;
    EvarianceValues(numMCC) = Evariance(numMCC)/perSpin;


    if(numMCC%100 == 0){

        VarE(num) = sum(Evariance);
        VarM(num) = sum(Mvariance);
        num += 1;

        Evariance.fill(0);
        Mvariance.fill(0);
    }

    return outputValues;
}


void print(vec anaValues, vec outputVales, double T){
    cout << "Temp: " << T << endl;
    cout << "    " << endl;

    cout << "num <E>: " << outputVales(0) << endl;
    cout << "num <Cv>: " << outputVales(1) << endl;
    cout << "num <M>: " << outputVales(2) << endl;
    cout << "num <X>: " << outputVales(3) << endl;
    cout << "num <|M|>: " << outputVales(4) << endl;
    cout << "num <|X|>: " << outputVales(5) << endl;
    cout << "    " << endl;

    cout << "ana <E>: " << anaValues(1) << endl;
    cout << "ana <Cv>: " << anaValues(2) << endl;
    cout << "ana <M>: " << anaValues(3) << endl;
    cout << "ana <X>: " << anaValues(4) << endl;
    cout << "ana <|M|> " << anaValues(5) << endl;
    cout << "ana <|X|> " << anaValues(6) << endl;
    cout << "    " << endl;
    return;
}


void WriteToFileT(ofstream &analyticalFile, ofstream &numericalFile, vec anaValues,
                  vec outputVales, double T, double MCC, int count){
    //Printed as: E, Cv, M, X, |M|, |X|, T, #MCC

    analyticalFile << setprecision(8) << anaValues(1) << setw(25) <<
             anaValues(2) << setw(25) << anaValues(3) <<
             setw(25) << anaValues(4) << setw(25) <<
             anaValues(5) << setw(25) << anaValues(6) <<
             setw(25) << T << setw(25) << MCC << endl;

    numericalFile << setprecision(8) << outputVales(0) << setw(25) <<
             outputVales(1) << setw(25) << outputVales(2) <<
             setw(25) << outputVales(3) << setw(25) <<
             outputVales(4) << setw(25) << outputVales(5) <<
             setw(25) << T << setw(25) << MCC << setw(15) << count << endl;

    return;
}


void initializeSpin(double T, int N, mat& s, long& idum,
                   double& E, double& M){
    if( T <= 1.5){
        s.fill(1); //spin matrix
        E = -2.0*N*N; //init value for low T
        M = N*N; // init value for low T
    }

    else{

        for( int x=0 ; x<N ; x++ ){
            for( int y=0 ; y<N ; y++ ){

                int component;
                if(ran0(&idum) < 0.5){
                    component = 1;
                }
                else{
                    component = -1;
                }
                s(x,y) = component;
            }
        }

        for( int x=0 ; x<N ; x++ ){
            for( int y=0 ; y<N ; y++ ){

                E -= (double) s(x,y)*
                        (s(x,(y+N-1)%N)+
                         s((x+N-1)%N,y));//change in energy
                M += s(x,y);

            }
        }
   }

    return;
}
/*
void SaveValues(int N, int cycle, vec numValues, vec& meanE, vec& Evalues,
                vec& MC_cycles, double E){
    //Saves the mean values from every monte
    //MC cycle into a matrix

    double norm = 1.0/((double) cycle); //divide by tot. num. of cycles
    double perSpin = 1.0/((double) N*(double) N);
    double Eavg = numValues(0)*norm;
    double E2avg = numValues(2)*norm;

    meanE(cycle) = Eavg*perSpin;//Eavg
    meanEvariance(cycle) = ((E2avg - Eavg*Eavg)/(T*T))*perSpin;//CVavg
    MC_cycles(cycle) = cycle;
    return;
}
*/





/*

void SaveAllValues(int N, int cycle, double T, vec& numValues, vec& meanE, vec& meanCv,
                   vec& meanM, vec& meanX, vec& MC_cycles, double E, double M){
    \\Saves the mean values from every monte
     \\MC cycle into a matrix
    numValues(0) += E; numValues(1) += M;
    numValues(2) += E*E; numValues(3) += M*M;
    numValues(4) += abs(M);

    double norm = 1.0/((double) cycle); //divide by tot. num. of cycles
    double perSpin = 1.0/((double) N*(double) N);
    double Eavg = numValues(0)*norm;
    double Mavg = numValues(1)*norm;
    double E2avg = numValues(2)*norm;
    double M2avg = numValues(3)*norm;
    double absMavg = numValues(4)*norm;

    meanE(cycle) = Eavg*perSpin;//Eavg
    meanCv(cycle) = ((E2avg - Eavg*Eavg)/(T*T))*perSpin;//CVavg
    meanM(cycle) = absMavg*perSpin;//absMavg
    meanX(cycle) = ((M2avg - absMavg*absMavg)/T)*perSpin;//absXavg
    MC_cycles(cycle) = cycle;
    return;
}


void WriteToFileMCC(int MCC, ofstream &numericalFile, vec& meanE, vec& meanCv,
                    vec& meanM, vec& meanX, vec& MC_cycles, vec count){
    //Printed as: E, Cv, M, X, |M|, |X|, #MCC
    for(int cycles=0 ; cycles<MCC ; cycles++){
        numericalFile << setprecision(8) << meanE(cycles) << setw(25) <<
                 meanCv(cycles) << setw(25) << meanM(cycles) <<
                 setw(25) << meanX(cycles) << setw(25) <<
                 MC_cycles(cycles) << setw(25) << count(cycles) << endl;
    }

    return;
}


            //WriteToFileMCC(MCC, numericalFile, meanE, meanCv, meanM,
              //             meanX, MC_cycles, countVec);
                //SaveAllValues(N, cycle, T, numValues, meanE, meanCv,
                  //            meanM, meanX, MC_cycles, E, M);

*/








