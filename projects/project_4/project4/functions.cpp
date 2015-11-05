#include <iostream>
#include <armadillo>
#include <cmath>
#include "lib.h"
#include <cstdlib>
using namespace std;
using namespace arma;

vec anaMeanValues(double T){
    /*All the values in this function
    is defined to be dimensionless.
    Find expectation values of energy,
    magnetization, specific heat capacity
    and susceptibility (analytically calc.)*/

    double beta = 1.0/T;
    double B = 8.0*beta;

    double Z = 4.0*cosh(B)+12.0;
    double E = (8.0*sinh(B))/(cosh(B) + 3.0);
    double Cv = -(((64.0)/(cosh(B) + 3.0))*
                  (cosh(B) + ((sinh(B)*sinh(B))/
                              (cosh(B) + 3.0))))/(T*T);
    double M = 0.0;
    double X = ((8.0*(exp(B) + 1.0))/(cosh(B) + 3.0))/(T);

    vec values(5);
    values[0] = Z; values[1] = E; values[2] = Cv;
    values[3] = M; values[4] = X;

    return values;
}


void Metropolis(int nSpin, mat& s, double T, double& E, double& M){
    /*Flip random spins, one at a time, and uses
     Metropolis to decide if we want to keep that
     change*/

    long idum = time(0);
    for(int i=0 ; i<(nSpin*nSpin) ; i++){
        //pick random spin
        int x = (int) (ran1(&idum)*(double)nSpin);
        int y = (int) (ran1(&idum)*(double)nSpin);
        s(x,y) *= -1; //flip spin
        double deltaE = 2*s(x,y)*
                (s(x,(y+nSpin+1)%nSpin)+
                 s(x,(y+nSpin-1)%nSpin)+
                 s((x+nSpin+1)%nSpin,y)+
                 s((x+nSpin-1)%nSpin,y));//change in energy

        if(deltaE>(-1.0)){
            double r = ran1(&idum);//what num?????
            double w = exp(-deltaE/T);
            if(r>w){
                s(x,y) *= -1;
                exit(1);
            }
        }
        M += 2.0*s(x,y);
        E += deltaE;
    }

    return;
}



vec numMeanValues(double E, double M, vec& values){
    /*Find expectation values of energy,
    magnetization, specific heat capacity
    and susceptibility (numerical calc.)*/
    values(0) += E; values(1) += M;
    values(2) += E*E; values(3) += M*M;
    values(4) += abs(M);
    return values;
}


vec outputValues(double T, int N, int MCc, vec& numValues){
    double norm = 1.0/((double) MCc); //divide by tot. num. of cycles
    double perSpin = 1.0/N*N;
    double E2avg = numValues(2)*norm*perSpin;
    double M2avg = numValues(3)*norm*perSpin;
    double Eavg = numValues(0)*norm*perSpin;
    double Mavg = numValues(1)*norm*perSpin;
    double absMavg = numValues(4)*norm*perSpin;

    vec outputValues(5);

    outputValues(0) = Eavg;//Eavg
    outputValues(1) = (E2avg - Eavg*Eavg)/(T*T);//CVavg
    outputValues(2) = Mavg;//Mavg
    outputValues(3) = (M2avg - Mavg*Mavg)/T;//Xavg
    outputValues(4) = absMavg;//absMavg

    return outputValues;
}















