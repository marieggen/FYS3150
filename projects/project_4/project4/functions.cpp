#include <iostream>
#include <armadillo>
#include <cmath>
#include <cstdlib>
#include "lib.cpp"
#include <fstream>
#include <iomanip>

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
    values[0] = Z; values[1] = E; values[2] = Cv;
    values[3] = M; values[4] = X; values[5] = absM;
    values[6] = absX;

    return values;
}

int hei = 0;
void Metropolis(int nSpin, mat& s, double T, double& E, double& M,long& idum){
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

        double r = ran0(&idum);//test random num
        double w = exp(-deltaE/T);
        if((deltaE>(0.0) && (r<w)) || deltaE<=0){

            s(x,y) *= -1; //flip spin
            M += 2.0*s(x,y);
            E += deltaE;
        }
    }
    return;
}



void numMeanValues(double E, double M, vec& values){
    /*Find expectation values of energy,
    magnetization, specific heat capacity
    and susceptibility (numerical calc.)*/
    values(0) += E; values(1) += M;
    values(2) += E*E; values(3) += M*M;
    values(4) += abs(M);
    return;
}


vec outputValues(double T, int N, int MCc, vec& numValues){
    double norm = 1.0/((double) MCc); //divide by tot. num. of cycles
    double perSpin = 1.0/N*N;
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
    outputValues(5) = ((M2avg - absMavg*absMavg)/T)*perSpin;;//absXavg

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


void WriteToFile(ofstream &analyticalFile, ofstream &numericalFile, vec anaValues, vec outputVales, double T){
    //Printed as: E, Cv, M, X, |M|, |X|, T

    analyticalFile << setprecision(8) << anaValues(1) << setw(15) <<
             anaValues(2) << setw(15) << anaValues(3) <<
             setw(15) << anaValues(4) << setw(15) <<
             anaValues(5) << setw(15) << anaValues(6) <<
             setw(15) << T << endl;

    numericalFile << setprecision(8) << outputVales(0) << setw(15) <<
             outputVales(1) << setw(15) << outputVales(2) <<
             setw(15) << outputVales(3) << setw(15) <<
             outputVales(4) << setw(15) << outputVales(5) <<
             setw(15) << T << endl;

    return;
}

































