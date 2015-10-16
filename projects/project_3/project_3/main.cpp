#include <iostream>
#include <cstdlib>
#include <ctime>
#include <armadillo>
#include "mpi.h"
#include "functions.h"
#include <cmath>
#include "lib.h"
#include <time.h>
using namespace std;
using namespace arma;
#define pi M_PI
int main()
{
    //steps of integration
    const int N = 20;

    //Legendre:
    //arrays for integration points and weights
    vec x(N);
    vec w(N);

    //Limits on integral
    const double a = -3.0;
    const double b = 3.0;

    //set up the mesh points and weights
    //    legendre_points_weights(x,w,N,a,b);
    //    double num_value = legendre_sum(x,w,N);

    //    //Compare computed and real value
    //    double real_value = (5*pi*pi)/(16.0*16.0);

    //    cout << num_value << endl;
    //    cout << real_value << endl;

    //Laguerre:
    //Const. (charge of particle)
    double alf = 2.0;
    double alpha = 2.0;

    //arrays for integration points and weights
    vec xu(N);
    vec wu(N);

    vec xTheta(N);
    vec wTheta(N);

    vec xPhi(N);
    vec wPhi(N);

    //Limits on integrals
    const double thetaLim1 = 0; //1.0;
    const double thetaLim2 = pi; //-1.0;

    const double phiLim1 = 0.0;
    const double phiLim2 = 2.0*pi;

//    gauss_laguerre(xu,wu, N, alf);//W(x)=x^(alf)*exp(-x)
//    legendre_points_weights(xTheta,wTheta,N,thetaLim1,thetaLim2);//W(x=1)
//    legendre_points_weights(xPhi,wPhi,N,phiLim1,phiLim2);//W(x)=1

//    double num = spherical_sum(xu, wu, xTheta, wTheta,
//                               xPhi, wPhi, N);
//    //double numSph_value = num/pow(2.0*alpha,6);
//    double numSph_value = num/pow(2.0*alpha,5);
//    //double realSph_value = (pi*pi)/64.0;
//    double realSph_value = (5*pi*pi)/(16.0*16.0);

//    cout << numSph_value << endl;
//    cout << realSph_value << endl;



    //Monte Carlo brute force
    double real_valueMC = (5*pi*pi)/(16.0*16.0);
    int n = 10000;
    double fx;
    double preMean = 0.0;
    double preVar = 0.0;
    long idum = time(0);

    for(int i=0 ; i<=n ; i++){
        double x1 = ran0(&idum);
        double y1 = ran0(&idum);
        double z1 = ran0(&idum);
        double x2 = ran0(&idum);
        double y2 = ran0(&idum);
        double z2 = ran0(&idum);

        fx = int_function(x1,y1,z1,x2,y2,z2);
        preMean += fx;
        preVar += fx*fx;
    }

    double mean = preMean/(double)n;
    double mean_sq = preVar/(double)n;
    double var = mean_sq - mean*mean;
    double std_dev = sqrt(var);

    cout << "standard deviation = " << std_dev << endl;
    cout << "integral = " << mean << endl;
    cout << "real value = " << real_valueMC << endl;

    //Improved Monte Carlo


    return 0;
}

#undef pi
















