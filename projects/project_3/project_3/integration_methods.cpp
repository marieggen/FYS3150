#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <armadillo>
#include "functions.h"
#include <cmath>
#include "lib.h"
#include <time.h>

using namespace std;
using namespace arma;
#define pi M_PI

void legendre(const int N, const double a, const double b){

    //N = steps of integration (20)
    //a = lower limit on integral (-3.0)
    //b = upper limit on integral (3.0)

    //arrays for integration points and weights
    vec x(N);
    vec w(N);

    //Find mesh points and weights
    legendre_points_weights(x,w,N,a,b);
    double num_value = legendre_sum(x,w,N);

    //Analytically computed value
    double real_value = (5*pi*pi)/(16.0*16.0);

    cout << "Legendre value = " << num_value << endl;
    cout << "Real value = " << real_value << endl;
    cout << "Relative error = " << abs(num_value - real_value)/real_value
         << endl;

    return;
}



void laguerre(const int N){

    //N = steps of integration (20)

    //Constants
    double alf = 2.0; //Laguerre const.
    double alpha = 2.0; //charge particle

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

    //Find mesh points and weights
    gauss_laguerre(xu,wu, N, alf);//W(x)=x^(alf)*exp(-x)
    legendre_points_weights(xTheta,wTheta,N,thetaLim1,thetaLim2);//W(x=1)
    legendre_points_weights(xPhi,wPhi,N,phiLim1,phiLim2);//W(x)=1

    //Compute sum
    double num = spherical_sum(xu, wu, xTheta, wTheta,
                               xPhi, wPhi, N);

    //Integral when 1/length(r1-r2) = 1
    //double numSph_value = num/pow(2.0*alpha,6);
    //double realSph_value = (pi*pi)/64.0;

    //Integral when 1/length(r1-r2) != 1
    double numSph_value = num/pow(2.0*alpha,5);
    double realSph_value = (5*pi*pi)/(16.0*16.0);

    cout << "Laguerre value = " << numSph_value << endl;
    cout << "Real Value = " << realSph_value << endl;
    cout << "Relative error = " <<
            abs(numSph_value-realSph_value)/realSph_value << endl;

    return;
}


void bf_MC(int n, double y_max, double y_min){

    //n = number of simulations (10e6)
    //y_max = higher limit on integral (5.0)
    //y_min = lower limit on integral (-5.0)


    //Analytically calc. value
    double real_value = (5*pi*pi)/(16.0*16.0);

    //Def. of variables
    int num_var = 6;
    double interval = y_max-y_min;
    double jacobidet = pow(interval,num_var);
    vec y(num_var);

    double fy;
    double preMean = 0.0;
    double preVar = 0.0;
    long idum = time(0); //random seed

    //Run simulations
    for(int i=0 ; i<n ; i++){
        for(int j=0 ; j<num_var ; j++){
            double x_value = ran0(&idum);
            y(j) = (x_value*interval) + y_min;
        }

        fy = int_function(y(0),y(1),y(2),y(3),y(4),y(5));
        preMean += fy;
        preVar += fy*fy;
    }

    //Find mean, std.deviation and variance
    double mean = preMean/((double)n);
    double mean_sq = (preVar/(double)n);
    double variance = mean_sq - mean*mean;
    double std_dev = sqrt(variance);

    cout << setiosflags(ios::showpoint | ios::uppercase);
    cout << "Standard deviation = " << setw(10) <<
            setprecision(8) << jacobidet*(std_dev/sqrt((double)n)) << endl;
    cout << "MC integral = " << setw(10) <<
            setprecision(8) << jacobidet*mean << endl;
    cout << "Real value = " << setw(10) <<
            setprecision(8) << real_value << endl;
    cout << "Diff = " << setw(10) <<
            setprecision(8) << jacobidet*mean - real_value << endl;

    return;
}

void is_MC(int n){

    //n = number of simulations (10e6)

    //Analytically calc. value
    double real_value = (5*pi*pi)/(16.0*16.0);

    //Def. of constants and limits
    double alfa = 2.0;
    int num_var = 6;
    long idum = time(0); //random seed
    double theta_max = pi;
    double theta_min = 0;
    double phi_max = 2*pi;
    double phi_min = 0;

    double thetaInterval = theta_max-theta_min;
    double phiInterval = phi_max-phi_min;
    double thetaJacobi = pow(thetaInterval,2);
    double phiJacobi = pow(phiInterval,2);
    double jacobidetSph = (thetaJacobi*phiJacobi)/(pow(2.0*alfa,5));

    double r1, r2;
    vec x_value(num_var);

    double fr;
    double preSphMean = 0.0;
    double preSphVar = 0.0;

    //Run simulations
    for(int i=0 ; i<n ; i++){
        for(int j=0 ; j<num_var ; j++){
            x_value(j) = ran0(&idum);
        }

        double theta1 = (x_value(0)*thetaInterval) + theta_min;
        double theta2 = (x_value(1)*thetaInterval) + theta_min;
        double phi1= (x_value(2)*phiInterval) + phi_min;
        double phi2= (x_value(3)*phiInterval) + phi_min;
        double r1 = (-log(1-x_value(4)))/2*alfa;
        double r2 = (-log(1-x_value(5)))/2*alfa;

        double fr = int_MCFunction(r1,r2,theta1,theta2,phi1,phi2);
        preSphMean += fr;
        preSphVar += fr*fr;
    }

    //Find mean, std.deviation and variance
    double meanSph = preSphMean/((double)n);
    double mean_sqSph = (preSphVar/(double)n);
    double varianceSph = mean_sqSph - meanSph*meanSph;
    double std_devSph = sqrt(varianceSph);

    cout << setiosflags(ios::showpoint | ios::uppercase);
    cout << "Standard deviation = " << setw(10) <<
            setprecision(8) << jacobidetSph*(std_devSph/sqrt((double)n)) << endl;
    cout << "MC integral = " << setw(10) <<
            setprecision(8) << jacobidetSph*meanSph << endl;
    cout << "Real value = " << setw(10) <<
            setprecision(8) << real_value << endl;
    cout << "Diff = " << setw(10) <<
            setprecision(8) << jacobidetSph*meanSph - real_value << endl;

    return;
}

#undef pi
