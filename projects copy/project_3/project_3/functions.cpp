#include <cmath>
#include <iostream>
#include <armadillo>
#include "functions.h"
using namespace std;
using namespace arma;

void legendre_points_weights(vec& x, vec& w, const int n,
                             const double x1, const double x2)
{
    for(int ii = 1; ii <= (n+1)/2; ii++)
    {
        double z = cos(M_PI * (ii - 0.25)/(n + 0.5));
        double z1;
        double pp;
        do {
            double p1 = 1.0;
            double p2 = 0.0;


            for(int jj = 1; jj <= n; jj++)
            {
                double p3 = p2;
                p2 = p1;
                p1 = ((2.0 * jj - 1.0) * z * p2 - (jj - 1.0) * p3)/jj;
            }

            pp = n * (z * p1 - p2)/(z * z - 1.0);
            z1 = z;
            z  = z1 - p1/pp;                   // Newton's method
        } while(abs(z - z1) > 1e-10);

        double xm = 0.5 * (x2 + x1); // change of variables
        double xl = 0.5 * (x2 - x1);

        x(ii - 1) = xm - xl *z;
        x(n - ii) = xm + xl*z;

        w(ii - 1) = 2.0 * xl/((1.0 - z * z) * pp * pp);
        w(n - ii) = w(ii-1);
    }
    //cout << x << endl;
    //cout << w << endl;
}

double legendre_sum(vec& x, vec& w, const int N){
    //evaluate the integral with the Gauss-Legendre method
    //Note that we initialize the sum
    double int_gauss = 0.0;
    //six-double loops
    for (int i=0;i<N;i++){
        for (int j = 0;j<N;j++){
            for (int k = 0;k<N;k++){
                for (int l = 0;l<N;l++){
                    for (int m = 0;m<N;m++){
                        for (int n = 0;n<N;n++){
                            int_gauss+=w(i)*w(j)*w(k)*w(l)*w(m)*w(n)
                                    *int_function(x(i),x(j),x(k),x(l),x(m),x(n));
                        }}}}}}
    return int_gauss;
}

double spherical_sum(vec& xr, vec& wr, vec& xTheta, vec& wTheta,
                     vec& xPhi, vec& wPhi, const int N){
    //evaluate the integral with the Gauss-Legendre method
    //Note that we initialize the sum
    double int_gauss = 0.0;
    //six-double loops
    for (int i=0;i<N;i++){
        for (int j = 0;j<N;j++){
            for (int k = 0;k<N;k++){
                for (int l = 0;l<N;l++){
                    for (int m = 0;m<N;m++){
                        for (int n = 0;n<N;n++){
                            int_gauss+=wr(i)*wr(j)*wTheta(k)*wTheta(l)*wPhi(m)*wPhi(n)
                            *int_sphFunction(xr(i),xr(j),xTheta(k),xTheta(l),xPhi(m),xPhi(n));
                        }}}}}}
    return int_gauss;
}


//The integrand functions
double int_function(double x1, double y1, double z1, double x2, double y2, double z2)
{
    double alpha = 2.0;
    //Evaluate different parts of integrand function
    double exp1=-2*alpha*sqrt(x1*x1+y1*y1+z1*z1);
    double exp2=-2*alpha*sqrt(x2*x2+y2*y2+z2*z2);
    double deno=sqrt(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2));
    //Test if denominator is zero
    if(deno < pow(10.,-6.)) { return 0;}
    else return exp(exp1+exp2)/deno;
}

double int_sphFunction(double r1, double r2, double theta1,
                       double theta2, double phi1, double phi2)
{
    //Evaluate different parts of integrand function
    double numerator = sin(theta1)*sin(theta2);
    double cosBeta = (cos(theta1)*cos(theta2)) +
            (sin(theta1)*sin(theta2)*cos(phi1-phi2));
    double root = (r1*r1)+(r2*r2)-(2*r1*r2*cosBeta);
    double deno = sqrt(root);
    //Test if denominator is zero
    if(deno < pow(10.,-6.) || root < 0.0) { return 0;}
    else return numerator/deno;
}

double int_MCFunction(double r1, double r2, double theta1,
                       double theta2, double phi1, double phi2)
{
    //Evaluate different parts of integrand function
    double numerator = r1*r1*r2*r2*sin(theta1)*sin(theta2);
    double cosBeta = (cos(theta1)*cos(theta2)) +
            (sin(theta1)*sin(theta2)*cos(phi1-phi2));
    double root = (r1*r1)+(r2*r2)-(2*r1*r2*cosBeta);
    double deno = sqrt(root);
    //Test if denominator is zero
    if(deno < pow(10.,-6.) || root < 0.0) { return 0;}
    else return numerator/deno;
}

#define EPS 3.0e-14
#define MAXIT 10
void gauss_laguerre(vec& x, vec& w, int n, double alf)
{
    int i,its,j;
    double ai;
    double p1,p2,p3,pp,z,z1;

    for (i=1;i<=n;i++) {
        if (i == 1) {
            z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
        } else if (i == 2) {
            z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
        } else {
            ai=i-2;
            z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
                  (1.0+3.5*ai))*(z-x[i-3])/(1.0+0.3*alf);
        }
        for (its=1;its<=MAXIT;its++) {
            p1=1.0;
            p2=0.0;
            for (j=1;j<=n;j++) {
                p3=p2;
                p2=p1;
                p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
            }
            pp=(n*p1-(n+alf)*p2)/z;
            z1=z;
            z=z1-p1/pp;
            if (fabs(z-z1) <= EPS) break;
        }
        if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
        x(i-1)=z;
        w(i-1) = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
    }
}

double gammln( double xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
                          24.01409824083091,-1.231739572450155,
                          0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}


#undef EPS
#undef MAXIT
