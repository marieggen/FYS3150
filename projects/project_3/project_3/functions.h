#ifndef FUNCTIONS
#define FUNCTIONS

#endif // FUNCTIONS

#include <armadillo>
using namespace arma;

double gammln(double);

double spherical_sum(vec& xr, vec& wr, vec& xTheta, vec& wTheta,
                     vec& xPhi, vec& wPhi, const int N);

void gauss_laguerre(vec& x, vec& w, int n, double alf);

double int_function(double x1, double y1, double z1,
                    double x2, double y2, double z2);

void legendre_points_weights(vec& x, vec& w, const int n,
                             const double x1, const double x2);

double legendre_sum(vec& x, vec& w, const int N);

double int_sphFunction(double r1, double r2,
                       double theta1, double theta2, double phi1, double phi2);
