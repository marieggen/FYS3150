#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <armadillo>
#include "functions.h"
#include <cmath>
#include "lib.h"
#include "time.h"
#include "integration_methods.h"

clock_t start, finish;

int main()
{

    start = clock();
    legendre(40,-3.0,3.0);
    finish = clock();
    double time_legendre = (finish - start)/double(CLOCKS_PER_SEC);
    cout << "Time spent by legendre = " << time_legendre << endl;
/*
    start = clock();
    laguerre(30);
    finish = clock();
    double time_laguerre = (finish - start)/double(CLOCKS_PER_SEC);
    cout << "Time spent by laguerre = " << time_laguerre << endl;

    start = clock();
    bf_MC(10e7,5.0,-5.0);
    finish = clock();
    double time_bfMC = (finish - start)/double(CLOCKS_PER_SEC);
    cout << "Time spent by brute force Monte Carlo = "
         << time_bfMC << endl;

    start = clock();
    is_MC(10e7);
    finish = clock();
    double time_isMC = (finish - start)/double(CLOCKS_PER_SEC);
    cout << "Time spent by importance sampling Monte Carlo = "
         << time_isMC << endl;
*/
    return 0;
}

