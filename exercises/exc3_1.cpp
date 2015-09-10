#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
using namespace std;

//typedef t double

ofstream ofile;

double f(double);
double df2(double, double);
double df3(double, double);
double df_real(double);


int main(int args, char* argv[])
{
	char *outfilename;

	int i, N = 10;
	double x = sqrt(2), h = 10.0;
	double error2[N-1], error3[N-1];
	double numFirstDerivative2;
	double numFirstDerivative3;
	double anaFirstDerivative;

	if( args <= 1 )
	{
		cout << "Bad Usage: " << argv[0] <<
		" ,read also two output files on same line." << endl;
		exit(1);
	}
	else
	{
		outfilename = argv[1];
	}


	ofile.open(outfilename);

	for(i = 1 ; i<=N ; i++)
	{
		double step = pow(h,-i);
		numFirstDerivative2 = df2(x,step);
		numFirstDerivative3 = df3(x,step);
		anaFirstDerivative = df_real(x);

		error2[i-1] = abs(anaFirstDerivative - numFirstDerivative2);
		error3[i-1] = abs(anaFirstDerivative - numFirstDerivative3);

		ofile << setw(15) << setprecision(8) << step << setw(15) << error2[i-1] << setw(15) << error3[i-1] << endl;
		cout << i << endl;
	} 
	ofile.close();

	return 0;
}


/*
ofile << "RESULTS:" << endl;
ofile << setiosflags(ios::showpoint | ios::uppercase);
ofile <<"R_min = " << setw(15) << setprecision(8) <<r_min <<endl; ofile <<"R_max = " << setw(15) << setprecision(8) <<r_max <<endl; ofile <<"Number of steps = " << setw(15) << max_step << endl; ofile << "Five lowest eigenvalues:" << endl;
for(i = 0; i < 5; i++) {
ofile << setw(15) << setprecision(8) << d[i] << endl;
*/



double f(double x)
{
	return atan(x);
}

double df2(double x, double h)
{
	double plussh = x + h;
	return (f(plussh) - f(x))/h;
}

double df3(double x, double h)
{
	double plussh = x + h;
	double minush = x - h;
	return (f(plussh)-f(minush))/(2.0*h);
}

double df_real(double x)
{
	return 1.0/(1+(x*x));
}







