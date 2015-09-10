#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "time.h"
using namespace std;

ofstream ofile;
clock_t start, finish;
double f(double x);
double u_analytic(double x);
void solver_double_derivative(int n, char* outfilename, char* outfilename_error);

int main(int args, char** argv)
{

		if( args <= 3 )
	{
		cout << "Bad Usage " << argv[0] <<
		": read also output file on same line." << endl;
		exit(1);
	}
	else
	{
		char* outfilename = argv[1];
		char*outfilename_error = argv[2];
		int n = atoi(argv[3]); // atoi = ascii to int, atof = ascii to float
		solver_double_derivative(n, outfilename, outfilename_error);

	}
		return 0;
}




double f(double x)
{
	//analytical expression f(x)
	return 100.0*exp(-10.0*x);
}

double u_analytic(double x)
{
	//analytical expression u(x)
	return 1-(1-exp(-10))*x-exp(-10*x);
}

void solver_double_derivative(int n, char*outfilename, char*outfilename_error)
{
	//function that writes solution data of u(x) and x to a file

	int i;
	double h = (1.0/(n+1));
	double h_sq = h*h;

	double* a = new double[n];
	double* b = new double[n];
	double* c = new double[n];
	double* x = new double[n+2];
	double* y = new double[n];
	double* u = new double[n+2];
	double* u_ana = new double[n+2];

	double* errorContainer = new double[n+2];
	//double* steps = new double[n+2];

	for(i=0;i<=(n-1);i++)
		{
			//fill vectors
			x[i+1] = 0.0 + (i+1)*h;
			a[i] = 1.0;
			b[i] = -2.0;
			c[i] = 1.0;
			y[i] = -h_sq*f(x[i+1]);
		}
		start = clock();
		//boundary conditions on x
		x[0] = 0;
		x[n+1] = 1;

		for(i=1;i<=(n-1);i++)
		{
			//forward substitution
			double factor = (1.0/b[i-1]);
			b[i] = b[i] - factor;
			y[i] = y[i] - (factor*y[i-1]);
			a[i-1] = 0;
		}

		//boundary conditions on u(x)
		u[n+1] = 0;
		u[0] = 0;

		for(i=n;i>=1;i--)
		{
			//fill u(x) by backward substitution
			u[i] = (y[i-1]-u[i+1])/b[i-1];
		}
		finish = clock();
		double time_spent = (finish - start)/(double)CLOCKS_PER_SEC;

		for(i=0;i<=(n+1);i++)
		{
			//fill vector with analytical values of u(x)
			u_ana[i] = u_analytic(x[i]);
		}

		ofile.open(outfilename);
		for(i=0;i<=(n+1);i++)
		{
			//write data to file
			ofile << setprecision(8) << x[i] << setw(15) << u[i] << setw(15) << u_ana[i] << setw(15) << time_spent << endl;
		}
		ofile.close();


		for(i=1;i<=(n+1);i++)
		{
			double absError = abs((u[i]-u_ana[i])/u_ana[i]);
			errorContainer[i] = log10(absError);
		}

		ofile.open(outfilename_error);
		for(i=1;i<=(n+1);i++)
		{
			//write data to file
			//cout << errorContainer[i] << endl;
			ofile << setprecision(8) << errorContainer[i] << setw(15) << log10(h) << endl;
		}
		ofile.close();



}


/*
void error_analysis(char* inputfile)
{

}
*/























