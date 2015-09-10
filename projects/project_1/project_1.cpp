#include <iostream>
//#include<armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
using namespace std;
//using namespace arma;

ofstream ofile;
double f(double x);
double u_analytic(double x);

int main(int args, char** argv)
{

		if( args <= 1 )
	{
		cout << "Bad Usage " << argv[0] <<
		": read also output file on same line." << endl;
		exit(1);
	}
	else
	{
		char* outfilename = argv[1];
		int n, i;

		cout << "Enter an integer n for the nxn-matrix:" << endl;
		cin >> n;

		double h = (1.0/(n+1));
		double h_sq = h*h;

		float* a = new float[n];
		float* b = new float[n];
		float* c = new float[n];
		double* x = new double[n+2];
		double* y = new double[n];
		double* u = new double[n+2];
		double* u_ana = new double[n+2];

		for(i=0;i<=(n-1);i++)
		{
			//fill vectors
			x[i+1] = 0.0 + (i+1)*h;
			a[i] = 1.0;
			b[i] = -2.0;
			c[i] = 1.0;
			y[i] = -h_sq*f(x[i+1]);
		}

		//boundary conditions on x
		x[0] = 0;
		x[n+1] = 1;

		for(i=1;i<=(n-1);i++)
		{
			//forward substitution
			double factor = (a[i-1]/b[i-1]);
			b[i] = b[i] - (factor*c[i-1]);
			y[i] = y[i] - (factor*y[i-1]);
			a[i-1] = 0;
		}

		//boundary conditions on u(x)
		u[n+1] = 0;
		u[0] = 0;

		for(i=n;i>=1;i--)
		{
			//fill u(x) by backward substitution
			u[i] = (1/b[i-1])*(y[i-1]-(c[i-1]*u[i+1]));
		}

		for(i=0;i<=(n+1);i++)
		{
			//fill vector with analytical values of u(x)
			u_ana[i] = u_analytic(x[i]);
		}

		ofile.open(outfilename);

		for(i=0;i<=(n+1);i++)
		{
			//write data to file
			ofile << setw(15) << setprecision(8) << x[i] << setw(15) << u[i] << endl;
		}

		ofile.close();
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








