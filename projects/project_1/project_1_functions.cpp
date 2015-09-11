#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <vector>
#include "lib.cpp"
using namespace std;

ofstream ofile;
clock_t start, finish;
double f(double x);
double u_analytic(double x);
void solver_double_derivative(int n, char* outfilename);
void error_analysis(char*inFilename,char*outFilename, int n);

// typedef vector<vector<double> > mat;
// typedef vector<double> vec;

int main(int args, char** argv)
{

		if( args <= 4 )
	{
		cout << "Bad Usage " << argv[0] <<
		": read also output file on same line." << endl;
		exit(1);
	}
	else
	{
		char* dataSolver = argv[1];
		char*dataError = argv[2];
		char*time_lu = argv[3];
		int n = atoi(argv[4]); // atoi = ascii to int, atof = ascii to float
		double h = (1.0/(n+1.0));
		double h_sq = h*h;
		int i,j;

		solver_double_derivative(n, dataSolver);
		error_analysis(dataSolver, dataError, n);

		//time taking of LU-decomposision. Functions from lib.cpp.

		
		double* x = new double[n];
		double* y = new double[n];
		// mat A(n, vec(n,0));
		
		double** A = new double*[n];
		for (i = 0; i <= (n-1); i++) {
			A[i] = new double[n];
		}

		for(i=0;i<n;i++)
		{	
			//fill vectors and matrix
			x[i] = (i+1)*h;
			y[i] = h_sq*f(x[i]);

			for(j=0;j<n;j++)
			{
				A[i][j] = 0;

				if(i == j){
					A[i][j] = 2;
				}
				if(i == (j-1)){
					A[i][j] = -1;
				}
				if(i == (j+1)){ 
					A[i][j] = -1;
				}
			}
		}

		start = clock();
		int *indx = new int[n];
		double d;
		ludcmp(A,n,indx,&d);
		lubksb(A,n,indx,y);
		finish = clock();	
		double time_spent = (finish - start)/double(CLOCKS_PER_SEC);
		ofile.open(time_lu);
		for(i=0;i<n;i++)
		{
			ofile << setprecision(16) << y[i] << setw(25) << log10(h) << setw(25) << time_spent << endl;
		}
		ofile.close();

		delete y;
		delete x;
		delete []A;
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

void solver_double_derivative(int n, char*outfilename)
{
	//function that writes solution data of u(x) and x to a file

	int i;
	double h = (1.0/(n+1.0));
	double h_sq = h*h;

	double* a = new double[n];
	double* b = new double[n];
	double* c = new double[n];
	double* x = new double[n+2];
	double* y = new double[n];
	double* u = new double[n+2];
	double* u_ana = new double[n+2];
	double* errorContainer = new double[n+2];

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
			ofile << setprecision(16) << x[i] << setw(25) << u[i] << setw(25) << u_ana[i] << setw(25) << time_spent << endl;
		}
		ofile.close();
}

void error_analysis(char*inFilename,char*outFilename, int n)
{	
	//Estimates relative error of a numerical solution and writes 
	//the error data to file.

	int i;
	int count = 0;
	double h = (1.0/(n+1.0));
	double x,numerical,analytical,time_spent;
	vector<double> errorContainer;
	vector<double> xValues;
  	vector<double> numericalValues;
 	vector<double> analyticalValues;


	ifstream inFile(inFilename);
 	// Make sure the file stream is good
  	if(!inFile) 
  	{
    cout << endl << "Failed to open file " << inFilename;
    }
 
    
    while (inFile >> x >> numerical >> analytical >> time_spent)
    {
    	xValues.push_back(x);
    	numericalValues.push_back(numerical);
    	analyticalValues.push_back(analytical);
    	if(analytical != 0)
		{
    		double absError = abs((numerical-analytical)/analytical);
    		errorContainer.push_back(log10(absError));
    		count++;
    	}
  	}

	ofile.open(outFilename);
	for(i=0;i<=count;i++)
	{
		//write data to file
		ofile << setprecision(8) << errorContainer[i] << setw(15) << log10(h) << endl;
	}
	ofile.close();

}





 

























