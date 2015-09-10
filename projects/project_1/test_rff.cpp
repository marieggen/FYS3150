#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;



int main()
{
	const char* filename = "data_test.txt";
	ifstream ifile(filename);

	if(!ifile) 
	{
    	cout << endl << "Failed to open file " << filename;
    	return 1;
  	}

  	cout << "Reached first step" << endl;
  	
	int i;
	int lines = 12;
	double* c1 = new double[lines];
	double* c2 = new double[lines];
	double* c3 = new double[lines];
	ifile >> lines;

	cout << "Reached second step" << endl;

	for(i=1;i<=lines;i++)
	{
		cout << ifile << endl;
		ifile >> c1[i] >> c2[i] >> c3[i];
	}

	for(i=1;i<=lines;i++)
	{
		cout << c1[i] << "   " << c2[i] << "   " << c3[i] << endl;
	}

	cout << "Reached last step" << endl;

	return 0;
}