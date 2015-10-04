#include <iostream>
#include <armadillo>
#include <cmath>
using namespace std;
using namespace arma;

int main(){
	int n = 4;
	double eps = 10e-8;
	double rho_min = 0.0;
	double rho_max = 1.0;
	double h = (rho_max-rho_min)/(n);
	double h_sq = h*h;

	for(int i=0;i<=n;i++){

		double rho = rho_min + i*h;
		cout << rho << endl;
	}


	return 0;
}
