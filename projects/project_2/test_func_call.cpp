#include <iostream>
#include <armadillo>
#include <cmath>
using namespace std;
using namespace arma;

void trans(mat A, mat& C);

int main(int argc,char** argv)
{
	//Initialize matrices
	int n = 3,i,j;
	mat A(n,n);
	for(i=0 ; i<n ; i++){
		for(j=0 ; j<n ; j++){
			A(i,j) = (double) j+ (double) i;
		}
	}
	mat C(n,n);	
	cout << C << endl;
	trans(A,C);
	cout << A << endl;
	cout << C << endl;

	return 0;
}

void trans(mat A, mat& C)
{
	C = strans(A);
	return;
}

