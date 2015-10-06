#include <iostream>
#include <armadillo>
#include <cmath>
using namespace std;
using namespace arma;

// double maxNonDiag(mat& A, int* k, int* l, int n);
// void rotate(mat& A, mat& S, int k, int l, int n);
// void jacobi(int eps, mat& A, mat& S, int n);

double maxNonDiag(mat& A, int* k, int* l, int n);
void rotate(mat& A, mat& S, int k, int l, int n);
void jacobi(int eps, mat& A, mat& S, int n, int m);

int main(int argc,char** argv)
{
	/*
	int n = 10;
	int m = n+2;
	double eps = 10e-8;
	double rho_min = 0.0;
	double rho_max = 5.0;
	double h = (rho_max-rho_min)/(m);
	double h_sq = h*h;
	double k1 = -1.0/h_sq;

	//Initialize matrices and vectors
	mat A = eye<mat>(n,n);  //tridiag vector
	vec V(n);				//Potential = rho^2
	mat S = eye<mat>(n,n);	//Identity matrix

	for(int i=0 ; i<n ; i++){
		V(i) = (i+1)*h*(i+1)*h; 

		A(i,i) = (2.0/h_sq)+V(i);
		if(i<n-1){

			A(i,i+1) = k1;
			A(i+1,i) = k1;
		}
	}
	//cout << A << endl;
	//Checking eig.vals of A
	mat eigval = eig_sym(A);
	cout << eigval(0) << " , " << eigval(1) << " , " << eigval(2) << endl;	

	//Run jacobi-method
	jacobi(eps,A,S,n,m);
	cout << A(0,0) << " , " << A(1,1) << " , " << A(2,2) << endl;
*/
	
	int n = 4;
	int m = n+2;
	double eps = 10e-8;
	mat S = eye<mat>(n,n);

	mat A(n,n);
	/*
	for(int i=0 ; i<n ; i++){
		for(int j=0 ; j<n ; j++){
			A(i,j) = (double) j+2* (double) i;
		}
	}*/

	for(int i=0 ; i<n ; i++){
		for(int j=0 ; j<n ; j++){
			if(i==j){
				A(i,j) = 2;
			}
			else{
				A(i,j) = 3;
			}
		}
	}

	cout << A << endl;

	mat eigval = eig_sym(A);
	cout << eigval << endl;

	jacobi(eps,A,S,n,m);
	cout << A(0,0) << " , " << A(1,1) << " , " << A(2,2) << "," << A(3,3) << endl;
	

	return 0;
}



void jacobi(int eps, mat& A, mat& S, int n, int m)
{
	int iterations = 0, k=0, l=0;
	double maxIterations = (double) m * (double) m * (double) m;
	double maxElement = 1.0;

	while(fabs(maxElement) > eps && (double) iterations < maxIterations){
		maxElement = maxNonDiag(A, &k, &l, n);
		//cout << A << endl;
		//cout << "max element in matrix: "<< maxElement << endl;
		rotate(A,S,k,l,n);
		iterations++;
	}

	//cout << "number of rotations: "<< iterations-1 << endl;
	//cout << "max iterations: " << maxIterations << endl;

	return;
}

double maxNonDiag(mat& A, int* k, int* l, int n)
{
	double max=0.0;
	int a = 0;

	for(int i=0 ; i<(n-1) ; i++){
		a++;
		for(int j=a ; j<n ; j++){
			if(fabs(A(i,j))>max){

				max = fabs(A(i,j));
				*k = i;
				*l = j;

			}
		}
	}
	cout << "i=" << *k << "and j=" << *l << endl;
	return max;
}


// Function to find the values of cos and sin
void rotate(mat& A, mat& S, int k, int l, int n) 
{
	double s, c, t;
	if ( A(k,l) != 0.0 ) {
		
		double tau = (A(l,l) - A(k,k))/(2*A(k,l)); 
		//cout << "tau: " << tau << endl;

		if ( tau > 0 ) {
    		t = 1.0/(tau + sqrt(1.0 + tau*tau));
    		//t = -tau + sqrt(1.0 + tau*tau);
    		//cout << "(tau > 0) tan: " << t << endl; 
		} 
		else {
			t = -1.0/( -tau + sqrt(1.0 + tau*tau));
			//t = -tau - sqrt(1.0 + tau*tau);
			//cout << "(tau <= 0) tan: " << t << endl;

		}
		c = 1.0/sqrt(1+(t*t));
		s = c*t;
	} 
	
	else {c = 1.0; s = 0.0;}

	//cout << "cos and sin: "<< c << " , " << s << endl;
	double a_kk, a_ll, a_ik, a_il, s_ik, s_il;
	a_kk = A(k,k);
	a_ll = A(l,l);

	//Change the matrix elements with indices k and l 
	A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll; 
	A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll; 
	A(k,l) = 0.0; 
	A(l,k) = 0.0;

	//Change the remaining elements
	for (int i = 0 ; i < (n) ; i++){
		if (i != k && i != l){

			a_ik = A(i,k);
			a_il = A(i,l);
			A(i,k) = c*a_ik - s*a_il;  
			A(i,l) = c*a_il + s*a_ik; 
			A(k,i) = A(i,k);
			A(l,i) = A(i,l);

		}
		// Finally, we compute the new eigenvectors
   		s_ik = S(i,k);
   		s_il = S(i,l);
   		S(i,k) = c*s_ik - s*s_il;
   		S(i,l) = c*s_il + s*s_ik;
     		
	}
	//cout << "The number of rotations is: " << numRotations << endl;
	return;
}

