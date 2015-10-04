#include <iostream>
#include <armadillo>
#include <cmath>
using namespace std;
using namespace arma;

double maxNonDiag(mat& A, int* k, int* l, int n);
void rotate(mat& A, mat& S, int k, int l, int n);
void jacobi(int eps, mat& A, mat& S, int n);


int main(int argc,char** argv)
{
	int n = 10,i,j;
	double eps = 10e-8;
	//Initialize mat&rices
	mat A(n,n);
	for(i=0 ; i<n ; i++){
		for(j=0 ; j<n ; j++){
			A(i,j) = j+i;
		}
	}
	mat S = eye<mat>(n,n);	
	mat S_T = strans(S);

	jacobi(eps,A,S,n);
	cout << A << endl;

	return 0;
}


void jacobi(int eps, mat& A, mat& S, int n)
{
	int iterations = 0, k=0, l=0;
	double maxIterations = (double) n * (double) n * (double) n;
	double maxElement = maxNonDiag(A, &k, &l, n);

	while(fabs(maxElement) > eps && (double) iterations < maxIterations){
		maxElement = maxNonDiag(A, &k, &l, n);
		rotate(A,S,k,l,n);
		iterations++;
		//cout << "max element in matrix: "<< maxElement << endl;
		//cout << A << endl;
	}

	//cout << "number of iterations: "<< iterations << endl;
	//cout << "max iterations: " << maxIterations << endl;

	return;
}

double maxNonDiag(mat& A, int* k, int* l, int n)
{
	int i,j;
	double max=0.0;
	for(i=0 ; i<n ; i++){
		for(j=0 ; j<n ; j++){
			if(i!=j){
				if(fabs(A(i,j))>max){
					max = fabs(A(i,j));
					*k = i;
					*l = j;
				}
			}
		}
	}
	//cout << "j=" << *k << "and i=" << *l << endl;
	return max;
}


// Function to find the values of cos and sin
void rotate(mat& A, mat& S, int k, int l, int n) 
{
	int i,j,numRotations;
	double s, c, t, tau;
	if ( A(k,l) != 0.0 ) {
		numRotations++;
		tau = (A(l,l) - A(k,k)/(2*A(k,l))); 

		if ( tau > 0 ) {
    		t = 1.0/(tau + sqrt(1.0 + tau*tau));
		} 
		else {
			t = -1.0/( -tau + sqrt(1.0 + tau*tau));
		}
		c = 1/sqrt(1+t*t);
		s = c*t;
	} 
	
	else {c = 1.0; s = 0.0;}


	double a_kk, a_ll, a_ik, a_il, s_ik, s_il;
	a_kk = A(k,k);
	a_ll = A(l,l);

	//Change the mattrix elements with indices k and l 
	A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll; 
	A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll; 
	A(k,l) = 0.0; 
	A(l,k) = 0.0;

	//Change the remaining elements
	for (i = 0; i < n; i++){
		if (i != k && i != l){
			a_ik = A(i,k);
			a_il = A(i,l);
			A(i,k) = c*a_ik - s*a_il; 
			A(k,i) = A(i,k); 
			A(i,l) = c*a_il + s*a_ik; 
			A(l,i) = A(i,l);
		}
		// Finally, we compute the new eigenvectors
   		s_ik = S(i,k);
   		s_il = S(i,l);
   		S(i,k) = c*s_ik - s*s_il;
   		S(i,l) = c*s_il + s*s_ik;
	}
	cout << "The number of rotations is: " << numRotations << endl;
	return;
}











