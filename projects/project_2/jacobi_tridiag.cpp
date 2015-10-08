#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "time.h"
using namespace std;
using namespace arma;

// double maxNonDiag(mat& A, int* k, int* l, int n);
// void rotate(mat& A, mat& S, int k, int l, int n);
// void jacobi(int eps, mat& A, mat& S, int n);

ofstream ofile;
clock_t start, finish;
double maxNonDiag(mat& A, int* k, int* l, int n);
void rotate(mat& A, mat& S, int k, int l, int n);
void jacobi(int eps, mat& A, mat& S, int n, int m);
int unit_tests();

int main(int argc,char** argv)
{
	if( argc <= 4 )
	{
		cout << "Bad Usage of " << argv[0] << ": write also"
		 "1: outputfile eigenvalues 2:outputfile eigenvectors"
		 "3:matrix size 4:strength of oscillator potential, "
		 "on the same command line." << endl;
		exit(1);
	}
	else
	{
		char* eigval = argv[1];
		char* eigvec = argv[2];
		int n = atoi(argv[3]); // atoi = ascii to int, atof = ascii to float
		float omega_r = atof(argv[4]);

		int result = unit_tests();
		cout << result << endl;

		/*
		int m = n+2;
		double eps = 10e-8;
		double rho_min = 0.0;
		double rho_max = 4.0;
		double h = (rho_max-rho_min)/(m);
		double h_sq = h*h;
		double k1 = -1.0/h_sq;
		ofile.open("rho.txt");
		ofile << setprecision(5) << rho_min << setw(20) << rho_max << endl;
		ofile.close();

		//Initialize matrices and vectors
		mat A = eye<mat>(n,n);  //tridiag vector
		vec V(n);				//Potential = rho^2
		mat S = eye<mat>(n,n);	//Identity matrix

		for(int i=0 ; i<n ; i++){
			double rho = (i+1)*h;
			double rho_sq = rho*rho;
			//V(i) = rho_sq; 
			V(i) = (omega_r*omega_r*rho_sq) + (1.0/rho);

			A(i,i) = (2.0/h_sq)+V(i);
			if(i<n-1){

				A(i,i+1) = k1;
				A(i+1,i) = k1;
			}
		}

		//Checking eig.vals of A
		start = clock();
		mat test_eigval = eig_sym(A);
		finish = clock();
		double time_spent_arma = (finish - start)/double(CLOCKS_PER_SEC);
		cout << "When n = " << n << " the time spent with Armadillo is " << time_spent_arma << " sec." <<endl;
		//cout << eigval << endl;	

		//Run jacobi-method and measure the time spent
		start = clock();
		jacobi(eps,A,S,n,m);
		finish = clock();
		double time_spent = (finish - start)/double(CLOCKS_PER_SEC);
		cout << "When n = " << n << " the time spent is " << time_spent << " sec." <<endl;
		
		//Write eigenvalues and eigenvec's to files
		ofile.open(eigval);
		for(int i=0;i<n;i++)
		{
			ofile << setprecision(16) << A(i,i) << endl;
		}
		ofile.close();

		ofile.open(eigvec);
		ofile << n << endl;
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++)
			{

			ofile << setprecision(16) << S(i,j) << endl;

			}
		}
		ofile.close();
		*/
	}

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
	//cout << "i=" << *k << "and j=" << *l << endl;
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

int unit_tests()
{
	//Test of function maxNonDiag
	int n = 3;
	mat A(n,n);
	mat B(n,n);
	for(int i=0 ; i<n ; i++)
	{
		for(int j=0 ; j<n ; j++)
		{
			A(i,j) = i+j;
			B(i,j) = i-j;
		}
	}

	int k = 0;
	int l = 0;
	double max_A = maxNonDiag(A, &k, &l, n);
	double max_B = maxNonDiag(B, &k, &l, n);

	if(max_A != 3)
	{
		cout << "maxNonDiag() did not pass for positive element." << endl;
		exit(0);
	}

	if(max_B != 2)
	{
		cout << "maxNonDiag() did not pass for negative element." << endl;
		exit(0);
	}

	//Test of the jacobi-method
	int m = 2;
	int p = m+2;
	double eps = 10e-8;
	mat R = eye<mat>(m,m);	//Identity matrix
	mat C(m,m);
	for(int i=0 ; i<m ; i++)
	C(0,0) = 2;
	C(1,1) = 2;
	C(0,1) = 3;
	C(1,0) = 3;

	jacobi(eps, C, R, m, p);
	int l1 = C(0,0);
	int l2 = C(1,1);
	if(l1 == 5 || l2 == -1)//WHAAAT?????????????????????
	{
		int a = C(0,0) != 5;
		cout << "jacobi() did not pass." << endl;
		cout << a << C(0,0) << endl;
		exit(0);
	}

	//Test if eigvecs in matrix is orthogonal.
	mat S = eye<mat>(n,n);
	int q = n+2;
	jacobi(eps, A, S, n, q);

	int ind_1 = rand() % (n+1);
	int ind_2 = rand() % (n+1);
	cout << ind_1 << ind_2 << endl;
	cout << S << endl;




	return 0;
}













