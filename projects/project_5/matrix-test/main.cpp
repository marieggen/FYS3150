#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;

int main()
{

    int n = 200;
    double delta_x = 0.02;
    double delta_t = 0.2;
    double a = delta_t/(delta_x*delta_x);

    mat I(n,n,fill::eye);
    mat B(n,n,fill::zeros);
    for(int i = 0 ; i < n ; i++){
        B(i,i) = 2;
        if(i<(n-1)){
            B(i,i+1) = -1;
            B(i+1,i) = -1;
        }
    }



    mat A_bkwrd(n,n);
    A_bkwrd = I + a*B;
    mat A_bkwrd_inv = inv(A_bkwrd);
    vec eigval_bkwrd = eig_sym(A_bkwrd_inv);
    //cout << eigval_bkrwd << endl;

    mat A_forw(n,n);
    A_forw = I - a*B;
    vec eigval_forw = eig_sym(A_forw);
    //cout << eigval_forw << endl;

    mat A_CrN(n,n);
    mat A_CrN1(n,n);
    mat A_CrN2(n,n);
    A_CrN1 = I + a*B;
    A_CrN2 = I - a*B;
    mat A_CrN1_inv = inv(A_CrN1);
    mat A_CrN2_inv = inv(A_CrN2);
    A_CrN = A_CrN1_inv*A_CrN2_inv;

    vec eigval_CrN = eig_sym(A_CrN);
    cout << eigval_CrN << endl;

    return 0;
}

