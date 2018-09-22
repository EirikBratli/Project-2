#ifndef JACOBI_METHOD_H
#define JACOBI_METHOD_H

#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace arma;

double offdiag(mat A, int &k, int &l, int n){
    double max=0.0;
    for (int i=0; i<n; i++){
        for (int j=i+1; j<n; j++){      //only look at off-diagonal elements
            if(fabs(A(i, j)) > max){
                l = i;
                k = j;
                max = fabs(A(i, j));
                cout << max <<endl;
                cout << "k = " << k << endl;
                cout << "l = " << l << endl;
            }
        }
    }
return max;
}

void Rotate( mat &A, mat &R, int k, int l, int n){
  double s, c;
  if (A(k,l) != 0.0){
    double t, tau;
    tau = (A(l,l) - A(k,k))/(double (A(k,l)*2.0));

    if (tau >= 0){
      t = 1.0/(tau + sqrt(1+tau*tau));
    }
    else {
      t = -1.0/(-tau + sqrt(1+tau*tau));
    }

    c = 1.0/(sqrt(1+t*t));
    s = c*t;
  }
  else{
    c = 1.0;
    s = 0.0;
  }
//approximate off-diagonal elements A(k, l), A(l, k) to zero
  double a_kk, a_ll, a_ik, a_il, r_ik, r_il, theta;
  a_kk = A(k,k);
  a_ll = A(l,l);
  A(k,k) = a_kk*c*c - 2.0*A(k,l)*c*s + a_ll*s*s;
  A(l,l) = a_ll*c*c + 2.0*A(k,l)*c*s + a_kk*s*s;
  A(k,l) = 0.0;
  A(l,k) = 0.0;
  for (int i=0; i<n-1; i++){
    if (i != k && i != l){
      a_ik = A(i,k);
      a_il = A(i,l);
      A(i,k) = a_ik*c - a_il*s;
      A(k,i) = A(i,k);
      A(i,l) = a_il*c - a_ik*s;
      A(l,i) = A(i,l);
    }
    // Calculate the eigenvectors
    r_ik = R(i,k);
    r_il = R(i,l);

    R(i,k) = r_ik*c + r_il*s;
    R(i,l) = r_il*c + r_ik*s;
  }
  theta = asin(s);
//  cout << "theta = " << theta << endl;
}


#endif // JACOBI_METHOD_H
