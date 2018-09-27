//#define CATCH_CONFIG_MAIN
//#include "catch.h"

#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>

using namespace std;
using namespace arma;



vec get_eigenvalues(mat A, int n){
  vec eigen(n);
  for (int i=0; i<n; i++){
    eigen(i) = A(i,i);
  }
  eigen(n-1) = A(n-1,n-1);
  sort(eigen.begin(), eigen.begin()+n);
  return eigen;
}

mat get_eigenvectors(mat A, mat V, int n){
  vec eigenvals=get_eigenvalues(A,n);
  mat vecs(n,n);
  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++){
      if (A(j,j) == eigenvals(i)){
        for (int k=0; k<n; k++){
          vecs(i,k) = V(k,j);
        }
      }
    }
  }
  return vecs;
}


double offdiag(mat A, int &k, int &l, int n){
  double max = 0.0;
  for (int i=0; i<n; i++){
    for (int j=i+1; j<n; j++){
      double aij = fabs(A(i,j));
      if (aij >= max){
        max = aij;
        k = i;
        l = j;
      }
    }
  }
  return max;
}

void Rotate( mat& A, mat& R, int k, int l, int n){
//mat Rotate( mat& A, mat& R, int k, int l, int n){
  double s, c;
  if (A(k,l) != 0.0){
    double t, tau;
    tau = (A(l,l) - A(k,k))/(A(k,l)*2.0);

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

  double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
  a_kk = A(k,k);
  a_ll = A(l,l);
  A(k,k) = a_kk*c*c - 2*A(k,l)*c*s + a_ll*s*s;
  A(l,l) = a_ll*c*c + 2*A(k,l)*c*s + a_kk*s*s;
  A(k,l) = 0.0;
  A(l,k) = 0.0;


  for (int i=0; i<n-1; i++){
    if (i != k && i != l) {
      a_ik = A(i,k);
      a_il = A(i,l);
      A(i,k) = a_ik*c - a_il*s;
      A(k,i) = A(i,k);
      A(i,l) = a_il*c + a_ik*s;
      A(l,i) = A(i,l);
    }
    // Calculate the eigen vectors
    r_ik = R(i,k);
    r_il = R(i,l);

    R(i,k) = r_ik*c - r_il*s;
    R(i,l) = r_il*c + r_ik*s;
    if (i==n-2){
      r_ik = R(i+1,k);
      r_il = R(i+1,l);

      R(i+1,k) = r_ik*c - r_il*s;
      R(i+1,l) = r_il*c + r_ik*s;
    }
  }

  //return A; // if not void function
}
