#ifndef JACOBI_METHOD_H
#define JACOBI_METHOD_H

#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <catch2.h>

using namespace std;
using namespace arma;

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

//caclulate eigenvectors and eigenvalues
vector<double> get_eigenvals(mat A,int n){
    vector<double> eigen;
    for(int i=0;i<n;i++){
        eigen.push_back(A(i,i));        //add element to vector
    }
    sort (eigen.begin(), eigen.begin()+n);
    return eigen;
}

mat get_eigenvecs(mat a, mat v, int n){
    vector<double>eigenvals=get_eigenvals(a,n);
    mat vecs(n,n);
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(a(j,j)==eigenvals[i]){
                for(int k=0;k<n;k++){
                      vecs(i,k)=v(k,j);
                }
             }
         }
    }
    return vecs;
}

//test eigenvectors and orthogonality:
double tolerance = 1e-15;
double orthogonality(mat ortho, mat oto, int n){
for (int i=0; i<n; i++){
    for(int j=0; j<n; j++){
        if (i==j){
            if(ortho(i, j) != Approx(0.9996).epsilon(0.001) && oto(i, j) != Approx(1e00)){
            cout << "eigvec " << i << " - Failed" << endl;
            return 1;}
        }
        else{
            if(fabs(ortho(i, j)) != Approx(1e-16).margin(tolerance) && oto(i, j) != Approx(1e-17).margin(tolerance)){
            cout << "eigvec " << i << " - Failed" << endl;
            return 1;}
            }
        }
   cout << "eigvec " << i  << " - Passed"<< endl;
}
}
/*
TEST_CASE("Testing eigenvector orthogonality"){
  int n;
  mat ortho; mat eigenvec;
  ortho = eigenvec.t()*eigenvec;
  for(int i=0; i<n; i++){
      for(int j=0; j<n; j++){
          if(i==j){
              REQUIRE(ortho(i, j)==Approx(1e+00).epsilon(0.001));
          }
          else{
              REQUIRE(ortho(i, j)==Approx(1e-16).margin(1e-14));
          }
      }
  }
}
*/
#endif // JACOBI_METHOD_H
