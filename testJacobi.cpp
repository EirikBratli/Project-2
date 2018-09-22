//#define CATCH_CONFIG_MAIN

#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "Jacobi_method.hpp"
//#include "catch.h"

using namespace std;
using namespace arma;


int main (int argc, char* argv[]) {
  //int n = 4;
  int n = atoi(argv[1]);
  int i; int j; int k, l;
  double h, Rmax, d, a;
  Rmax = 1.0;
  h = Rmax/n;
  d = 2.0/(h*h); // Diagonal elements
  a = -1.0/(h*h);   // Off Diagonal elements

  double epsilon = 1.0e-10;
  vec eigval; mat eigvec;
  mat A = zeros(n,n);
  mat R = zeros(n,n);
  // Define the matrix A as a tridiagonal matrix:
  for (int i=0; i<n; i++){
    for (int j=0; j<n-1; j++){
      if (i==j){
        A(i,i) = d;
        A(i,j+1) = a;
        A(i+1,j) = a;
      }
    }
  }
  A(n-1,n-1) = d;

  mat Arot = A;
  double max_number_iterations = (double)n*(double)n*(double)n;
  int iteration = 0;
  double maxoff = offdiag(A, k, l, n);
  cout << maxoff << endl;
  while (fabs(maxoff) > epsilon && (double) iteration < max_number_iterations){
    max:maxoff = offdiag(A, k, l, n);
    //mat Rot = Rotate(A, R, k, l, n);    // No void function for Rotate
    //A = Rot;
    Rotate(A, R, k, l, n);
    maxoff = offdiag(A, k, l, n);
    iteration++;
    //Arot = Rot.t()*Arot*Rot;


  }
  cout << "Number of iteration: " << iteration << endl;
  A.print("A:");
  //R.print("R:");

  // test of the maxoffdiag function, and see the difference from armadillo function.
  /*mat C = A;
  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++){
      if (i<j){
        C(i,j) = maxoff;
      }
      if (i > j){
        C(i,j) = maxoff;
      }
    }
  }

  */
  // Test rotation:


  //mat B = diagmat(A);
  //B.print("B:");
  eig_sym(eigval, eigvec, A);
  eigvec.print("Eigenvectors");
  eigval.print("Eigenvalues");


  // analytic Eigenvalues;
  double pi = acos(-1.0);
  double aa = 2*a; double fac2 = pi/(n+1);
  for (int i=0; i<n; i++){
    double EigValExcact = d + aa*cos((i+1)*fac2);
    //cout << EigValExcact << endl;
  }
  // Test orthogonality:
  mat V = eigvec.t()*eigvec;
  V.print("V:");
  
  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++){
      if (i==j){
        if (V(i,j) >= 1.0-epsilon){
          cout << "OK, Diagonal elements = 1, i,j= " <<i <<j<< endl;
        }
        else{
          cout <<"i,j="<<i<<", "<<j<<" "<< setw(10)<<setprecision(15)<<V(i,j) <<" Too big difference from 1"<< endl;
        }
      }
      if (i != j){
        if (V(i,j) < epsilon){
          cout << "Ok, off-diagonal elements is 0, i,j= "<<i<<j << endl;
        }
        else{
          cout <<"i,j="<<i<<", "<<j<<" "<< setw(10)<<setprecision(15)<<V(i,j) <<" Too big difference from 0" << endl;
        }
      }
    }
  }



  return 0;
}


/*
double offdiag(mat A, int &k, int &l, int n){
  double max = 0.0;
  //cout << n << " " << max << endl;
  for (int i=0; i<n; i++){
    //cout << i << " " << max << endl;
    for (int j=i+1; j<n; j++){
      double aij = fabs(A(i,j));
      //cout << i << " second loop" << j << endl;
      if (aij >= max){
        max = aij;
        k = i;
        l = j;
      }
    }
  }
  //cout << "max=" << max << " k=" << k << " l=" << l << endl;
  return max;
}

void Rotate( mat& A, mat& R, int k, int l, int n){
//mat Rotate( mat& A, mat& R, int k, int l, int n){
  //A.print("A_before:");
  double s, c;
  //cout << "A_kl=" << A(k,l) << endl;
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
  //cout << "c="<< c <<" s=" << s << "k,l"<< k<<l<< endl;

  double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
  a_kk = A(k,k);
  a_ll = A(l,l);
  A(k,k) = a_kk*c*c - 2*A(k,l)*c*s + a_ll*s*s;
  A(l,l) = a_ll*c*c + 2*A(k,l)*c*s + a_kk*s*s;
  A(k,l) = 0.0;
  A(l,k) = 0.0;
  //A.print("A:");

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
  }
  //A.print("A_after:");
  //return A; // Not orthogonal!!!
}
*/