#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace arma;
//double offdiag(double, int *, int *, int );

double offdiag(mat A, int *k, int *l, int n){
  double max = 0.0;
  for (int i=0; i<n; i++){
    for (int j=i+1; j<n; j++){
      if (fabs(A(i,j)) > max){
        max = fabs(A(i,j));
        *l = i;
        *k = j;
      }
    }
  }
  return max;
}

void Rotate( mat A, mat R, int k, int l, int n){
  double s, c;
  if (A(k,l) != 0.0){
    double t, tau;
    tau = (A(l,l) - A(k,k))/(A(k,l)*2.0);

    if (tau >= 0){
      t = 1.0/(tau + sqrt(1+tau));
    }
    else {
      t = -1.0/(-tau + sqrt(1+tau));
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
      A(i,l) = a_il*c - a_ik*s;
      A(l,i) = A(i,l);
    }
    // Calculate the eigen vectors
    r_ik = R(i,k);
    r_il = R(i,l);

    R(i,k) = r_ik*c + r_il*s;
    R(i,l) = r_il*c + r_ik*s;
  }

  return;
}


int main (int n, char* argv[]) {
  n = atoi(argv[1]);
  int i; int j; int k, l;
  double epsilon = 1.0e-10;
  vec eigval; mat eigvec;
  mat A = zeros(n,n);
  mat R = zeros(n,n);
  for (int i=0; i<n; i++){
    for (int j=0; j<n-1; j++){
      if (i==j){
        A(i,i) = 2.0;
        A(i,j+1) = -1.0;
        A(i+1,j) = -1.0;
      }
    }
  }
  A(n-1,n-1) = 2.0;

  double max_number_iterations = (double)n*(double)n*(double)n;
  int iteration = 0;
  double maxoff = offdiag(A, &k, &l, n); //A.print("A:");
  while (fabs(maxoff) > epsilon && (double) iteration < max_number_iterations){
    max:maxoff = offdiag(A, &k, &l, n);
    Rotate(A, R, k, l, n);
    iteration++;
  }
  cout << "Number of iteration: " << iteration << endl;

  // test of the maxoffdiag function, and see the difference from armadillo function.
  mat C = A;
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

  mat B = diagmat(A);
  eig_sym(eigval, eigvec, B);
  //eigval.print("Eigenvalues");
  B.print("Diag:");
  C.print("C:");
  return 0;
}
