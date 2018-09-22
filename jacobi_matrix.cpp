#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <jacobi_method.h>

using namespace std;
using namespace arma;


int main (int argc, char* argv[]) {
  int i, j, k, l;
  double h; double epsilon = 1.0e-8;
  //n = atoi(argv[1]);
  int n= 4;
  h = 1/((double) n);
  vec eigval; mat eigvec;
  mat A = zeros(n,n);// A.randu(n, n);
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
A = 1/((double) pow((double)h, 2))*A;


//norm of matrix:
double norm=0;
for(int m=0; m<n; m++){
    norm += pow(fabs(A(m, m)), 2);
}
cout << "norm" << norm<< endl;
//test of off-diagonal elements
  double max_number_iterations = pow(n, 3);
  int iteration = 0;
  double maxoff = offdiag(A, k, l, n); //A.print("A:");

int iterations =0;
  while (fabs(maxoff) > epsilon && (double) iteration < max_number_iterations){
    maxoff = max(offdiag(A, k, l, n));
    //mat Rot = Rotate(A, R, k, l, n);
    Rotate(A, R, k, l, n);
//    A = Rot;
    iteration++;
  }
  cout << "Number of iterations: " << iteration << endl;

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

A.print("A");
/*mat B;
  mat B = diagmat(A);
  eig_sym(eigval, eigvec, B);     //eigen decomposition of matrix B
  eigval.print("Eigenvalues");
  B.print("Diag:");
  C.print("C:");

  return 0;
*/}
