#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "Jacobi_method.hpp"
#include "Testfunc.hpp"
//#include "catch.h"

using namespace std;
using namespace arma;



int main (int argc, char* argv[]) {
  //int n = 4;
  int n = atoi(argv[1]);
  double mass_e = 9.10938e-31;     // kg
  double hbar = 1.05457e-34;       // J s
  double omega = 1.0;             // s^-1     // 0.00038
  double omega_r = 0.01; // 0.5, 1.0, 5.0
  double alpha = pow(hbar*hbar/(mass_e*mass_e*omega*omega),0.25);
  cout << "Infinity = " << alpha << " omega="<<omega << endl;
  cout << n <<endl;
  int i; int j; int k, l;
  double h, rhomax, rho0, d, a;
  double infinity = 1.0e20;

  rhomax = 5.0; rho0 = 0.0;
  h = (rhomax)/((double) n);
  d = 2.0/(h*h);    // Diagonal elements
  a = -1.0/(h*h);   // Off Diagonal elements
  cout << "h= " << h << endl;
  double epsilon = 1.0e-10;

  vec eigval; mat eigvec; //vec lmbd=zeros(n);
  mat A = zeros(n,n);
  mat R = zeros(n,n);
  // Define the matrix A as a tridiagonal matrix:
  for (int i=0; i<n-1; i++){
    double rho = (i+1)*h;
    A(i,i) = d + omega_r*rho*rho + 1.0/rho ;
    A(i,i+1) = a;
    A(i+1,i) = a;
    R(i,i) = 1.0;
    //lmbd(i) = 4.0*i + 3.0;      // exact eigen values

  }
  //lmbd(n-1) = 4.0*(n-1) + 3.0;

  A(n-1,n-1) = d + omega_r*(n-1)*h*(n-1)*h + 1.0/((n-1)*h);
  A.print("A:");
  R(n-1,n-1) = 1.0;
  //R.print("R:");

  mat Arot = A;
  double max_number_iterations = (double)n*(double)n*(double)n;
  int iteration = 0;
  double maxoff = offdiag(A, k, l, n);
  cout << "max off diagonal="<< maxoff << endl;
  while (fabs(maxoff) > epsilon && (double) iteration < max_number_iterations){
    max:maxoff = offdiag(A, k, l, n);
    //mat Rot = Rotate(A, R, k, l, n);    // No void function for Rotate
    //A = Rot;
    Rotate(A, R, k, l, n);
    maxoff = offdiag(A, k, l, n);
    iteration++;
    //A = Rot.t()*A*Rot;


  }

  cout << "Number of iteration: " << iteration << endl;
  A.print("A:");
  //R.print("R:");

  // Test rotation:
  vec eigvals;
  eigvals = get_eigenvalues(A, n);
  eigvals.print("Eigenvalues");
  /*
  double lmbd;
  for (int i=0; i<n; i++){
    if (i < 10){
      lmbd = 4.0*i + 3.0;
      cout << "i=" << i <<" Computed eigenvalue= " <<eigvals(i) << " exact= " << lmbd << endl;
    }
  }
  mat V;
  V = get_eigenvectors(A, R, n);
  */
  //TestEigenvalues(eigvals, n);


  //V.print("Eigen vectors:");
  //(V.t()*V).print("Orthogonality?");
  return 0;
}
