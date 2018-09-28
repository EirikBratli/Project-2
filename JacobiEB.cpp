//#define CATCH_CONFIG_MAIN

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

ofstream ofile;     // create file for output


int main (int argc, char* argv[]) {
  int n;
  string filename;
  if (argc <=1){
    cout << "Not enough arguments, need output filename and integration points N" << endl;
    exit(1);
  }
  else{
    filename = argv[1];
    n = atoi(argv[2]);
  }
  // Update filename:
  string outdata = filename;
  outdata.append(".txt");

  int i; int j; int k, l;
  double h, Rmax, d, a;
  Rmax = 5.0;
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
        R(i,j) = 1.0;
      }
    }
  }
  A(n-1,n-1) = d;
  R(n-1,n-1) = 1.0;
  //R.print("R:");
  mat Arot = A;
  double max_number_iterations = (double)n*(double)n*(double)n;
  int iteration = 0;
  double maxoff = offdiag(A, k, l, n);
  cout << maxoff << endl;
  cout << "Max off diagonal= " << maxoff << endl;
  clock_t start, finish;
  start =clock();  // start timing
  while (fabs(maxoff) > epsilon && (double) iteration < max_number_iterations){
    max:maxoff = offdiag(A, k, l, n);
    //mat Rot = Rotate(A, R, k, l, n);    // No void function for Rotate
    //A = Rot;
    Rotate(A, R, k, l, n);
    maxoff = offdiag(A, k, l, n);
    iteration++;
    //A = Rot.t()*A*Rot;

  }
  finish =clock();   // end timing
  double time_used = (double)(finish - start)/(CLOCKS_PER_SEC );
  cout << "Max off diagonal= " << maxoff << " and epsilon= " << epsilon << endl;
  cout << "Number of iteration: " << iteration << endl;
  cout << setprecision(10) << "Time used: " << time_used << " s at n=" << n << endl;

  A.print("A:");
  //R.print("R:");

  // Test rotation:
  vec eigvals;
  eigvals = get_eigenvalues(A, n);
  eigvals.print("Eigenvalues");

  mat V;
  V = get_eigenvectors(A, R, n);
  V.print("Eigen vectors:");


  Orthogonality(R, V, n, epsilon);
  TestOffdiagonal(maxoff, epsilon);


  // Write to file, including: eigenvectors
  /*
  ofile.open(outdata);
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  int ex0 = 0; int ex1 = 1; int ex2 = 2;
  for(int i=0; i<n; i++){
    //for (int j=0; j<2; j++){
    ofile << setw(15) << setprecision(8) << V(i,ex0) << "  ";
    ofile << setw(15) << setprecision(8) << V(i,ex1) << "  ";
    ofile << setw(15) << setprecision(8) << V(i,ex2) << "  ";
    //}
    ofile << "\n";
  }
  ofile.close();
  */
  vec Eval;
  eig_sym(Eval, eigvec, A);    // Armadillo function

  // analytic Eigenvalues;
  double pi = acos(-1.0);
  double aa = 2*a; double fac2 = pi/(n+1);
  for (int i=0; i<n; i++){
    double EigValExcact = d + aa*cos((i+1)*fac2);
    //cout << EigValExcact << " & " << eigvals(i)<< " & " << Eval(i) << endl;
  }

  return 0;
}
