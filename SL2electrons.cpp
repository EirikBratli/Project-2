/*
Solve Schrodinger equation for two electron i harmonic oscillator potential and
with repulsion
*/

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
  int n; int i; int j; int k, l; int w;

  string filename;
  if (argc <=1){
    cout << "Not enough arguments, need output filename and integration points N" << endl;
    exit(1);
  }
  else{
    filename = argv[1];
    n = atoi(argv[2]);
    //omega_r = atof(argv[3])
  }

  vec wr = zeros(4); //= [0.01, 0.5, 1.0, 5.0];
  wr(0) = 0.01; wr(1) = 0.5; wr(2) = 1.0; wr(3) = 5.0;
  double mass_e = 9.10938e-31;     // kg
  double hbar = 1.05457e-34;       // J s

  double h, rhomax, rho0, d, a;
  rhomax = 5.0; rho0 = 0.0;
  h = (rhomax)/((double) n);
  d = 2.0/(h*h);    // Diagonal elements
  a = -1.0/(h*h);   // Off Diagonal elements
  cout << "h= " << h << endl;
  double epsilon = 1.0e-10;

  // loop over different omega_r values
  for (int w=0; w<4; w++){
    // Define filename to output files:
    string outdata = filename;
    string argument = to_string(w+1);

    outdata.append("omega_r");
    outdata.append(argument);
    outdata.append(".txt");

    double omega_r = wr(w); // 0.5, 1.0, 5.0  // s^-1
    cout << "omega_r=" << omega_r << endl;

    double ww_r = omega_r*omega_r;    // Square of omega_r
    vec eigval; mat eigvec;
    mat A = zeros(n,n);
    mat R = zeros(n,n);
    // Define the matrix A as a tridiagonal matrix, and R:
    for (int i=0; i<n-1; i++){
      double rho = (i+1)*h;
      A(i,i) = d + ww_r*rho*rho + 1.0/rho;
      A(i,i+1) = a;
      A(i+1,i) = a;
      R(i,i) = 1.0;

    }

    A(n-1,n-1) = d + ww_r*(n-1)*h*(n-1)*h + 1.0/((n-1)*h);
    R(n-1,n-1) = 1.0;

    double max_number_iterations = (double)n*(double)n*(double)n;
    int iteration = 0;
    double maxoff = offdiag(A, k, l, n);
    cout << "max off diagonal="<< maxoff << endl;

    clock_t start, finish;
    start = clock();  // start timing
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
    cout << "Number of iteration: " << iteration << endl;
    cout << setprecision(10) << "Time used: " << time_used << " s at n=" << n << endl;
    //A.print("A:");
    //R.print("R:");

    // Test rotation:
    vec eigvals;
    eigvals = get_eigenvalues(A, n);
    //eigvals.print("Eigenvalues");

    mat V;
    V = get_eigenvectors(A, R, n);
    // test orthonormality:
    // test orthogonality:
    //(V.t()*V).print("orthogonality?")
    //(sum(V)).print("orthonormality?");

    ofile.open(outdata);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    int ex0 = 0; int ex1 = 1; int ex2 = 2;
    for(int i=0; i<n; i++){
      ofile << setw(15) << setprecision(8) << i*h << "  ";
      ofile << setw(15) << setprecision(8) << V(ex0, i) << "  ";
      ofile << "\n";
    }
    ofile.close();


  }
  return 0;
}
