/*
Solve Schrodinger equation for one electron i harmonic oscillator potential
Can be used for solving for two electrons with no repulsion.
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

/*
To do:
Try to get more precise eigenvalues by more integration points and/or adjust rhomax higher.
*/

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
  }
  vec wr = zeros(4); //= [0.01, 0.5, 1.0, 5.0];
  wr(0) = 0.01; wr(1) = 0.5; wr(2) = 1.0; wr(3) = 5.0;
  double mass_e = 9.10938e-31;     // kg
  double hbar = 1.05457e-34;       // J s
  //double omega = 1.0;             // s^-1     // 0.00038
  //double alpha = pow(hbar*hbar/(mass_e*mass_e*omega*omega),0.25);
  double h, rhomax, rho0, d, a;
  double infinity = 1.0e20;

  rhomax = 5.0; rho0 = 0.0;
  h = (rhomax)/((double) n);
  d = 2.0/(h*h);    // Diagonal elements
  a = -1.0/(h*h);   // Off Diagonal elements
  cout << "h= " << h << endl;
  double epsilon = 1.0e-10;

  for (int w=0; w<4; w++){

    string outdata = filename;
    string argument = to_string(w+1);
    outdata.append(argument);
    outdata.append(".txt");


    double omega_r = wr(w); // 0.5, 1.0, 5.0  // s^-1
    cout << "omega_r=" << omega_r << endl;
    double wwr = omega_r*omega_r;    // Square of omega_r

    vec eigval; mat eigvec; //vec lmbd=zeros(n);
    mat A = zeros(n,n);
    mat R = zeros(n,n);
    // Define the matrix A as a tridiagonal matrix:
    for (int i=0; i<n-1; i++){
      double rho = (i+1)*h;
      A(i,i) = d + rho*rho*wwr;
      A(i,i+1) = a;
      A(i+1,i) = a;
      R(i,i) = 1.0;
      //lmbd(i) = 4.0*i + 3.0;      // exact eigen values

    }
    //lmbd(n-1) = 4.0*(n-1) + 3.0;

    A(n-1,n-1) = d + (n-1)*h*(n-1)*h*wwr;
    //A.print("A:");
    R(n-1,n-1) = 1.0;
    //R.print("R:");

    mat Arot = A;
    double max_number_iterations = (double)n*(double)n*(double)n;
    int iteration = 0;
    double maxoff = offdiag(A, k, l, n);
    cout << "max off diagonal="<< maxoff << endl;

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
    cout << "Number of iteration: " << iteration << endl;
    cout << setprecision(10) << "Time used: " << time_used << " s at n=" << n << endl;
    //A.print("A:");
    //R.print("R:");

    // Test rotation:
    vec eigvals;
    eigvals = get_eigenvalues(A, n);
    //eigvals.print("Eigenvalues");

    double lmbd;
    for (int i=0; i<5; i++){
      if (i < 10){
        lmbd = 4.0*i + 3.0;
        cout << "i=" << i <<" Computed eigenvalue= " <<eigvals(i) << " exact= " << lmbd << endl;
      }
    }

    mat V; mat O;
    V = get_eigenvectors(A, R, n);
    //O = (V.t()*V);
    //O.print("orthogonal?");
    //V.print("Eigenvectors");

    // Write to file:

    ofile.open(outdata);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    int ex0 = 0; int ex1 = 1; int ex2 = 2;
    for(int i=0; i<n; i++){
      ofile << setw(15) << setprecision(8) << i*h << "  ";
      ofile << setw(15) << setprecision(8) << V(i,ex0) << "  "; // ground state
      ofile << setw(15) << setprecision(8) << V(i,ex1) << "  ";
      ofile << setw(15) << setprecision(8) << V(i,ex2) << "  ";
      ofile << "\n";
    }
    ofile.close();


    //TestEigenvalues(eigvals, n);
    //V.print("Eigen vectors:");
    //(V.t()*V).print("Orthogonality?");
  }
  return 0;
}

/*
For n = 400, rhomax=5
i=0 Computed eigenvalue= 2.99995 exact= 3
i=1 Computed eigenvalue= 6.99975 exact= 7
i=2 Computed eigenvalue= 10.9996 exact= 11
i=3 Computed eigenvalue= 15.0037 exact= 15
i=4 Computed eigenvalue= 19.0616 exact= 19
Compouted eigenvalue2.99995 equal exact = 3
Computed eigenvalue 6.99975 are differ to much from exact:7
Computed eigenvalue 10.9996 are differ to much from exact:11
Computed eigenvalue 15.0037 are differ to much from exact:15
Computed eigenvalue 19.0616 are differ to much from exact:19
Computed eigenvalue 23.3844 are differ to much from exact:23
Computed eigenvalue 28.3062 are differ to much from exact:27
Computed eigenvalue 34.0249 are differ to much from exact:31
Computed eigenvalue 40.5836 are differ to much from exact:35
Computed eigenvalue 47.9773 are differ to much from exact:39
Computed eigenvalue 56.1949 are differ to much from exact:43
Computed eigenvalue 65.2274 are differ to much from exact:47
Computed eigenvalue 75.0684 are differ to much from exact:51

*/
