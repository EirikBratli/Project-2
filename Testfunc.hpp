//#include "Jacobi_method.hpp"
#include <armadillo>
#include <iostream>

void TestEigenvalues(vec v , int n){
  /*
  Test eigenvalues for SL with one electron.
  */
  double diff, exact;
  double eps = 0.001;
  for (int i=0; i<n; i++){
    exact = 4.0*i + 3.0;
    diff = fabs(exact - v(i));
    if (diff > eps){
      cout << "Computed eigenvalue "<< v(i) << " are differ to much from exact:" << exact << endl;
    }
    else {
      cout << "Compouted eigenvalue" << v(i) << " equal exact = " << exact << endl;
    }
  }
}

void Orthogonality(mat R, mat V, int n, double eps){
  /*
  Using the vector/matrix properties and armadillo functions.
  Since V^T V = identity matrix
  */
  cout << "------------------" << endl;
  cout << "Testing orthogonality of the eigenvectors" << endl;
  int element1 =0; int element2 = 0;
  mat test1 = R.t()*R;
  mat test2 = V.t()*V;
  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++){
      if (i != j){
        if (test1(i,j) < eps){
            element1 += 1;
        }
        if (test2(i,j) < eps){
          element2 += 1;
        }
      }
    }
  }
  if (element1 == n*n - n && element2 == n*n - n){
    cout << "Passed test: Orthogonal eigenvectors" << endl;
    test2.print("V^T V = I");
  }
  else{
    cout << "Not passed: Eigenvectors are not orthogonal" << endl;
  }
}

void TestOffdiagonal(double maxoff, double eps){
  /*
  Test if the max off diagonal elements of matrix A are smaller than epsilon
  */
  cout << "------------------" << endl;
  cout << "Testing maximum value of off diagonal elements" << endl;
  cout << "Max off diagonal element = " << maxoff << endl;
  if (maxoff > eps){
    cout << "Not passed: Off diagonal elements of A are not zero" << endl;
  }
  else{
    cout << "Passed test: Off dialonal elements are zero" << endl;
  }
}
