//#include "Jacobi_method.hpp"
#include <armadillo>
#include <iostream>

void TestEigenvalues(vec v , int n){
  double diff, exact;
  double eps = 0.0001;
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

/*
mat TestOrthogonality(mat A, mat V, int n){
  int m = 4; double epsilon = 1e-10;
  mat Check = V.t()*V;

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
        if (fabs(V(i,j)) < epsilon){
          cout << "Ok, off-diagonal elements = 0, i,j= "<<i<<j << endl;
        }
        else{
          cout <<"i,j="<<i<<", "<<j<<" "<< setw(10)<<setprecision(15)<<V(i,j) <<" Too big difference from 0" << endl;
        }
      }
    }
  }

  return V;
}
*/
