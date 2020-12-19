#pragma once
#ifndef SOLVER_H
#define SOLVER_H

#include <armadillo>
#include <cmath>

#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <iostream>

using namespace std;
using namespace arma;

class Solver {
public:
	void EulerForward(vec u,vec unew, int tsteps, int n, double dx, double dt, double alpha,int t1, string filename);
	void ForwardStep(vec& unew, vec u, double alpha, int n);
	void EulerForward2D(mat& v, int tsteps, int n, double dx, double dt, double alpha, int t1, string filename);
	void EulerBackward(vec u, vec unew, int tsteps, int n, double dx, double dt, double alpha, int t1, string filename);
	void CrankNicolson(vec u, vec unew, int tsteps, int n, double dx, double dt, double alpha, int t1, string filename);
	void write2file(vec u, string filename, int n, double dx, double dt, int tsteps);
	void write2file2D(mat v, string filename, int n, double dx, double dt, int tsteps);
	void tridiag(double a, double b, double c, vec u, vec& unew, int n);
	void D1_analytic(vec& analytic, double dt, double dx, int tsteps, int n);
	void D2_analytic(mat& analytic, double dt, double dx, int tsteps, int n);
	clock_t start, end;
	double pi = 3.14159;
};


#endif //SOLVER_H