#include "solver.h"

int main(int argc, char* argv[]) {
	if (argc < 3) {
		cout << "Bad Usage: " << argv[0] <<
			"Too few arguments. Enter filename, function to run and size of n " << endl;
		cout << "(1 for Forward Euler, 2 for Euler Backward, 3 for Crank Nicolson, 4 for Forward Euler in two dimensions" << endl;
		exit(1);
	}

	string filename = argv[1];
	int func = atoi(argv[2]);
	int n = atoi(argv[3]);
	int tsteps = 500000; double dx = 1 / (double(n)); double dt = (1 / double(tsteps)); double alpha = dt / (dx * dx); int t1 = 10000;
	vec u = vec(n+1, fill::zeros);
	vec unew =  vec(n+1, fill::zeros);
	u(n) = 1.;
	Solver diffuse;
	
	if (func == 1) {
		diffuse.EulerForward(u, unew, tsteps, n+1, dx, dt, alpha, t1, filename);
	}
	else if (func == 2) {
		diffuse.EulerBackward(u, unew, tsteps, n + 1, dx, dt, alpha, t1, filename);
	}
	else if (func == 3) {
		diffuse.CrankNicolson(u, unew, tsteps, n + 1, dx, dt, alpha, t1, filename);
	}
	else if (func == 4) {
		mat v = mat(n+2, n+2, fill::zeros);
		v(n, n) = 1.;
		diffuse.EulerForward2D(v, tsteps, n + 1, dx, dt, alpha, t1, filename);
	}
	return 0;
}