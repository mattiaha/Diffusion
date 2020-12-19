#include "solver.h"
ofstream ofile;


void Solver::EulerBackward(vec u, vec unew, int tsteps, int n, double dx, double dt, double alpha, int t1, string filename) {
	string filename1 = "EB1" + filename;
	string filename2 = "EB2" + filename;
	double a, b, c;
	start = clock();

	
	a = c = -alpha;
	b = 1+2 * alpha;
	for (int t = 1; t < t1; t++) {
		tridiag(a, b, c, u, unew, n);
		unew(n - 1) = 1.;
		unew(0) = 0;
		u = unew;
	}
	vec u1 = unew;
	for (int t = t1; t <= tsteps; t++) {
		tridiag(a, b, c, u, unew, n);
		unew(n - 1) = 1.;
		unew(0) = 0;
		u = unew;
	}
	end = clock();
	cout << scientific << "EulerForward CPU time (sec) : " << ((double)end - (double)start) / CLOCKS_PER_SEC << endl;
	write2file(u1, filename1, n , dx, dt,t1);
	write2file(unew, filename2, n , dx, dt,tsteps);
}
void Solver::CrankNicolson(vec u, vec unew, int tsteps, int n, double dx, double dt, double alpha,int t1, string filename) {
	double a, b, c;
	alpha = alpha / 2.;
	a = c = -alpha;
	b = 1 + 2* alpha;
	string filename1 = "CN1" + filename;
	string filename2 = "CN2" + filename;
	unew(n - 1) = 1.;

	start = clock();
	for (int t = 1;t <t1;t++){
		ForwardStep(u, unew, alpha, n);
		u(0) = 0;
		u(n - 1) = 1.;
		tridiag(a, b, c, u, unew, n);
		unew(n - 1) = 1.;
		unew(0) = 0;
	}
	vec u1 = unew;
	for (int t = t1; t <= tsteps; t++) {
		ForwardStep(u, unew, alpha, n);
		u(0) = 0;
		u(n - 1) = 1.;
		tridiag(a, b, c, u, unew, n);
		unew(n - 1) = 1.;
		unew(0) = 0;
	}

	end = clock();
	cout << scientific << "EulerForward CPU time (sec) : " << ((double)end - (double)start) / CLOCKS_PER_SEC << endl;
	write2file(u1, filename1, n, dx,dt,t1);
	write2file(unew, filename2, n, dx, dt,tsteps);
}
void Solver::EulerForward(vec u, vec unew, int tsteps, int n, double dx, double dt, double alpha, int t1, string filename) {
	string filename1 = "EF1" + filename;
	string filename2 = "EF2" + filename;
	start = clock();


	if (alpha > 0.5) {
		cout << "Choose larger dx or lower dt" << endl;
		exit(1);

	}
	else {
		for (int t = 1; t < t1; t++) {
			ForwardStep(unew, u, alpha, n);
			unew(n - 1) = 1.;
			u = unew;
		}
		vec u1 = unew;

		for (int t = t1; t <= tsteps; t++) {
			ForwardStep(unew, u, alpha, n);
			unew(n - 1) = 1.;
			u = unew;
		}
		end = clock();
		cout << scientific << "EulerForward CPU time (sec) : " << ((double)end - (double)start) / CLOCKS_PER_SEC << endl;
		write2file(u1, filename1, n, dx, dt, t1);
		write2file(unew, filename2, n, dx, dt, tsteps);

	}
}
void Solver::EulerForward2D(mat& v, int tsteps, int n, double dx, double dt, double alpha, int t1, string filename) {
	string filename1 = "2dEF1" + filename;
	start = clock();
	mat v_new(n+1, n+1, fill::zeros);
	if (alpha > 0.5) {
		cout << "Choose larger dx or lower dt" << endl;
		exit(1);

	}
	else {
		
		
		for (int t = 1; t <= tsteps; t++) {

			for (int i = 1; i < n; i++) {
				
				for (int j = 1; j < n; j++) {
					v_new(i, j) = v(i, j) + alpha * (v(i + 1, j) + v(i - 1, j) + v(i, j + 1) + v(i, j - 1) - 4 * v(i, j));
				}
				v_new(n-1, n-1 ) = 1.;
				v_new(0, 0) = 0;
			}
			v_new(n-1 , n-1 ) = 1.;
			v_new(0, 0) = 0;
			v = v_new;
		}
		end = clock();
		cout << scientific << "EulerForward CPU time (sec) : " << ((double)end - (double)start) / CLOCKS_PER_SEC << endl;
		cout << v_new << endl;

	}
	write2file2D(v, filename, n, dx, dt, tsteps);
}

void Solver::write2file(vec u, string filename, int n, double dx, double dt, int tsteps) {
	vec analytic=vec(n,fill::zeros);
	D1_analytic(analytic, dt, dx, tsteps, n);
	ofile.open(filename);
	ofile << setiosflags(ios::showpoint | ios::uppercase);
	ofile << "     x:             analytic:               u(x):" << endl;
	for (int i = 0; i < n; i++) {
		ofile << setw(15) << setprecision(8) << dx * i ;
		ofile << setw(15) << setprecision(8) << analytic(i);
		ofile << setw(15) << setprecision(8) << u(i) << endl;
	}
	ofile.close();
}
void Solver::write2file2D(mat v, string filename, int n, double dx, double dt, int tsteps) {
	mat analytic(n,n, fill::zeros);
	D2_analytic(analytic, dt, dx, tsteps, n);
	cout << analytic << endl;
	ofile.open(filename);
	ofile << setiosflags(ios::showpoint | ios::uppercase);
	for (int i = 0; i < n+1; i++) {
		for (int j = 0; j < n+1; j++) {
			ofile << setw(15) << setprecision(8) << v(i, j) ;
			
		}
		cout << endl;
	}
	ofile.close();
	int m = (n-1) * (n-1);
	vec error(m);
	int a = 0;
	for (int i = 1; i < n; i++) {
		for (int j = 1; j < n; j++) {
			error(a) = abs((v(i,j)-analytic(i,j)/analytic(i,j)));
			a += 1;
		}
	}
	cout << "Mean error between numerical and analytical solution is:" << endl;
	cout << mean(error) << endl;
}

void Solver::ForwardStep(vec& unew, vec u, double alpha, int n) {
	for (int i = 1; i < n - 1; i++) {
		unew(i) = alpha * u(i - 1) + (1 - 2. * alpha) * u(i) + alpha * u(i + 1);
	}
}


void Solver::tridiag(double a, double b, double c, vec u, vec& unew, int n) {
	vec gam = vec(n);
	double bet;
	unew(0) = u(0) / b;
	bet = b;
	
	for (int i = 1; i < n; i++) {
		gam(i) = c / bet;
		
		bet = b - a * gam(i);
		unew(i) = (u(i) - a * unew(i - 1)) / bet;
		
	}
	
	for (int i = (n - 2); i >= 0; i--) {
		unew(i) -= gam(i + 1) * unew(i + 1);
	}
	
}

void Solver::D1_analytic(vec& analytic, double dt, double dx, int tsteps, int n) {
	double time = dt*tsteps;
	for (int i=0;i<n;i++){
		for (int m = 1; m < 10; m++) {
			analytic(i) += (2. / (m*pi)) * pow(-1, m) * sin(pi * dx * i * m)* exp(-(time * pow(m, 2) * pow(pi, 2)));
			
		}
		analytic(i) += dx * i;
	}
}

void Solver::D2_analytic(mat& analytic, double dt, double dx, int tsteps, int n) {
	double time = dt * tsteps;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int m = 1; m < 10; m++) {
				for (int k = 1; k < 10; k++) {
					analytic(i, j) += (4/(double(m*k))) *pow(-1,(m+k))*sin(m*pi*dx*i)*sin(k*pi*dx*j)*exp(-((m*m+k*k)*pi*pi*time));
				}
			}
			analytic(i, j) = analytic(i, j) / (pi * pi) + dx*i*dx*j;
			
		}
	}


}