//mathFunctions.cpp
// edited by TWS, November 2019, March 2026
#include <math.h>
#include <stdio.h>
#include "nrutil.h"

double simpsons(double* f, double dx, int Npoints) {
	double out = 0.;
	int k, N = Npoints;
	if (N >= 9) {	// Use simpsons rule if there are sufficient points
		out = 17 * f[0] + 59 * f[1] + 43 * f[2] + 49 * f[3]
			+ 49 * f[Npoints - 4] + 43 * f[Npoints - 3] + 59 * f[Npoints - 2] + 17 * f[Npoints - 1];
		for (k = 4; k < Npoints - 4; k++) out = out + 48 * f[k];
		out = out * dx / 48;
	}
	else {// Use trapezoid rule if there are fewer points
		for (k = 1; k < Npoints - 1; k++) out = out + f[k];
		out = out * 2 + f[0] + f[Npoints - 1];
		out = out * dx / 2;
	}
	return out;
}

double doubleIntegral(double** f, double dx, double dy, int Nx, int Ny) {
	int k;
	double value;
	double* integrand2;
	integrand2 = dvector(0, Nx - 1);
	for (k = 0; k < Nx; k++) integrand2[k] = simpsons(f[k], dy, Ny);
	value = simpsons(integrand2, dx, Nx);
	free_dvector(integrand2, 0, Nx - 1);
	return value;
}
