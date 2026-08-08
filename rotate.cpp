// rotate.cpp
// edited by TWS, November 2019

void rotate_fiber_to_prolate_ff(double **Aff, double **A22, double **A23, double **A33) {
	extern int Nmu, Nnu;
	extern double **sps, **cps, **sps_sq, **cps_sq;

	int k, j;
	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
		A22[k][j] = sps_sq[k][j] * Aff[k][j];
		A23[k][j] = -cps[k][j] * sps[k][j] * Aff[k][j];
		A33[k][j] = cps_sq[k][j] * Aff[k][j];
	}
}

void rotate_fiber_to_prolate(double **Afib11, double **Afib12, double **Afib13,
	double **Afib22, double **Afib23, double **Afib33, double **A11, double **A12, double **A13,
	double **A22, double **A23, double **A33) {
	extern int Nmu, Nnu;
	extern double **sps, **cps, **sps_sq, **cps_sq;

	int k, j;
	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
		A11[k][j] = Afib11[k][j];
		A12[k][j] = sps[k][j] * Afib13[k][j] - cps[k][j] * Afib12[k][j];
		A13[k][j] = -sps[k][j] * Afib12[k][j] - cps[k][j] * Afib13[k][j];
		A22[k][j] = cps_sq[k][j] * Afib22[k][j] + sps[k][j] * (sps[k][j] * Afib33[k][j] - 2 * cps[k][j] * Afib23[k][j]);
		A23[k][j] = sps[k][j] * cps[k][j] * Afib22[k][j] + (cps_sq[k][j] - sps_sq[k][j])*Afib23[k][j]
			- sps[k][j] * cps[k][j] * Afib33[k][j];
		A33[k][j] = sps_sq[k][j] * Afib22[k][j] + cps[k][j] * (2 * sps[k][j] * Afib23[k][j] + cps[k][j] * Afib33[k][j]);
	}
}
