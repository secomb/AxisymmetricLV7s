// updateMapping.cpp
#include <math.h>
#include <complex>
#include "nrutil.h"

using namespace std;
// Complex variable type definition "CX"
typedef std::complex<double> CX;

double simpsons(double* f, double dx, int Npoints);
double doubleIntegral(double** f, double dx, double dy, int Nx, int Ny);


void updateMapping() {
	extern int Nmu, Nnu;
	extern double a0, a1, a2, muin0, Vlv;
	extern double *nuVec, *vol_integrand;
	extern double **snu, **cnu, **snu_sq, **cnu_sq, **shmu0, **chmu0, **sum_sq0;
	extern double **shmu, **chmu, **sh_sq, **ch_sq, **sum_sq, **root_sum_sq, **dmu_dmu0, **dmu_dnu, **dmu_da1, **dmu_da2;
	extern double **ddmu_dnuda1, **ddmu_dnuda2, **ddmu_dmu0da1, **ddmu_dmu0da2;
	extern double *vIntegrand_da1, *vIntegrand_da2;
	extern double eta_a1, eta_a2;
	extern double** integrand, third, r3, chin0;

	int k, j, n;
	double a = a0 + a1;
	double asqr = DSQR(a);
	double a0sqr = DSQR(a0);
	double acube = DCUB(a);
	double a0cube = DCUB(a0);
	double ff2, ff3, dff2_da1, dff2_da2, dff3_da1, dff3_da2, ff1, ffac;
	double C0, F, RHS, b, LHS, Xmu, S, err;
	CX Sc, root3, muCX;

	// Update mu values using explicit solution for root of cubic
	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
		C0 = third - cnu_sq[k][j];
		F = a2 * a0sqr * (chin0 * (third * DSQR(chin0) - cnu_sq[k][j]) - C0);
		RHS = a0cube * (chmu0[k][j] * (third * DSQR(chmu0[k][j]) - cnu_sq[k][j]) - C0) - F;
		b = -cnu_sq[k][j];
		S = -(RHS / acube + C0);	// Intermediate root calculations are complex
		Xmu = chmu[k][j];
		LHS = Xmu * (DSQR(Xmu) * third + b) + S;
		err = abs(LHS);
		if (err > 0.01) {		//use exact solution
			Sc = S;
			root3 = pow(sqrt(4. * DCUB(b) + 9. * Sc * Sc) - 3. * Sc, third);
			chmu[k][j] = real(root3 / r3 - r3 * b / root3); // The imaginary part is 0, but still need to take only real part
		}
		else{				//use Newton iteration if possible
			n = 0;
			do {
				Xmu -= LHS / (DSQR(Xmu) + b);
				LHS = Xmu * (DSQR(Xmu) * third + b) + S;
				err = abs(LHS);
				n++;
			} while (err > 1.e-10 && n < 100);
			chmu[k][j] = Xmu;
		}
		if (chmu[k][j] < 1.) {		//Newton failed, use exact solution
			Sc = S;
			root3 = pow(sqrt(4. * DCUB(b) + 9. * Sc * Sc) - 3. * Sc, third);
			chmu[k][j] = real(root3 / r3 - r3 * b / root3);
		}
		if (chmu[k][j] < 1.)
			chmu[k][j] = 1.;	//just in case
		shmu[k][j] = sqrt(DSQR(chmu[k][j]) - 1.);
	}


	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
		sh_sq[k][j] = DSQR(shmu[k][j]);
		ch_sq[k][j] = DSQR(chmu[k][j]);
		sum_sq[k][j] = sh_sq[k][j] + snu_sq[k][j];
		root_sum_sq[k][j] = sqrt(sum_sq[k][j]);
		ffac = 2. * cnu[k][j] * snu[k][j];
		ff1 = a0cube * shmu0[k][j] * sum_sq0[k][j];					//eq. (S14)
		ff2 = ffac * (acube * (1. - chmu[k][j]) + a0cube * (chmu0[k][j] - 1.) - a2 * a0sqr * (chin0 - 1.));	//eq. (S15)
		ff3 = acube * shmu[k][j] * sum_sq[k][j];	//eq. (S16)

		dmu_dmu0[k][j] = ff1 / ff3;		//eq. (S13)
		dmu_dnu[k][j] = ff2 / ff3;		//eq. (S13)
		dmu_da1[k][j] = asqr * (1. - DCUB(chmu[k][j]) + (3. * chmu[k][j] - 3.) * cnu_sq[k][j]) / ff3;	//eq. (S24)
		dmu_da2[k][j] = a0sqr * (1. - DCUB(chin0) + (3. * chin0 - 3.) * cnu_sq[k][j]) / (3. * ff3);	//eq. (S25)
		dff2_da1 = ffac * (3. * asqr * (1. - chmu[k][j]) - acube * shmu[k][j] * dmu_da1[k][j]);	//eq. (S26)
		dff3_da1 = 3. * asqr * shmu[k][j] * sum_sq[k][j] + acube * chmu[k][j] * (3. * DSQR(shmu[k][j]) + DSQR(snu[k][j])) * dmu_da1[k][j];	//eq. (S27)
		dff2_da2 = -ffac * (acube * shmu[k][j] * dmu_da2[k][j] + a0sqr * (chin0 - 1.));		//eq. (S28)
		dff3_da2 = acube * chmu[k][j] * (3. * DSQR(shmu[k][j]) + DSQR(snu[k][j])) * dmu_da2[k][j];		//eq. (S29)
		ddmu_dnuda1[k][j] = (dff2_da1 * ff3 - ff2 * dff3_da1) / DSQR(ff3);		//eq. (S30)
		ddmu_dnuda2[k][j] = (dff2_da2 * ff3 - ff2 * dff3_da2) / DSQR(ff3);		//eq. (S30)
		ddmu_dmu0da1[k][j] = -ff1 / DSQR(ff3) * dff3_da1;						//eq. (S30)
		ddmu_dmu0da2[k][j] = -ff1 / DSQR(ff3) * dff3_da2;						//eq. (S30)
	}
	double dnu = nuVec[1] - nuVec[0];
	double c_in, coshmu, dcoshmu;
	for (j = 0; j < Nnu; j++) {
		c_in = chmu[0][j];
		dcoshmu = (c_in - 1.) / (Nmu - 1);
		for (k = 0; k < Nmu; k++) {
			coshmu = 1. + k * dcoshmu;
			integrand[k][j] = dcoshmu * acube * snu[0][j] * (DSQR(coshmu) - cnu_sq[0][j]);
		}
	}
	Vlv = 2. * PI * doubleIntegral(integrand, 1., dnu, Nmu, Nnu);

	eta_a1 = 0.;
	for (j = 0; j < Nnu; j++)
		vIntegrand_da2[j] = a0sqr * ((1. - DCUB(chin0)) / 3. + (chin0 - 1.) * cnu_sq[0][j]) * snu[0][j];	//Eq. (S71)
	eta_a2 = 2. * PI * simpsons(vIntegrand_da2, dnu, Nnu);
}
