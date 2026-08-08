// initialgeometry.cpp
// TWS, December 2019

#include <math.h>
#include "nrutil.h"

double doubleIntegral(double** f, double dx, double dy, int Nx, int Ny);

void initialgeometry(int Nmu, int Nnu, int Nr) {
	extern double nu_up, muin0, muout0, a0, psi_in_b0, psi_out_b0, rin0, rout0, chin0;
	extern double rin0_la, rout0_la, rin0_rv, rout0_rv, rin0_ra, rout0_ra;
	extern double *mu0Vec, *nuVec, *psi_eq, *psiStar_eq, **omega;
	extern double **mu0, **nu;
	extern double **cps, **sps, **cps_sq, **sps_sq;
	extern double **snu, **cnu, **shmu0, **chmu0;
	extern double **snu_sq, **cnu_sq, **sum_sq0, **root_sum_sq0;
	extern double **ell0, **gmu0, **gnu0, **gphi0, **integratingFactor;
	extern double *r0_la, *r0_rv, *r0_ra;

	int k, j;
	double denom;

	chin0 = cosh(muin0);
	// Calculate mu and nu vectors of length Nmu, Nnu
	for (k = 0; k < Nnu; k++) nuVec[k] = nu_up + k * (PI - nu_up) / (Nnu - 1);
	for (k = 0; k < Nmu; k++) mu0Vec[k] = muin0 + k *(muout0 - muin0) / (Nmu - 1);
	double dmu = mu0Vec[1] - mu0Vec[0];
	double dnu = nuVec[1] - nuVec[0];
	// Set up initial spherical chamber values
	for (k = 0; k < Nr; k++) {
		r0_la[k] = rin0_la + k * (rout0_la - rin0_la) / (Nr - 1);
		r0_rv[k] = rin0_rv + k * (rout0_rv - rin0_rv) / (Nr - 1);
		r0_ra[k] = rin0_ra + k * (rout0_ra - rin0_ra) / (Nr - 1);
	}

	// Calculate fiber angles at nu = pi/2
	for (k = 0; k < Nmu; k++) {
		psi_eq[k] = (psi_in_b0 * (muout0 - mu0Vec[k]) + psi_out_b0 * (mu0Vec[k] - muin0)) / (muout0 - muin0);
		if (psi_eq[k] < 0) psiStar_eq[k] = PI + psi_eq[k];
		else psiStar_eq[k] = psi_eq[k];
	}

	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
		omega[k][j] = 1 / (tan(psiStar_eq[k]) * tanh(mu0Vec[k]));				//eq. (S45)
		mu0[k][j] = mu0Vec[k];
		nu[k][j] = nuVec[j];
		snu[k][j] = sin(nu[k][j]);
		cnu[k][j] = cos(nu[k][j]);
		snu_sq[k][j] = DSQR(snu[k][j]);
		cnu_sq[k][j] = DSQR(cnu[k][j]);
		shmu0[k][j] = sinh(mu0[k][j]);
		chmu0[k][j] = cosh(mu0[k][j]);
		sum_sq0[k][j] = DSQR(shmu0[k][j]) + snu_sq[k][j];
		root_sum_sq0[k][j] = sqrt(sum_sq0[k][j]);
		denom = sqrt(sum_sq0[k][j] + DSQR(shmu0[k][j] * snu_sq[k][j] * omega[k][j]));	//eq. (S46.3)
		cps[k][j] = shmu0[k][j] * snu_sq[k][j] * omega[k][j] / denom;	//eq. (S42)
		sps[k][j] = root_sum_sq0[k][j] / denom;	//eq. (S43)
		cps_sq[k][j] = DSQR(cps[k][j]);
		sps_sq[k][j] = DSQR(sps[k][j]);
		ell0[k][j] = a0 * denom;	//denominator in stretch ratio. eq. (S46.3)
		gmu0[k][j] = a0 * root_sum_sq0[k][j];
		gnu0[k][j] = gmu0[k][j];
		gphi0[k][j] = a0 * shmu0[k][j] * snu[k][j];
		integratingFactor[k][j] = gmu0[k][j] * gnu0[k][j] * gphi0[k][j];
	}
	// Test of integration
	//  Vwall_lv = 2. * PI * doubleIntegral(integratingFactor, dmu, dnu, Nmu, Nnu);
	// Exact result
	//  Vwall_lv = 2. / 3. * PI * DCUB(a0) *
	//	((cos(nu_up) + 1.) * cosh(muout0) * DSQR(sinh(muout0)) + (cosh(muout0) - 1.) * cos(nu_up) * DSQR(sin(nu_up))
	//		- (cos(nu_up) + 1.) * cosh(muin0) * DSQR(sinh(muin0)) - (cosh(muin0) - 1.) * cos(nu_up) * DSQR(sin(nu_up)));
}
