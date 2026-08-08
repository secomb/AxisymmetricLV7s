// updateLVstrain.cpp
// edited by TWS, November 2019
// equation numbers refer to Moulton et al., 2017, including supplementary material (S)
#include <math.h>
#include <stdio.h>
#include "nrutil.h"

void updateDefGradTensor() {
	extern int Nmu, Nnu;
	extern double a0, a1, a2, a3;
	extern double **shmu0, **chmu0, **root_sum_sq0;
	extern double **shmu, **chmu, **sum_sq, **root_sum_sq;
	extern double **dmu_dmu0, **dmu_da1, **dmu_da2, **dmu_dnu;
	extern double **ddmu_dmu0da1, **ddmu_dmu0da2, **ddmu_dnuda1, **ddmu_dnuda2;
	extern double **snu_sq;
	extern double **F11, **F12, **F22, **F23, **F33;
	extern double **dF11_da1, **dF12_da1, **dF22_da1, **dF23_da1, **dF33_da1;
	extern double **dF11_da2, **dF12_da2, **dF22_da2, **dF23_da2, **dF23_da3, **dF33_da2;

	int k, j;
	double g, dg_da1, dg_da2;
	double a = a0 + a1;
	double gmu, gmu0, gphi, gphi0;		//note gnu = gmu

	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {		//deformation gradient tensor
		gmu0 = a0 * root_sum_sq0[k][j];
		gmu = a * root_sum_sq[k][j];
		gphi0 = a0 * shmu0[k][j];		//not including factor of sin(nu)
		gphi = a * shmu[k][j] * snu_sq[k][j];
		g = gmu / gmu0;
		dg_da1 = (sum_sq[k][j] + a * chmu[k][j] * shmu[k][j] * dmu_da1[k][j]) / (gmu0 * root_sum_sq[k][j]);	//eq. (S31)
		dg_da2 = a * chmu[k][j] * shmu[k][j] * dmu_da2[k][j] / (gmu0 * root_sum_sq[k][j]);	//eq. (S32)

		F11[k][j] = g * dmu_dmu0[k][j];	//eq. (S8)
		dF11_da1[k][j] = dg_da1 * dmu_dmu0[k][j] + g * ddmu_dmu0da1[k][j];
		dF11_da2[k][j] = dg_da2 * dmu_dmu0[k][j] + g * ddmu_dmu0da2[k][j];

		F12[k][j] = g * dmu_dnu[k][j];	//eq. (S9)
		dF12_da1[k][j] = dg_da1 * dmu_dnu[k][j] + g * ddmu_dnuda1[k][j];
		dF12_da2[k][j] = dg_da2 * dmu_dnu[k][j] + g * ddmu_dnuda2[k][j];

		F22[k][j] = g;					//eq. (S10)
		dF22_da1[k][j] = dg_da1;
		dF22_da2[k][j] = dg_da2;

		F23[k][j] = -a3 * gphi / gmu0;	//eq. (S11)
		dF23_da1[k][j] = -a3 * snu_sq[k][j] * (shmu[k][j] + a * chmu[k][j] * dmu_da1[k][j]) / gmu0;	//eq. (S33)
		dF23_da2[k][j] = -a * a3 * chmu[k][j] * snu_sq[k][j] * dmu_da2[k][j] / gmu0;	//eq. (s33)
		dF23_da3[k][j] = -gphi / gmu0;

		F33[k][j] = a * shmu[k][j] / gphi0;		//eq. (S12)
		dF33_da1[k][j] = (shmu[k][j] + a * chmu[k][j] * dmu_da1[k][j]) / gphi0;	//eq. (S34)
		dF33_da2[k][j] = a * chmu[k][j] * dmu_da2[k][j] / gphi0;	//eq. (S34)
	}
}

void updateCauchyTensor() {	//right Cauchy-Green tensor
	extern int Nmu, Nnu;
	extern double **F11, **F12, **F22, **F23, **F33;
	extern double **C11, **C12, **C22, **C23, **C33;
	extern double **Cinv11, **Cinv12, **Cinv13, **Cinv22, **Cinv23, **Cinv33;

	int k, j;
	double detC = 0.;
	double byDetC = 0.;

	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
		C11[k][j] = DSQR(F11[k][j]);		//C = F^T F
		C12[k][j] = F11[k][j] * F12[k][j];
		C22[k][j] = DSQR(F12[k][j]) + DSQR(F22[k][j]) + DSQR(F23[k][j]);
		C23[k][j] = F33[k][j] * F23[k][j];
		C33[k][j] = DSQR(F33[k][j]);
		detC = C11[k][j] * C22[k][j] * C33[k][j] - C11[k][j] * DSQR(C23[k][j]) - DSQR(C12[k][j]) * C33[k][j];
		byDetC = 1. / detC;
		Cinv11[k][j] = byDetC * (C22[k][j] * C33[k][j] - DSQR(C23[k][j]));
		Cinv12[k][j] = byDetC * (-C12[k][j] * C33[k][j]);
		Cinv13[k][j] = byDetC * (C12[k][j] * C23[k][j]);
		Cinv22[k][j] = byDetC * (C11[k][j] * C33[k][j]);
		Cinv23[k][j] = byDetC * (-C11[k][j] * C23[k][j]);
		Cinv33[k][j] = byDetC * (C11[k][j] * C22[k][j] - DSQR(C12[k][j]));
	}
}

void updateGreenStrain() {	//Green-Lagrange strain tensor
	extern int Nmu, Nnu;
	extern double **E11, **E12, **E22, **E23, **E33;
	extern double **C11, **C12, **C22, **C23, **C33;
	extern double **Cinv11, **Cinv12, **Cinv13, **Cinv22, **Cinv23, **Cinv33;
	extern double **F11, **F12, **F22, **F23, **F33;
	extern double **dF11_da1, **dF12_da1, **dF22_da1, **dF23_da1, **dF33_da1;
	extern double **dF11_da2, **dF12_da2, **dF22_da2, **dF23_da2, **dF23_da3, **dF33_da2;
	extern double **dE11_da1, **dE12_da1, **dE22_da1, **dE23_da1, **dE33_da1;
	extern double **dE11_da2, **dE12_da2, **dE22_da2, **dE22_da3, **dE23_da2, **dE23_da3, **dE33_da2;
	extern double **sps, **cps, **cps_sq, **sps_sq;
	extern double **Ess, **Esn, **Enn, **Esf, **Eff, **Enf;

	int k, j;

	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
		E11[k][j] = 0.5*(C11[k][j] - 1);	//E = (1/2) (C - 1)
		E12[k][j] = 0.5*C12[k][j];
		E22[k][j] = 0.5*(C22[k][j] - 1);
		E23[k][j] = 0.5*C23[k][j];
		E33[k][j] = 0.5*(C33[k][j] - 1);
		dE11_da1[k][j] = F11[k][j] * dF11_da1[k][j];
		dE11_da2[k][j] = F11[k][j] * dF11_da2[k][j];
		dE12_da1[k][j] = 0.5*(F12[k][j] * dF11_da1[k][j] + F11[k][j] * dF12_da1[k][j]);
		dE12_da2[k][j] = 0.5*(F12[k][j] * dF11_da2[k][j] + F11[k][j] * dF12_da2[k][j]);
		dE22_da1[k][j] = F12[k][j] * dF12_da1[k][j] + F22[k][j] * dF22_da1[k][j] + F23[k][j] * dF23_da1[k][j];
		dE22_da2[k][j] = F12[k][j] * dF12_da2[k][j] + F22[k][j] * dF22_da2[k][j] + F23[k][j] * dF23_da2[k][j];
		dE22_da3[k][j] = F23[k][j] * dF23_da3[k][j];
		dE23_da1[k][j] = 0.5*(F23[k][j] * dF33_da1[k][j] + F33[k][j] * dF23_da1[k][j]);
		dE23_da2[k][j] = 0.5*(F23[k][j] * dF33_da2[k][j] + F33[k][j] * dF23_da2[k][j]);
		dE23_da3[k][j] = 0.5*F33[k][j] * dF23_da3[k][j];
		dE33_da1[k][j] = F33[k][j] * dF33_da1[k][j];
		dE33_da2[k][j] = F33[k][j] * dF33_da2[k][j];
	}
	//for loops repaired by TWS, November 2019
	// Rotate E_prolate to E_fiber for use in the elastic stress computation
	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
		Ess[k][j] = E11[k][j];
		Esn[k][j] = -cps[k][j] * E12[k][j]; // E13 is zero 
		Esf[k][j] = sps[k][j] * E12[k][j]; // E13 is zero 
		Enn[k][j] = cps_sq[k][j] * E22[k][j]
			+ sps[k][j] * (sps[k][j] * E33[k][j] + 2 * cps[k][j] * E23[k][j]);
		Enf[k][j] = -sps[k][j] * cps[k][j] * E22[k][j] + (cps_sq[k][j] - sps_sq[k][j])*E23[k][j]
			+ sps[k][j] * cps[k][j] * E33[k][j];
		Eff[k][j] = sps_sq[k][j] * E22[k][j] + cps[k][j] * (cps[k][j] * E33[k][j] - 2 * sps[k][j] * E23[k][j]);
	}
}

void updateFiberStrain(double t) {
	extern int Nmu, Nnu;
	extern double a0, a1, a2, a3, Tc, Ta, dtMax;
	extern double **shmu, **chmu, **sh_sq, **snu_sq, **sum_sq, **ell, **ell0, **eps_f, **eps_fed, **dmu_dnu;
	extern double **dell_da1, **dell_da2, **dell_da3, **omega;
	extern double **dmu_da1, **dmu_da2, **ddmu_dnuda1, **ddmu_dnuda2;

	int k, j;
	double a = a0 + a1;
	double denom, diffac, omegafac;

	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
		diffac = 1. + DSQR(dmu_dnu[k][j]);
		omegafac = DSQR(a3 + omega[k][j]);
		denom = sqrt(diffac * sum_sq[k][j] + omegafac * DSQR(shmu[k][j] * snu_sq[k][j]));	//eq. (S46.3)
		ell[k][j] = a * denom;																//eq. (S46.3)
		eps_f[k][j] = 0.5 * (DSQR(ell[k][j]) / DSQR(ell0[k][j]) - 1.);	//fiber strain, after eq. (S70)
		dell_da1[k][j] = denom + a * (sum_sq[k][j] * dmu_dnu[k][j] * ddmu_dnuda1[k][j]
			+ (diffac + omegafac * DSQR(snu_sq[k][j])) * shmu[k][j] * chmu[k][j] * dmu_da1[k][j]) / denom;
		dell_da2[k][j] = a * (sum_sq[k][j] * dmu_dnu[k][j] * ddmu_dnuda2[k][j]
			+ (diffac + omegafac * DSQR(snu_sq[k][j]))	* shmu[k][j] * chmu[k][j] * dmu_da2[k][j]) / denom;
		dell_da3[k][j] = a * (a3 + omega[k][j]) * DSQR(snu_sq[k][j]) * sh_sq[k][j] / denom;
	}
	// Store end-diastolic strain at appropriate time 
	double truncT = t - (int)(t / Tc) * Tc;
	if (fabs(truncT - Tc + Ta) < dtMax) for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) eps_fed[k][j] = eps_f[k][j];
}
