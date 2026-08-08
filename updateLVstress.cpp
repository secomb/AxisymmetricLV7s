// updateLVstress.cpp
// edited by TWS, November 2019
#include <math.h>
#include <stdio.h>
#include "nrutil.h"

void rotate_fiber_to_prolate_ff(double **Aff, double **A22, double **A23, double **A33);
void rotate_fiber_to_prolate(double **Afib11, double **Afib12, double **Afib13,
	double **Afib22, double **Afib23, double **Afib33, double **A11, double **A12, double **A13,
	double **A22, double **A23, double **A33);
double simpsons(double* f, double dx, int Npoints);
double doubleIntegral(double** f, double dx, double dy, int Nx, int Ny);

void computeElastic() {
	extern int Nmu, Nnu;
	extern double c1, bff, bxx, bfx;
	extern double **Ess, **Esn, **Enn, **Esf, **Eff, **Enf;
	extern double **Sess, **Sesn, **Senn, **Sesf, **Seff, **Senf, **eW;
	extern double **Se11, **Se12, **Se13, **Se22, **Se23, **Se33;

	int k, j;
	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
		eW[k][j] = exp(bff * DSQR(Eff[k][j]) + bxx*(DSQR(Enn[k][j])
			+ DSQR(Ess[k][j]) + 2. * DSQR(Esn[k][j]))
			+ bfx*(2. * DSQR(Esf[k][j]) + 2. * DSQR(Enf[k][j])));
		Sess[k][j] = c1 * eW[k][j] * bxx * Ess[k][j];
		Senn[k][j] = c1 * eW[k][j] * bxx * Enn[k][j];
		Seff[k][j] = c1 * eW[k][j] * bff * Eff[k][j];
		Sesn[k][j] = c1 * eW[k][j] * bxx * Esn[k][j];
		Sesf[k][j] = c1 * eW[k][j] * bfx * Esf[k][j];
		Senf[k][j] = c1 * eW[k][j] * bfx * Enf[k][j];
	}
	rotate_fiber_to_prolate(Sess, Sesn, Sesf, Senn, Senf, Seff, Se11, Se12, Se13, Se22, Se23, Se33);
}

void computeViscous() {
	extern int Nmu, Nnu;
	extern double kv;
	extern double **Cinv11, **Cinv12, **Cinv13, **Cinv22, **Cinv23, **Cinv33;
	extern double **dE11_da1, **dE12_da1, **dE22_da1, **dE23_da1, **dE33_da1;
	extern double **dE11_da2, **dE12_da2, **dE22_da2, **dE22_da3, **dE23_da2, **dE23_da3, **dE33_da2;
	extern double **Sv11_a1, **Sv11_a2, **Sv11_a3, **Sv12_a1, **Sv12_a2, **Sv12_a3;
	extern double **Sv13_a1, **Sv13_a2, **Sv13_a3, **Sv22_a1, **Sv22_a2, **Sv22_a3;
	extern double **Sv23_a1, **Sv23_a2, **Sv23_a3, **Sv33_a1, **Sv33_a2, **Sv33_a3;

	int k, j; 
	double kvtwo = 2. * kv;
	double cof11, cof12, cof22, cof23, cof33;
	// Computes Sv with a1 derivatives
	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
		// Calculate Sv11
		cof11 = DSQR(Cinv11[k][j]);
		cof12 = 2. * Cinv12[k][j] * Cinv11[k][j];
		cof22 = DSQR(Cinv12[k][j]);
		cof33 = DSQR(Cinv13[k][j]);
		cof23 = 2. * Cinv12[k][j] * Cinv13[k][j];

		Sv11_a1[k][j] = kvtwo * (cof11 * dE11_da1[k][j] + cof12 * dE12_da1[k][j]
			+ cof22 * dE22_da1[k][j] + cof33 * dE33_da1[k][j] + cof23 * dE23_da1[k][j]);
		Sv11_a2[k][j] = kvtwo * (cof11 * dE11_da2[k][j] + cof12 * dE12_da2[k][j]
			+ cof22 * dE22_da2[k][j] + cof33 * dE33_da2[k][j] + cof23 * dE23_da2[k][j]);
		Sv11_a3[k][j] = kvtwo * (cof22 * dE22_da3[k][j] + cof23 * dE23_da3[k][j]); // 0 terms left out

		// Calculate Sv12
		cof12 = DSQR(Cinv12[k][j]) + Cinv11[k][j] * Cinv22[k][j];
		cof11 = Cinv11[k][j] * Cinv12[k][j];
		cof22 = Cinv22[k][j] * Cinv12[k][j];
		cof23 = Cinv23[k][j] * Cinv12[k][j] + Cinv13[k][j] * Cinv22[k][j];
		cof33 = Cinv13[k][j] * Cinv23[k][j];

		Sv12_a1[k][j] = kvtwo * (cof11 * dE11_da1[k][j] + cof12 * dE12_da1[k][j]
			+ cof22 * dE22_da1[k][j] + cof23 * dE23_da1[k][j] + cof33 * dE33_da1[k][j]);
		Sv12_a2[k][j] = kvtwo * (cof11 * dE11_da2[k][j] + cof12 * dE12_da2[k][j]
			+ cof22 * dE22_da2[k][j] + cof23 * dE23_da2[k][j] + cof33 * dE33_da2[k][j]);
		Sv12_a3[k][j] = kvtwo * (cof22 * dE22_da3[k][j] + cof23 * dE23_da3[k][j]);

		// Calculate Sv13
		cof11 = Cinv11[k][j] * Cinv13[k][j];
		cof12 = Cinv11[k][j] * Cinv23[k][j] + Cinv12[k][j] * Cinv13[k][j];
		cof22 = Cinv12[k][j] * Cinv23[k][j];
		cof23 = Cinv23[k][j] * Cinv13[k][j] + Cinv12[k][j] * Cinv33[k][j];
		cof33 = Cinv33[k][j] * Cinv13[k][j];

		Sv13_a1[k][j] = kvtwo * (cof11 * dE11_da1[k][j] + cof12 * dE12_da1[k][j]
			+ cof22 * dE22_da1[k][j] + cof23 * dE23_da1[k][j] + cof33 * dE33_da1[k][j]);
		Sv13_a2[k][j] = kvtwo * (cof11 * dE11_da2[k][j] + cof12 * dE12_da2[k][j]
			+ cof22 * dE22_da2[k][j] + cof23 * dE23_da2[k][j] + cof33 * dE33_da2[k][j]);
		Sv13_a3[k][j] = kvtwo * (cof22 * dE22_da3[k][j] + cof23 * dE23_da3[k][j]);

		// Calculate Sv22
		cof11 = DSQR(Cinv12[k][j]);
		cof12 = 2. * Cinv22[k][j] * Cinv12[k][j];
		cof22 = DSQR(Cinv22[k][j]);
		cof23 = 2. * Cinv22[k][j] * Cinv23[k][j];
		cof33 = DSQR(Cinv23[k][j]); 

		Sv22_a1[k][j] = kvtwo * (cof11 * dE11_da1[k][j] + cof12 * dE12_da1[k][j]
			+ cof22 * dE22_da1[k][j] + cof23 * dE23_da1[k][j] + cof33 * dE33_da1[k][j]);
		Sv22_a2[k][j] = kvtwo * (cof11 * dE11_da2[k][j] + cof12 * dE12_da2[k][j]
			+ cof22 * dE22_da2[k][j] + cof23 * dE23_da2[k][j] + cof33 * dE33_da2[k][j]);
		Sv22_a3[k][j] = kvtwo * (cof22 * dE22_da3[k][j] + cof23 * dE23_da3[k][j]);

		// Calculate Sv23
		cof11 = Cinv12[k][j] * Cinv13[k][j];
		cof12 = Cinv13[k][j] * Cinv22[k][j] + Cinv12[k][j] * Cinv23[k][j];
		cof22 = Cinv22[k][j] * Cinv23[k][j];
		cof23 = DSQR(Cinv23[k][j]) + Cinv22[k][j] * Cinv33[k][j];
		cof33 = Cinv33[k][j] * Cinv23[k][j];

		Sv23_a1[k][j] = kvtwo * (cof11 * dE11_da1[k][j] + cof12 * dE12_da1[k][j]
			+ cof22 * dE22_da1[k][j] + cof23 * dE23_da1[k][j] + cof33 * dE33_da1[k][j]);
		Sv23_a2[k][j] = kvtwo * (cof11 * dE11_da2[k][j] + cof12 * dE12_da2[k][j]
			+ cof22 * dE22_da2[k][j] + cof23 * dE23_da2[k][j] + cof33 * dE33_da2[k][j]);
		Sv23_a3[k][j] = kvtwo * (cof22 * dE22_da3[k][j] + cof23 * dE23_da3[k][j]);

		// Calculate Sv33
		cof11 = DSQR(Cinv13[k][j]);
		cof12 = 2. * Cinv23[k][j] * Cinv13[k][j];
		cof22 = DSQR(Cinv23[k][j]);
		cof23 = 2. * Cinv23[k][j] * Cinv33[k][j];
		cof33 = DSQR(Cinv33[k][j]);

		Sv33_a1[k][j] = kvtwo * (cof11 * dE11_da1[k][j] + cof12 * dE12_da1[k][j]
			+ cof22 * dE22_da1[k][j] + cof23 * dE23_da1[k][j] + cof33 * dE33_da1[k][j]);
		Sv33_a2[k][j] = kvtwo * (cof11 * dE11_da2[k][j] + cof12 * dE12_da2[k][j]
			+ cof22 * dE22_da2[k][j] + cof23 * dE23_da2[k][j] + cof33 * dE33_da2[k][j]);
		Sv33_a3[k][j] = kvtwo * (cof22 * dE22_da3[k][j] + cof23 * dE23_da3[k][j]);
	}
}

void computeActive() {
	extern int Nmu, Nnu;
	extern double Ls0, Lsmax, Lsw, km, kav, kp, At;
	extern double **ell, **ell0, **dell_da1, **dell_da2, **dell_da3, **eps_fed, **integrand, **integratingFactor;
	extern double **Sf_ff_const, **Sf_ff_a1, **Sf_ff_a2, **Sf_ff_a3;
	extern double **Sf22_const, **Sf22_a1, **Sf22_a2, **Sf22_a3;
	extern double **Sf23_const, **Sf23_a1, **Sf23_a2, **Sf23_a3;
	extern double **Sf33_const, **Sf33_a1, **Sf33_a2, **Sf33_a3;
	extern double *mu0Vec, *nuVec;

	int k, j;
	double activationkm, activationkv, G, LS;
	double twoLsw_sq = 2. * DSQR(Lsw);
	double dmu = mu0Vec[1] - mu0Vec[0];
	double dnu = nuVec[1] - nuVec[0];

	// Compute active stress in fiber direction
	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
		// Precomputations
		LS = Ls0 * ell[k][j] / ell0[k][j];
		G = exp(-DSQR(LS - Lsmax) / twoLsw_sq);

		// activation with sarcomere length
		activationkm = At * G * km * (1. + kp * eps_fed[k][j]);
		activationkv = At * G * kav * (1. + kp * eps_fed[k][j]);

		// Compute active stress (PK2)
		Sf_ff_const[k][j] = activationkm * DSQR(ell0[k][j] / ell[k][j]);
		Sf_ff_a1[k][j] = activationkv * dell_da1[k][j] / ell[k][j];
		Sf_ff_a2[k][j] = activationkv * dell_da2[k][j] / ell[k][j];
		Sf_ff_a3[k][j] = activationkv * dell_da3[k][j] / ell[k][j];
	}

	// Rotate back to prolate coordinates
	rotate_fiber_to_prolate_ff(Sf_ff_const, Sf22_const, Sf23_const, Sf33_const);
	rotate_fiber_to_prolate_ff(Sf_ff_a1, Sf22_a1, Sf23_a1, Sf33_a1);
	rotate_fiber_to_prolate_ff(Sf_ff_a2, Sf22_a2, Sf23_a2, Sf33_a2);
	rotate_fiber_to_prolate_ff(Sf_ff_a3, Sf22_a3, Sf23_a3, Sf33_a3);
}

//void totalStress(double da1_dt, double da2_dt, double da3_dt) {
//	extern int Nmu, Nnu;
//	extern double **E11, **E12, **E22, **E23, **E33;		//not sure these are all computed!!
//	extern double **Se11, **Se12, **Se21, **Se13, **Se31, **Se22, **Se23, **Se32, **Se33;
//	extern double **Sv11_a1, **Sv11_a2, **Sv11_a3, **Sv12_a1, **Sv12_a2, **Sv12_a3;
//	extern double **Sv31_a1, **Sv31_a2, **Sv31_a3, **Sv21_a1, **Sv21_a2, **Sv21_a3;
//	extern double **Sv13_a1, **Sv13_a2, **Sv13_a3, **Sv22_a1, **Sv22_a2, **Sv22_a3;
//	extern double **Sv23_a1, **Sv23_a2, **Sv23_a3, **Sv33_a1, **Sv33_a2, **Sv33_a3;
//	extern double **Sv32_a1, **Sv32_a2, **Sv32_a3;
//	extern double **Sf11_const, **Sf11_a1, **Sf11_a2, **Sf11_a3;	//check wwhere these are computed
//	extern double **Sf12_const, **Sf12_a1, **Sf12_a2, **Sf12_a3;
//	extern double **Sf21_const, **Sf21_a1, **Sf21_a2, **Sf21_a3;
//	extern double **Sf13_const, **Sf13_a1, **Sf13_a2, **Sf13_a3;
//	extern double **Sf22_const, **Sf22_a1, **Sf22_a2, **Sf22_a3;
//	extern double **Sf23_const, **Sf23_a1, **Sf23_a2, **Sf23_a3;
//	extern double **Sf32_const, **Sf32_a1, **Sf32_a2, **Sf32_a3;
//	extern double **Sf31_const, **Sf31_a1, **Sf31_a2, **Sf31_a3;
//	extern double **Sf33_const, **Sf33_a1, **Sf33_a2, **Sf33_a3;
//	extern double ***E, ***Se, ***Sf, ***Sv;
//
//	int k, j, ii;
//	// The order of the storage is E_11, E_12, E_13, E_21, E_22, E_23, E_31, E_32, E_33
//	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
//		E[0][k][j] = E11[k][j];
//		E[1][k][j] = E12[k][j];
//		E[2][k][j] = 0.;		// E_13 = 0. by definition
//		E[3][k][j] = E12[k][j];
//		E[4][k][j] = E22[k][j];
//		E[5][k][j] = E23[k][j];
//		E[6][k][j] = 0.;		// E_31 = 0. by definition
//		E[7][k][j] = E23[k][j];
//		E[8][k][j] = E33[k][j];
//	}
//
//	// The order of the storage is S_11, S_12, S_13, S_21, S_22, S_23, S_31, S_32, S_33
//	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
//		Se[0][k][j] = Se11[k][j];
//		Sf[0][k][j] = Sf11_const[k][j] + Sf11_a1[k][j] * da1_dt + Sf11_a2[k][j] * da2_dt + Sf11_a3[k][j] * da3_dt;
//		Sv[0][k][j] = Sv11_a1[k][j] * da1_dt + Sv11_a2[k][j] * da2_dt + Sv11_a3[k][j] * da3_dt;
//		Se[1][k][j] = Se12[k][j];
//		Sf[1][k][j] = Sf12_const[k][j] + Sf12_a1[k][j] * da1_dt + Sf12_a2[k][j] * da2_dt + Sf12_a3[k][j] * da3_dt;
//		Sv[1][k][j] = Sv12_a1[k][j] * da1_dt + Sv12_a2[k][j] * da2_dt + Sv12_a3[k][j] * da3_dt;
//		Se[2][k][j] = Se13[k][j];
//		Sf[2][k][j] = Sf13_const[k][j] + Sf13_a1[k][j] * da1_dt + Sf13_a2[k][j] * da2_dt + Sf13_a3[k][j] * da3_dt;
//		Sv[2][k][j] = Sv13_a1[k][j] * da1_dt + Sv13_a2[k][j] * da2_dt + Sv13_a3[k][j] * da3_dt;
//		Se[3][k][j] = Se21[k][j];
//		Sf[3][k][j] = Sf21_const[k][j] + Sf21_a1[k][j] * da1_dt + Sf21_a2[k][j] * da2_dt + Sf21_a3[k][j] * da3_dt;
//		Sv[3][k][j] = Sv21_a1[k][j] * da1_dt + Sv21_a2[k][j] * da2_dt + Sv21_a3[k][j] * da3_dt;
//		Se[4][k][j] = Se22[k][j];
//		Sf[4][k][j] = Sf22_const[k][j] + Sf22_a1[k][j] * da1_dt + Sf22_a2[k][j] * da2_dt + Sf22_a3[k][j] * da3_dt;
//		Sv[4][k][j] = Sv22_a1[k][j] * da1_dt + Sv22_a2[k][j] * da2_dt + Sv22_a3[k][j] * da3_dt;
//		Se[5][k][j] = Se23[k][j];
//		Sf[5][k][j] = Sf23_const[k][j] + Sf23_a1[k][j] * da1_dt + Sf23_a2[k][j] * da2_dt + Sf23_a3[k][j] * da3_dt;
//		Sv[5][k][j] = Sv23_a1[k][j] * da1_dt + Sv23_a2[k][j] * da2_dt + Sv23_a3[k][j] * da3_dt;
//		Se[6][k][j] = Se31[k][j];
//		Sf[6][k][j] = Sf31_const[k][j] + Sf31_a1[k][j] * da1_dt + Sf31_a2[k][j] * da2_dt + Sf31_a3[k][j] * da3_dt;
//		Sv[6][k][j] = Sv31_a1[k][j] * da1_dt + Sv31_a2[k][j] * da2_dt + Sv31_a3[k][j] * da3_dt;
//		Se[7][k][j] = Se32[k][j];
//		Sf[7][k][j] = Sf32_const[k][j] + Sf32_a1[k][j] * da1_dt + Sf32_a2[k][j] * da2_dt + Sf32_a3[k][j] * da3_dt;
//		Sv[7][k][j] = Sv32_a1[k][j] * da1_dt + Sv32_a2[k][j] * da2_dt + Sv32_a3[k][j] * da3_dt;
//		Se[8][k][j] = Se33[k][j];
//		Sf[8][k][j] = Sf33_const[k][j] + Sf33_a1[k][j] * da1_dt + Sf33_a2[k][j] * da2_dt + Sf33_a3[k][j] * da3_dt;
//		Sv[8][k][j] = Sv33_a1[k][j] * da1_dt + Sv33_a2[k][j] * da2_dt + Sv33_a3[k][j] * da3_dt;
//	}
//}

//void internalWork() {
//	extern int Nmu, Nnu;
//	extern double a0, a1, diff_work_e, diff_work_f, diff_work_v;
//	extern double *mu0Vec, *nuVec;
//	extern double **integrand, **integratingFactor;
//	extern double ***E, ***Eprev, ***Se, ***Se_prev, ***Sf, ***Sf_prev, ***Sv, ***Sv_prev;
//
//	int ii, k, j;
//	double a = a0 + a1;
//	double dmu = mu0Vec[1] - mu0Vec[0];
//	double dnu = nuVec[1] - nuVec[0];
//	double integ;
//
//	// Internal elastic work integral
//	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
//		integ = 0.;
//		for (ii = 0; ii < 9; ii++)
//			integ += (E[ii][k][j] - Eprev[ii][k][j]) * (Se[ii][k][j] + Se_prev[ii][k][j]);
//		integrand[k][j] = integ * integratingFactor[k][j] / 2.;
//	}
//	diff_work_e = 2. * PI * doubleIntegral(integrand, dmu, dnu, Nmu, Nnu);
//
//	// Internal fiber work integral
//	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
//		integ = 0.;
//		for (ii = 0; ii < 9; ii++)
//			integ += (E[ii][k][j] - Eprev[ii][k][j]) * (Sf[ii][k][j] + Sf_prev[ii][k][j]);
//		integrand[k][j] = integ * integratingFactor[k][j] / 2.;
//	}
//	diff_work_f = 2. * PI * doubleIntegral(integrand, dmu, dnu, Nmu, Nnu);
//
//	// Internal viscous work integral
//	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
//		integ = 0.;
//		for (ii = 0; ii < 9; ii++)
//			integ += (E[ii][k][j] - Eprev[ii][k][j]) * (Sv[ii][k][j] + Sv_prev[ii][k][j]);
//		integrand[k][j] = integ * integratingFactor[k][j] / 2.;
//	}
//	diff_work_v = 2. * PI * doubleIntegral(integrand, dmu, dnu, Nmu, Nnu);
//}
