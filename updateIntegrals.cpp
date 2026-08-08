// equilibriumIntegrals.cpp
// edited by TWS, November 2019
#include <stdio.h>
#include <math.h>
#include <string>
#include "nrutil.h"

using namespace std;

double doubleIntegral(double** f, double dx, double dy, int Nx, int Ny);

void updateIntegrals(double t) {
	extern int Nmu, Nnu;
	extern double alpha_a1, alpha_a2, alpha_a3;
	extern double beta_a1, beta_a2, beta_a3;
	extern double gamma_a1, gamma_a2, gamma_a3;
	extern double kappa_a1, kappa_a2, kappa_a3;
	extern double *mu0Vec, *nuVec;
	extern double **Se11, **Se12, **Se13, **Se22, **Se23, **Se33;
	extern double **dE11_da1, **dE12_da1, **dE22_da1, **dE23_da1, **dE33_da1;
	extern double **dE11_da2, **dE12_da2, **dE22_da2, **dE22_da3, **dE23_da2, **dE23_da3, **dE33_da2;
	extern double **Sf22_const, **Sf23_const, **Sf33_const;
	extern double **Sf22_a1, **Sf23_a1, **Sf33_a1;
	extern double **Sf22_a2, **Sf23_a2, **Sf33_a2;
	extern double **Sf22_a3, **Sf23_a3, **Sf33_a3;
	extern double **Sv11_a1, **Sv11_a2, **Sv11_a3, **Sv12_a1, **Sv12_a2, **Sv12_a3;
	extern double **Sv13_a1, **Sv13_a2, **Sv13_a3, **Sv22_a1, **Sv22_a2, **Sv22_a3;
	extern double **Sv23_a1, **Sv23_a2, **Sv23_a3, **Sv33_a1, **Sv33_a2, **Sv33_a3;
	extern double **integrand, **integratingFactor, volume;
	extern double ** ell, ** ell0, lambda_mean;
	extern double** Sf_ff_const, ** Sf_ff_a1, ** Sf_ff_a2, ** Sf_ff_a3, Sf_mean, Sf1_mean, Sf2_mean, Sf3_mean;
	extern string CalcFiber;

	int k, j;
	double val;
	double dmu = mu0Vec[1] - mu0Vec[0];
	double dnu = nuVec[1] - nuVec[0];
	double integ;

	// KAPPA INTEGRALS
	// Compute kappa_a1 integral
	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
		integ = dE11_da1[k][j] * Se11[k][j]
			+ dE22_da1[k][j] * (Se22[k][j] + Sf22_const[k][j])
			+ dE33_da1[k][j] * (Se33[k][j] + Sf33_const[k][j])
			+ 2. * dE12_da1[k][j] * Se12[k][j]
			+ 2. * dE23_da1[k][j] * (Se23[k][j] + Sf23_const[k][j]);		//eq. (S77)
		integrand[k][j] = integ * integratingFactor[k][j];
	}
	val = doubleIntegral(integrand, dmu, dnu, Nmu, Nnu);
	kappa_a1 = 2.*PI*val;

	// Compute kappa_a2 integral
	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
		integ = dE11_da2[k][j] * Se11[k][j]
			+ dE22_da2[k][j] * (Se22[k][j] + Sf22_const[k][j])
			+ dE33_da2[k][j] * (Se33[k][j] + Sf33_const[k][j])
			+ 2. * dE12_da2[k][j] * Se12[k][j] 
			+ 2 * dE23_da2[k][j] * (Se23[k][j] + Sf23_const[k][j]);		//eq. (S77)
		integrand[k][j] = integ * integratingFactor[k][j];
	}
	val = doubleIntegral(integrand, dmu, dnu, Nmu, Nnu);
	kappa_a2 = 2.*PI*val;

	// Compute kappa_a3 integral
	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
		integ = dE22_da3[k][j] * (Se22[k][j] + Sf22_const[k][j])
			+ 2. * dE23_da3[k][j] * (Se23[k][j] + Sf23_const[k][j]);		//eq. (S77)
		integrand[k][j] = integ * integratingFactor[k][j];
	}
	val = doubleIntegral(integrand, dmu, dnu, Nmu, Nnu);
	kappa_a3 = 2.*PI*val;

	// ALPHA INTEGRALS
	// Compute alpha_a1 integral
	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
		integ = dE11_da1[k][j] * Sv11_a1[k][j]
			+ dE22_da1[k][j] * (Sv22_a1[k][j] + Sf22_a1[k][j])
			+ dE33_da1[k][j] * (Sv33_a1[k][j] + Sf33_a1[k][j])
			+ 2. * dE12_da1[k][j] * Sv12_a1[k][j] 
			+ 2 * dE23_da1[k][j] * (Sv23_a1[k][j] + Sf23_a1[k][j]);	//eq. (S77)
		integrand[k][j] = integ * integratingFactor[k][j];
	}
	val = doubleIntegral(integrand, dmu, dnu, Nmu, Nnu);
	alpha_a1 = 2.*PI*val;

	// Compute alpha_a2 integral
	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
		integ = dE11_da2[k][j] * Sv11_a1[k][j]
			+ dE22_da2[k][j] * (Sv22_a1[k][j] + Sf22_a1[k][j])
			+ dE33_da2[k][j] * (Sv33_a1[k][j] + Sf33_a1[k][j])
			+ 2. * dE12_da2[k][j] * Sv12_a1[k][j]
			+ 2. * dE23_da2[k][j] * (Sv23_a1[k][j] + Sf23_a1[k][j]);	//eq. (S77)
		integrand[k][j] = integ * integratingFactor[k][j];
	}
	val = doubleIntegral(integrand, dmu, dnu, Nmu, Nnu);
	alpha_a2 = 2.*PI*val;

	// Compute alpha_a3 integral
	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
		integ = dE22_da3[k][j] * (Sv22_a1[k][j] + Sf22_a1[k][j]) + 2 * dE23_da3[k][j] * (Sv23_a1[k][j] + Sf23_a1[k][j]);	//eq. (S77)
		integrand[k][j] = integ * integratingFactor[k][j];
	}
	val = doubleIntegral(integrand, dmu, dnu, Nmu, Nnu);
	alpha_a3 = 2.*PI*val;

	// BETA INTEGRALS
	// Compute beta_a1 integral
	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
		integ = dE11_da1[k][j] * Sv11_a2[k][j]
			+ dE22_da1[k][j] * (Sv22_a2[k][j] + Sf22_a2[k][j])
			+ dE33_da1[k][j] * (Sv33_a2[k][j] + Sf33_a2[k][j])
			+ 2. * dE12_da1[k][j] * Sv12_a2[k][j]
			+ 2. * dE23_da1[k][j] * (Sv23_a2[k][j] + Sf23_a2[k][j]);	//eq. (S77)
		integrand[k][j] = integ * integratingFactor[k][j];
	}
	val = doubleIntegral(integrand, dmu, dnu, Nmu, Nnu);
	beta_a1 = 2.*PI*val;

	// Compute beta_a2 integral
	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
		integ = dE11_da2[k][j] * Sv11_a2[k][j] + 2. * dE12_da2[k][j] * Sv12_a2[k][j] + dE22_da2[k][j] * (Sv22_a2[k][j] + Sf22_a2[k][j])
			+ 2. * dE23_da2[k][j] * (Sv23_a2[k][j] + Sf23_a2[k][j]) + dE33_da2[k][j] * (Sv33_a2[k][j] + Sf33_a2[k][j]);	//eq. (S77)
		integrand[k][j] = integ * integratingFactor[k][j];
	}
	val = doubleIntegral(integrand, dmu, dnu, Nmu, Nnu);
	beta_a2 = 2.*PI*val;

	// Compute beta_a3 integral
	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
		integ = dE22_da3[k][j] * (Sv22_a2[k][j] + Sf22_a2[k][j]) + 2. * dE23_da3[k][j] * (Sv23_a2[k][j] + Sf23_a2[k][j]);	//eq. (S77)
		integrand[k][j] = integ * integratingFactor[k][j];
	}
	val = doubleIntegral(integrand, dmu, dnu, Nmu, Nnu);
	beta_a3 = 2.*PI*val;

	// GAMMA INTEGRALS
	// Compute gamma_a1 integral
	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
		integ = dE11_da1[k][j] * Sv11_a3[k][j] + 2 * dE12_da1[k][j] * Sv12_a3[k][j] + dE22_da1[k][j] * (Sv22_a3[k][j] + Sf22_a3[k][j])
			+ 2 * dE23_da1[k][j] * (Sv23_a3[k][j] + Sf23_a3[k][j]) + dE33_da1[k][j] * (Sv33_a3[k][j] + Sf33_a3[k][j]);		//eq. (S77)
		integrand[k][j] = integ * integratingFactor[k][j];
	}
	val = doubleIntegral(integrand, dmu, dnu, Nmu, Nnu);
	gamma_a1 = 2.*PI*val;

	// Compute gamma_a2 integral
	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
		integ = dE11_da2[k][j] * Sv11_a3[k][j] + 2. * dE12_da2[k][j] * Sv12_a3[k][j] + dE22_da2[k][j] * (Sv22_a3[k][j] + Sf22_a3[k][j])
			+ 2. * dE23_da2[k][j] * (Sv23_a3[k][j] + Sf23_a3[k][j]) + dE33_da2[k][j] * (Sv33_a3[k][j] + Sf33_a3[k][j]);		//eq. (S77)
		integrand[k][j] = integ * integratingFactor[k][j];
	}
	val = doubleIntegral(integrand, dmu, dnu, Nmu, Nnu);
	gamma_a2 = 2.*PI*val;

	// Compute gamma_a3 integral
	for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
		integ = dE22_da3[k][j] * (Sv22_a3[k][j] + Sf22_a3[k][j]) + 2. * dE23_da3[k][j] * (Sv23_a3[k][j] + Sf23_a3[k][j]);		//eq. (S77)
		integrand[k][j] = integ * integratingFactor[k][j];
	}
	val = doubleIntegral(integrand, dmu, dnu, Nmu, Nnu);
	gamma_a3 = 2.*PI*val;

	// Compute mean stretch and active force on LV
	if (CalcFiber == "on") {
		volume = doubleIntegral(integratingFactor, dmu, dnu, Nmu, Nnu);
		for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++) {
			integ = ell[k][j] / ell0[k][j];
			integrand[k][j] = integ * integratingFactor[k][j];
		}
		lambda_mean = doubleIntegral(integrand, dmu, dnu, Nmu, Nnu) / volume;
		for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++)
			integrand[k][j] = Sf_ff_const[k][j] * integratingFactor[k][j];
		Sf_mean = doubleIntegral(integrand, dmu, dnu, Nmu, Nnu) / volume;
		for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++)
			integrand[k][j] = Sf_ff_a1[k][j] * integratingFactor[k][j];
		Sf1_mean = doubleIntegral(integrand, dmu, dnu, Nmu, Nnu) / volume;
		for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++)
			integrand[k][j] = Sf_ff_a2[k][j] * integratingFactor[k][j];
		Sf2_mean = doubleIntegral(integrand, dmu, dnu, Nmu, Nnu) / volume;
		for (k = 0; k < Nmu; k++) for (j = 0; j < Nnu; j++)
			integrand[k][j] = Sf_ff_a3[k][j] * integratingFactor[k][j];
		Sf3_mean = doubleIntegral(integrand, dmu, dnu, Nmu, Nnu) / volume;
	}
	else {
		lambda_mean = 0.;
		Sf_mean = 0.;
		Sf1_mean = 0.;
		Sf2_mean = 0.;
		Sf3_mean = 0.;
	}
}
