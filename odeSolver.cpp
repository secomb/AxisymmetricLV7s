// odeSolver.cpp
// edited by TWS, November 2019
#include <stdio.h>
#include "nrutil.h"
#include <math.h>
#include <string>

#include <chrono>

using namespace std;

int updateModel(double t);

int ODEsolver() {
	extern int Nvars, Nmat, Nr, Ncycles, Nt, Nta, iRecord, Nmu, Nnu, Nstore, Ncycles_show;
	extern double Tc, Ta, dtRecord, dtMax;
	extern double Vrv, Vlv, Vla, Vra;
	extern double Csa, Csv, Cpa, Cpv, Cpp;
	extern double Plv, Pla, Ppa, Psa, Psp, Psv, Pra, Ppv, Ppp, Prv;
	extern double Rsp, Rsv, Rra, Rpa, Rpp, Rpv, Zmiv, Zaov, Ztcv, Zpuv;
	extern double a0, a1, a2, a3, a4, a5, a6, At, At_la;
	extern double q_ra, q_rv, q_la, q_lv;
	extern double q_sp, q_sv, q_pp, q_pv;
	extern double eta_a1, eta_a2, nu_up;
	extern double eta_la, eta_rv, f_rv, eta_ra;
	extern double VA_aovo, VA_aovc, VA_mivo, VA_mivc, VA_puvo, VA_puvc, VA_tcvo, VA_tcvc;
	extern double Pext, Psp0, Psv0, Ppp0, Ppv0;
	extern double lambda_mean, Sf_mean, Sf1_mean, Sf2_mean, Sf3_mean;

	extern double *ytest, *dy, *r;
	extern double *y, *tvec, *a1vec, *a2vec, *a3vec, *a4vec, *a5vec, *a6vec;
	extern double **store, **Y, ** chmu, ** shmu, **cnu;
	extern double *echo, *trans, *q_lv_hist, *p_ao_hist, P_ao_0, q_lv_0, P_echo, Q_trans, Vaorta, Z0end;

	extern string SimpPulm;

	int i, k, error;
	int Nmu2, Nnu2, Nnu_mid;
	double tRecord, t = 0.;
	double Amin = 0.25;	//echo doppler sampling area (cm2) in calculation of blood velocity
	double Vtot = 0., dVtot = 0., dVaorta;

	//initial values
	//y = [0.1, 0., 0., -1., 2., 0., 0., 0., 0., 0., 0., 0., 10., 0.80, 2, 1., 0., 1., 0., 1.]
	for (k = 0; k < Nvars; k++) y[k] = 0.;
	y[0] = -0.5;
	y[1] = -0.2;
	y[3] = -1;
	y[4] = 2.;
	y[12] = Psp0;
	y[13] = Psv0;
	y[14] = Ppp0;
	if (SimpPulm == "on") {
		Ppv = Ppv0;
		y[16] = 1.;
		y[18] = 1.;
	}
	else {
		y[15] = Ppv0;
		y[17] = 1.;
		y[19] = 1.;
	}

	P_ao_0 = 0.;
	q_lv_0 = 0.;
	P_echo = 0.;
	Plv = 0.;
	Prv = 0.;
	Pla = 0.;
	Pra = 0.;

	Nmu2 = Nmu / 2 - 1;
	Nnu2 = Nnu / 2 - 1;
	Nnu_mid = int((PI / 2. - nu_up) / (PI - nu_up) * (Nnu - 1) + 0.5);		//find index on equator
	   
	for (iRecord = 0; iRecord < Nt; iRecord++) {
		for (i = 0; i < Nstore; i++) store[i][iRecord] = 0.;
		tvec[iRecord] = 0.;
		a1vec[iRecord] = 0.;
		a2vec[iRecord] = 0.;
		a3vec[iRecord] = 0.;
		a4vec[iRecord] = 0.;
		a5vec[iRecord] = 0.;
		a6vec[iRecord] = 0.;
	}
	tRecord = (Ncycles - Ncycles_show) * Tc;
	iRecord = 0;

	for (i = 0; i < Nta; i++) {			// initialize history
		q_lv_hist[i] = 0.;
		p_ao_hist[i] = 0.;
	}
	//////////////////// start of time integration loop //////////////////
	while (t < Tc * Ncycles) {
	
		error = updateModel(t);

		Vtot = Vla + Vlv + Vra + Vrv + Csa * Psp + Csv * Psv;
		dVtot = eta_a1 * dy[0] + eta_a2 * dy[1] + eta_la * dy[3] + f_rv * eta_rv * dy[4]
			+ eta_ra * dy[5] + Csa * dy[12] + Csv * dy[13];
		if (SimpPulm == "on") {
			Vtot += Cpp * Ppp;
			dVtot += Cpp * dy[14];
		}
		else {
			Vtot += Cpa * Ppp + Cpv * Ppv;
			dVtot += Cpa * dy[14] + Cpv * dy[15];
		}
		dVaorta = q_lv - Q_trans;				// d/dt aorta volume
		Vaorta += dVaorta * dtMax;				// update aorta volume
		Vtot += Vaorta;
		dVtot += dVaorta;

		//update q_lv history vector, calculate reflected pressure
		for (i = Nta - 1; i > 0; i--) {
			q_lv_hist[i] = q_lv_hist[i - 1];	//shift contents by one position
			p_ao_hist[i] = p_ao_hist[i - 1];
		}
		q_lv_hist[0] = q_lv;
		p_ao_hist[0] = Psp;
		P_echo = 0.;
		Q_trans = 0.;
		for (i = 0; i < Nta; i++) {
			P_echo += echo[i] * q_lv_hist[i] * dtMax;	//convolution, eq. (S98)
			Q_trans += trans[i] * q_lv_hist[i] * dtMax;	//convolution, eq. (S99)
		}
		P_ao_0 += (p_ao_hist[0] - p_ao_hist[Nta - 1]) / Nta;	//update pressure average over one cycle at downstream end of aorta
		q_lv_0 += (q_lv_hist[0] - q_lv_hist[Nta - 1]) / Nta;	//update LV outflow average over one cycle
		
		// Store results, include only the last Ncycles_show
		if (t >= tRecord) {		
			if (iRecord >= Nt) {
				printf("*** Error: storage arrays not big enough ***\n");
				return 1;
			}
			for (k = 0; k < Nvars; k++)	Y[k][iRecord] = y[k];
			tvec[iRecord] = t;
			a1vec[iRecord] = a1;
			a2vec[iRecord] = a2;
			a3vec[iRecord] = a3;
			a4vec[iRecord] = a4;
			a5vec[iRecord] = a5;
			a6vec[iRecord] = a6;
			store[0][iRecord] = Vtot;
			store[1][iRecord] = Vlv;
			store[2][iRecord] = Plv + Pext;
			store[3][iRecord] = Vla;
			store[4][iRecord] = Pla + Pext;
			store[5][iRecord] = Vrv;
			store[6][iRecord] = Prv + Pext;
			store[7][iRecord] = Vra;
			store[8][iRecord] = Pra + Pext;
			store[10][iRecord] = Psa + Pext;
			store[11][iRecord] = Psp + Pext;
			store[12][iRecord] = Psv + Pext;
			store[13][iRecord] = Ppa + Pext;
			store[14][iRecord] = Ppp + Pext;
			if (SimpPulm == "on") store[15][iRecord] = Ppv + Pext;
			else store[15][iRecord] = Ppv - Rpv * q_pv + Pext;
			store[16][iRecord] = q_la;
			store[17][iRecord] = q_lv;
			store[18][iRecord] = q_ra;
			store[19][iRecord] = q_rv;
			store[21][iRecord] = q_sp;
			store[22][iRecord] = q_sv;
			store[24][iRecord] = q_pp;
			store[25][iRecord] = q_pv;
			store[26][iRecord] = Zmiv;
			store[27][iRecord] = Zaov;
			store[28][iRecord] = At;
			store[29][iRecord] = At_la;
			store[30][iRecord] = lambda_mean;
			store[31][iRecord] = Sf_mean + dy[0] * Sf1_mean + dy[1] * Sf2_mean + dy[2] * Sf3_mean;
			store[32][iRecord] = Sf_mean;
			store[33][iRecord] = Ztcv;
			store[34][iRecord] = Zpuv;
			//internal base to apex length and time derivative
			store[35][iRecord] = (a0 + a1) * (chmu[0][0] * cnu[0][0] + chmu[0][Nnu - 1]);
			if (iRecord == 0) store[36][iRecord] = 0.;
			else store[36][iRecord] = -(store[35][iRecord] - store[35][iRecord - 1]) / (tvec[iRecord] - tvec[iRecord - 1]);
			store[37][iRecord] = P_echo + Pext;
			store[38][iRecord] = Q_trans;
			store[39][iRecord] = P_ao_0 + Pext;
			store[40][iRecord] = q_lv_0;
			
			//blood flow velocities in valves
				//assume fully open if q > 0, fully closed if q < 0
				//use of actual valve areas causes artifacts due to time lag of flow changes
				//set minimum value of area to Amin
			if (q_lv > 0.) store[41][iRecord] = q_lv / DMAX(Amin, VA_aovo);
			else store[41][iRecord] = q_lv / DMAX(Amin, VA_aovc);
			if (q_la > 0.) store[42][iRecord] = q_la / DMAX(Amin, VA_mivo);
			else store[42][iRecord] = q_la / DMAX(Amin, VA_mivc);
			if (q_rv > 0.) store[43][iRecord] = q_rv / DMAX(Amin, VA_puvo);
			else store[43][iRecord] = q_rv / DMAX(Amin, VA_puvc);
			if (q_ra > 0.) store[44][iRecord] = q_ra / DMAX(Amin, VA_tcvo);
			else store[44][iRecord] = q_ra / DMAX(Amin, VA_tcvc);
			//internal diameter of LV at equator
			store[45][iRecord] = 2. * (a0 + a1) * shmu[0][Nnu_mid];	//corrected July 2025
			//external diameter of LV at equator
			store[46][iRecord] = 2. * (a0 + a1) * shmu[Nmu - 1][Nnu_mid];	//corrected July 2025
			//Cauchy fiber stress, including force-velocity effect
			store[47][iRecord] = Vaorta;
			tRecord += dtRecord;		// Reset record variables
			iRecord++;
		}
		//integration using RK2 scheme
		for (k = 0; k < Nvars; k++) {
			ytest[k] = y[k];
			y[k] += dtMax * dy[k] / 2.;
		}
		for (k = Nvars - 4; k < Nvars; k++) {	//avoid overshoot of valve states
			if (y[k] > 1.) y[k] = 1.;
			if (y[k] < 0.) y[k] = 0.;
		}
		error = updateModel(t + dtMax / 2.);
		for (k = 0; k < Nvars; k++) y[k] = ytest[k] + dtMax * dy[k];
		for (k = Nvars - 4; k < Nvars; k++) {	//avoid overshoot of valve states
			if (y[k] > 1.) y[k] = 1.;
			if (y[k] < 0.) y[k] = 0.;
		}
		t += dtMax;
	}
	return 0;
}
