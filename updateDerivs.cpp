// updateDerivs.cpp
// TWS, December 2019
// Updated May 2022 to include only "Model2020C" with optional 1D aorta
#include <math.h>
#include <stdio.h>
#include <string>
#include "nrutil.h"

using namespace std;

void gaussj(double **a, int n, double *b);

int updateDerivs(double t) {
	extern int Nmat, Nvars;
	extern double **J, **J1;
	extern double *f, *X, *dy;
	extern double Plv, Psp, Psv, Pra, Prv, Ppp, Ppv, Pla, Psa, Ppa;
	extern double q_la, q_lv, q_ra, q_rv;
	extern double q_sp, q_sv, q_pp, q_pv;
	extern double Csa, Csv, Cpa, Cpv, Cpp;
	extern double alpha_a1, beta_a1, gamma_a1, eta_a1, kappa_a1;
	extern double alpha_a2, beta_a2, gamma_a2, eta_a2, kappa_a2;
	extern double alpha_a3, beta_a3, gamma_a3, kappa_a3;
	extern double eta_la, kappa_la, chi_la;
	extern double eta_rv, kappa_rv, chi_rv;
	extern double eta_ra, kappa_ra, chi_ra;
	extern double B_lv, B_la, B_rv, B_ra, L_lv, L_la, L_rv, L_ra, L_pv, L_sv;
	extern double Rsp, Rsv, Rpa, Rpp, Rpv, Rra, f_rv;
	extern double Vrv, km_rv;
	extern double VA_aovo, VA_aovc, VA_mivo, VA_mivc, VA_puvo, VA_puvc, VA_tcvo, VA_tcvc, VL_aov, VL_miv, VL_puv, VL_tcv, rho;
	extern double Zmiv, Zaov, Ztcv, Zpuv, Kaovo, Kaovc, Kmivo, Kmivc, Kpuvo, Kpuvc, Ktcvo, Ktcvc;
	extern double Z0start, Z0end, P_ao_0, q_lv_0, P_echo, Q_trans;
	extern string SimpPulm;

	int i, k, index;
	double linerror, linerror2;
	double VA_aov, VA_miv, VA_puv, VA_tcv;
	double decayrate = 100.;
	double R0 = 0.0;

	//system for LV
	for (k = 0; k < Nmat; k++) for (i = 0; i < Nmat; i++) J[i][k] = 0;
	J[0][0] = alpha_a1;
	J[0][1] = beta_a1;
	J[0][2] = gamma_a1;
	J[0][3] = -eta_a1;
	J[1][0] = alpha_a2;
	J[1][1] = beta_a2;
	J[1][2] = gamma_a2;
	J[1][3] = -eta_a2;
	J[2][0] = alpha_a3;
	J[2][1] = beta_a3;
	J[2][2] = gamma_a3;
	J[3][0] = -eta_a1;
	J[3][1] = -eta_a2;
	f[0] = -kappa_a1;
	f[1] = -kappa_a2;
	f[2] = -kappa_a3;
	f[3] = -q_la + q_lv;
	for (i = 0; i < 4; i++) {
		X[i] = f[i];
		for (k = 0; k < 4; k++) {
			J1[i][k] = J[i][k];
			//printf(" %12.8f", J[i][k]);
		}
		//printf("   %12.8f\n", f[i]);
	}
	gaussj(J1, 4, X);	//using Gaussian elimination to solve system of (S109-112)
	dy[0] = X[0]; // da1_dt 
	dy[1] = X[1]; // da2_dt 
	dy[2] = X[2]; // da3_dt 
	Plv = X[3];

	linerror2 = 0.;
	for (i = 0; i < 4; i++) {	//check that result is correct
		linerror = f[i];
		for (k = 0; k < 4; k++) linerror -= J[i][k] * X[k];
		linerror2 += DSQR(linerror);
	}
	if (linerror2 > 1.e-12 || linerror2 != linerror2) {
		printf("*** Error in linear solver: LV\n");
		return 1;
	}
	//system for LA
	dy[3] = (q_pv - q_la) / eta_la;				//eq. (S113)
	Pla = (kappa_la + chi_la * dy[3]) / eta_la;	//eq. (S123)
	//pulmonary vein inertia
	if (SimpPulm == "on") dy[6] = (Ppp - Pla - Rpv * q_pv) / L_pv;	//eq. (S120)
	else dy[6] = (Ppv - Pla - Rpv * q_pv) / L_pv;
	//system for RV
	dy[4] = (q_ra - q_rv) / (eta_rv * f_rv);	//eq. (S114)
	Prv = (kappa_rv + chi_rv * dy[4]) / eta_rv;	//eq. (S124)
	//system for RA
	dy[5] = (q_sv - q_ra) / eta_ra;				//eq. (S114)
	Pra = (kappa_ra + chi_ra * dy[5]) / eta_ra;	//eq. (S125)
	//systemic vein inertia
	dy[7] = (Psv - Pra - Rsv * q_sv) / L_sv;	//eq. (S121)
	//valve equations
	//Valve areas
	VA_aov = Zaov * VA_aovo + (1. - Zaov) * VA_aovc;	//eq. (S104)
	VA_miv = Zmiv * VA_mivo + (1. - Zmiv) * VA_mivc;
	VA_puv = Zpuv * VA_puvo + (1. - Zpuv) * VA_puvc;
	VA_tcv = Ztcv * VA_tcvo + (1. - Ztcv) * VA_tcvc;
	// inertia coefficients (convert dyn/cm2 to kPa)
	B_lv = 1.e-4 * rho / (2. * DSQR(VA_aov));			//eq. (S103)
	B_la = 1.e-4 * rho / (2. * DSQR(VA_miv));
	B_rv = 1.e-4 * rho / (2. * DSQR(VA_puv));
	B_ra = 1.e-4 * rho / (2. * DSQR(VA_tcv));
	L_lv = 1.e-4 * rho * VL_aov / VA_aov;				//eq. (S102)
	L_la = 1.e-4 * rho * VL_miv / VA_miv;
	L_rv = 1.e-4 * rho * VL_puv / VA_puv;
	L_ra = 1.e-4 * rho * VL_tcv / VA_tcv;
	//pressures downstream of valves
	Psa = P_ao_0 + P_echo + Z0start * 1.e-4 * (q_lv - q_lv_0);	//subtract DC component. eq. (S100)
	Ppa = Ppp + Rpa * q_rv;						//eq. (S126)
	//flow derivatives
	// d / dt[q_lv] - aortic valve
	if (VA_aov > 1e-3) dy[8] = 1. / L_lv * (Plv - Psa - B_lv * fabs(q_lv) * q_lv);		//eq. (S101)
	else dy[8] = -decayrate * q_lv;
	// d / dt[q_la] - mitral valve
	if (VA_miv > 1e-3) dy[9] = 1. / L_la * (Pla - Plv - B_la * fabs(q_la) * q_la);		//eq. (S101)
	else dy[9] = -decayrate * q_la;
	// d / dt[q_rv]
	if (VA_puv > 1e-3) dy[10] = 1. / L_rv * (Prv - Ppa - B_rv * fabs(q_rv) * q_rv);		//eq. (S101)
	else dy[10] = -decayrate * q_rv;
	// d / dt[q_ra]
	if (VA_tcv > 1e-3) dy[11] = 1. / L_ra * (Pra - Prv - B_ra * fabs(q_ra) * q_ra);		//eq. (S101)
	else dy[11] = -decayrate * q_ra;
	// d / dt[Psp]
	dy[12] = 1. / Csa * (Q_trans - (Psp - Psv) / Rsp);	//eq. (S116)
	// d / dt[Psv]
	dy[13] = 1. / Csv * ((Psp - Psv) / Rsp - q_sv);							//eq. (S117)
	q_sp = (Psp - Psv) / Rsp;
	index = 15;
	// d / dt[Ppp]
	if (SimpPulm == "on") {
		dy[14] = 1. / Cpp * (q_rv - q_pv);										//eq. (S118)
		Ppv = Ppp - Rpv * q_pv;
		q_pp = 0.;
	}
	else{
	// d / dt[Ppp]
	dy[14] = 1. / Cpa * (q_rv - (Ppp - Ppv) / Rpp);							//eq. (S118)
	// d / dt[Ppv]
	dy[15] = 1. / Cpv * ((Ppp - Ppv) / Rpp - q_pv);							//eq. (S119)
	q_pp = (Ppp - Ppv) / Rpp;
		index++;
	}
	//derivatives of valve states
	if (Plv > Psa) dy[index] = (1. - Zaov) * Kaovo * (Plv - Psa);				//eq. (S105)
	else dy[index] = Zaov * Kaovc * (Plv - Psa);
	index++;
	if (Pla > Plv) dy[index] = (1. - Zmiv) * Kmivo * (Pla - Plv);				//eq. (S106)
	else dy[index] = Zmiv * Kmivc * (Pla - Plv);
	index++;
	if (Prv > Ppa) dy[index] = (1. - Zpuv) * Kpuvo * (Prv - Ppa);				//eq. (S107)
	else dy[index] = Zpuv * Kpuvc * (Prv - Ppa);
	index++;
	if (Pra > Prv) dy[index] = (1. - Ztcv) * Ktcvo * (Pra - Prv);				//eq. (S108)
	else dy[index] = Ztcv * Ktcvc * (Pra - Prv);

	return 0;
}
