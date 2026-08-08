// updateParams - using values stored in Param array
// TWS, July 2025

#define _CRT_SECURE_NO_WARNINGS
#include <math.h>
#include <complex>
#include "nrutil.h"

using namespace std;

typedef std::complex<double> CX;

void updateParams(double* Param, int Nparam)
{
	extern int Nvars, Nmat, iRecord, Nstore, Nr, Ncycles, Nt, Nta, Ncycles_show, aortaVar;
	extern double Tc, dtRecord, dtMax;
	extern double Ta, Tca, Tca_shift, f_rv;
	extern double Tdlrv, Sdlrv;

	// LV parameters
	extern int Nmu, Nnu;
	extern double a0, nu_up, muin0, muout0;
	extern double psi_in_b0, psi_out_b0;
	extern double Ls0, Lsmax, Lsw;
	extern double km, kav, kp, kdfac;
	extern double kv, c1, bff, bfx, bxx;

	// Other chamber parameters
	extern double km_la, kav_la, kv_la, br_la, c1_la, bperp_la, rin0_la, rout0_la;
	extern double km_rv, kav_rv, kv_rv, br_rv, c1_rv, bperp_rv, rin0_rv, rout0_rv;
	extern double km_ra, kav_ra, kv_ra, br_ra, c1_ra, bperp_ra, rin0_ra, rout0_ra;
	extern double At_la;

	// Lumped model parameters
	extern double Vla, Vra, Vrv, Vlv;
	extern double Csa, Csv, Cpa, Cpv, Cpp;
	extern double Rsp, Rsv, Rpa, Rpp, Rpv, Rra;
	extern double q_la, q_lv, q_ra, q_rv;
	extern double q_sp, q_sv, q_pp, q_pv;

	// Valve/vein parameters
	extern double B_lv, B_la, B_rv, B_ra, L_lv, L_la, L_rv, L_ra, L_pv, L_sv, VA_pv, VA_sv, VL_pv, VL_sv;
	extern double Zmiv, Zaov, Ztcv, Zpuv;
	extern double Kaovo, Kaovc, Kmivo, Kmivc, Kpuvo, Kpuvc, Ktcvo, Ktcvc;
	extern double VA_aovo, VA_aovc, VA_mivo, VA_mivc, VA_puvo, VA_puvc, VA_tcvo, VA_tcvc, VL_aov, VL_miv, VL_puv, VL_tcv, rho;

	// Initial pressures, external tissue pressure
	extern double Psp0, Psv0, Ppp0, Ppv0, Pext;

	extern string SimpPulm;

	// used for cubic solving
	extern double third;
	double c0, pp, qq, acube, chout0;
	CX sdel, root;

	double LAtoSA, RefVolLV, RefVolLA, RefVolRV, RefVolRA, WallVolLV, WallVolLA, WallVolRV, WallVolRA, CSS;

	Ncycles = (int)Param[1];
	Ncycles_show = (int)Param[2];
	Nmu = (int)Param[3];
	Nnu = (int)Param[4];
	Nr = (int)Param[5];
	// Time steps
	dtMax = Param[6];
	dtRecord = Param[7];
	// density
	rho = Param[8];
	// Time parameters (s) 
	Tc = Param[9];
	Ta = Param[10];
	Tca = Param[11];
	Tca_shift = Param[12];

	// Reference configuration parameters (LV)
	a0 = Param[13];
	muin0 = Param[14];
	muout0 = Param[15];
	nu_up = Param[16];
	//Alternative geometric parameters
	RefVolLV = Param[17];
	WallVolLV = Param[18];
	LAtoSA = Param[19];
	// Muscle fiber angles
	psi_in_b0 = Param[20];
	psi_out_b0 = Param[21];
	// Muscle fiber sarcomere parameters
	Ls0 = Param[22];
	Lsmax = Param[23];
	Lsw = Param[24];
	// LV active force generation, activation-dependent viscous resistance
	km = Param[25];
	kav = Param[26];
	// LV activation with sarcomere length
	kp = Param[27];
	// LV, RV steepen activation onset (if kdfac < 1)
	kdfac = Param[28];
	// LV passive stress parameters: viscous, elastic
	kv = Param[29];
	c1 = Param[30];
	bff = Param[31];
	bxx = Param[32];
	bfx = Param[33];

	// Right ventricle parameters
	rin0_rv = Param[34];
	rout0_rv = Param[35];
	f_rv = Param[36];
	RefVolRV = Param[37];
	WallVolRV = Param[38];
	km_rv = Param[39];
	kav_rv = Param[40];
	kv_rv = Param[41];
	c1_rv = Param[42];
	br_rv = Param[43];
	bperp_rv = Param[44];

	// Left atrium parameters
	rin0_la = Param[45];
	rout0_la = Param[46];
	RefVolLA = Param[47];
	WallVolLA = Param[48];
	km_la = Param[49];
	kav_la = Param[50];
	kv_la = Param[51];
	c1_la = Param[52];
	br_la = Param[53];
	bperp_la = Param[54];

	// Right atrium parameters
	rin0_ra = Param[55];
	rout0_ra = Param[56];
	RefVolRA = Param[57];
	WallVolRA = Param[58];
	km_ra = Param[59];
	kav_ra = Param[60];
	kv_ra = Param[61];
	c1_ra = Param[62];
	br_ra = Param[63];
	bperp_ra = Param[64];

	// Valve parameters (Mynard dissertation) area VA in cm2, length VL in cm
	VA_aovo = Param[65];
	VL_aov = Param[66];
	VA_mivo = Param[67];
	VL_miv = Param[68];
	VA_puvo = Param[69];
	VL_puv = Param[70];
	VA_tcvo = Param[71];
	VL_tcv = Param[72];
	VA_aovc = Param[73];
	VA_mivc = Param[74];
	VA_puvc = Param[75];
	VA_tcvc = Param[76];
	//Valve rate constants,1/(kPa s) from Mynard in cm2/dyn/s
	Kaovo = Param[77];
	Kmivo = Param[78];
	Kpuvo = Param[79];
	Ktcvo = Param[80];
	Kaovc = Param[81];
	Kmivc = Param[82];
	Kpuvc = Param[83];
	Ktcvc = Param[84];
	// Vein inertia parameters
	VA_pv = Param[85];
	VL_pv = Param[86];
	VA_sv = Param[87];
	VL_sv = Param[88];

	// Systemic resistance
	Rsp = Param[96];
	Rsv = Param[97];
	// Pulmonary resistance
	Rpa = Param[98];
	if (SimpPulm != "on") Rpp = Param[99];
	Rpv = Param[100];
	// Compliances
	Csa = Param[101];
	Csv = Param[102];
	if (SimpPulm == "on") {
		Cpp = Param[103];
	}
	else {
		Cpa = Param[103];
		Cpv = Param[104];
	}
	// External tissue pressure
	Pext = Param[105];
	// Initial pressures
	Psp0 = Param[106];
	Psv0 = Param[107];
	Ppp0 = Param[108];
	Ppv0 = Param[109];
	// Delayed inactivation
	Tdlrv = Param[110];
	Sdlrv = Param[111];

	L_pv = 1.e-4 * rho * VL_pv / VA_pv;
	L_sv = 1.e-4 * rho * VL_sv / VA_sv;

	//Alternative geometric parameters
	if (RefVolLV > 0.) {
		muin0 = atanh((cos(nu_up) + 1.) / 2. / LAtoSA);
		CSS = (cos(nu_up) + 1.) * cosh(muin0) * DSQR(sinh(muin0))
			+ (cosh(muin0) - 1.) * cos(nu_up) * DSQR(sin(nu_up));
		a0 = pow(RefVolLV * 1.5 / PI / CSS, third);
		c0 = cos(nu_up);
		acube = DCUB(a0);
		pp = -(1. - c0 + c0 * c0);
		qq = -c0 * (1. - c0) - 3. * (RefVolLV + WallVolLV) / (2. * PI * acube * (1. + c0));
		sdel = sqrt(DSQR(qq / 2.) + DCUB(pp / 3.));
		root = pow(- qq / 2. + sdel, third) + pow(-qq / 2. - sdel, third);
		chout0 = real(root); // The imaginary part is 0, but still need to take only real part
		muout0 = acosh(chout0);
		//double VolSum = 2. / 3. * PI * acube * ((cos(nu_up) + 1.) * cosh(muout0) * DSQR(sinh(muout0))
		//	+ (cosh(muout0) - 1.) * cos(nu_up) * DSQR(sin(nu_up)));		// should equal RefVolLV + WallVolLV
		Param[13] = a0;
		Param[14] = muin0;
		Param[15] = muout0;
	}
	if (RefVolRV > 0.) {
		rin0_rv = pow(RefVolRV * 0.75 / PI / f_rv, third);
		rout0_rv = pow((RefVolRV + WallVolRV) * 0.75 / PI / f_rv, third);
		Param[34] = rin0_rv;
		Param[35] = rout0_rv;
	}
	if (RefVolLA > 0.) {
		rin0_la = pow(RefVolLA * 0.75 / PI, third);
		rout0_la = pow((RefVolLA + WallVolLA) * 0.75 / PI, third);
		Param[45] = rin0_la;
		Param[46] = rout0_la;
	}
	if (RefVolRA > 0.) {
		rin0_ra = pow(RefVolRA * 0.75 / PI, third);
		rout0_ra = pow((RefVolRA + WallVolRA) * 0.75 / PI, third);
		Param[55] = rin0_ra;
		Param[56] = rout0_ra;
	}
}

