// Main program for AxisymmetricLV7s, for parameter estimation
// Includes options for:
//  - parameter estimation using Levenberg-Marquardt algorithm
//  - uncertainty quanitfication using Monte Carlo Markov chain (MCMC) method
// Includes specification of estimated uncertainty in measured property values
// TWS, June 2026

#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <string>
#include <iomanip>
#include "nrutil.h"

using namespace std;

void readParameters(int run);
void readTargets(int run);
void LevMarq(double* targetVector, double* targetWeight, double* outputVector, int Ntarget,
	double* xx, double* xxlb, double* xxub, int Nopt, double * chisq, double tol, int Niter, int *ia, int run);
void setuparrays(int Nmu, int Nnu, int Nr, int Nmat, int Nvars, int Nt, int Nstore);
void freearrays(int Nmu, int Nnu, int Nr, int Nmat, int Nvars, int Nt, int Nstore);
int aorta(double* Param0, int* ParamVar, int Nparam, double* xx);
void runModel();
void makeGraphs(int run);
void outputFiles(int run);
void updateParams(double* Param, int Nparam);
void updateTargets(double* outputVector, int Ntarget);
void RunMCMC(int run);

int Nvars, Nmat, iRecord, Nstore, Nr, Ncycles, Nt, Nta, Ncycles_show;
double Tc, dtRecord, dtMax, eqDiff;
double Ta, Tca, Tca_shift, f_rv;
double *y, *ytest, *dy;
double *tvec, *a1vec, *a2vec, *a3vec, *a4vec, *a5vec, *a6vec, **Y, **store;

// LV parameters
int Nmu, Nnu;
double a0, nu_up, muin0, muout0;
double psi_in_b0, psi_out_b0;
double Ls0, Lsmax, Lsw;
double km, kav, kp, kdfac;
double kv, c1, bff, bfx, bxx;
double third = 1. / 3., r3 = pow(2., 1. / 3.), chin0;

// Other chambers
double km_la, kav_la, kv_la, br_la, c1_la, bperp_la, rin0_la, rout0_la;
double km_rv, kav_rv, kv_rv, br_rv, c1_rv, bperp_rv, rin0_rv, rout0_rv;
double km_ra, kav_ra, kv_ra, br_ra, c1_ra, bperp_ra, rin0_ra, rout0_ra;
double eta_la, kappa_la, chi_la;
double eta_rv, kappa_rv, chi_rv;
double eta_ra, kappa_ra, chi_ra;
double *r, *kappaj_integrand, *chij_integrand;

double *mu0Vec, *nuVec, **mu0, **nu;
double *psi_eq, *psiStar_eq, **omega;
double **cps, **sps, **cps_sq, **sps_sq;
double **snu, **cnu, **shmu0, **chmu0;
double **snu_sq, **cnu_sq, **sum_sq0, **root_sum_sq0;
double **ell0, **gmu0, **gnu0, **gphi0, **integratingFactor;
double *r0_la, *r0_rv, *r0_ra;

double *X, *f;
double **J, **J1;
double a1, a2, a3, a4, a5, a6;

// Inputs
double alpha_a1, alpha_a2, alpha_a3;
double beta_a1, beta_a2, beta_a3;
double gamma_a1, gamma_a2, gamma_a3;
double kappa_a1, kappa_a2, kappa_a3;
double eta_a1, eta_a2;

// Lumped model variables
double Vla, Vra, Vrv, Vlv;
double Csa, Csv, Cpa, Cpv, Cpp;
double Rsp, Rsv, Rpa, Rpp, Rpv, Rra;
double Ppa, Ppp, Psp, Pra, Ppv, Psa, Prv, Psv, Plv, Pla;
double q_la, q_lv, q_ra, q_rv;
double q_sp, q_sv, q_pp, q_pv;

// Subelement variables
double **shmu, **chmu;
double **sh_sq, **ch_sq, **sum_sq, **root_sum_sq;
double **dmu_dmu0, **dmu_dnu, **dmu_da1, **dmu_da2;
double **ddmu_dnuda1, **ddmu_dnuda2;
double **ddmu_dmu0da1, **ddmu_dmu0da2;

// Deformation gradient tensor (Symmetric)
double **F11, **F12, **F22, **F23, **F33;
double **dF11_da1, **dF11_da2;
double **dF12_da1, **dF12_da2;
double **dF22_da1, **dF22_da2;
double **dF23_da1, **dF23_da2, **dF23_da3;
double **dF33_da1, **dF33_da2;

// Cauchy deformation tensor (Symmetric)
double **C11, **C12, **C22, **C23, **C33;
double **Cinv11, **Cinv12, **Cinv13, **Cinv22, **Cinv23, **Cinv33;

// Green strain (Symmetric)
double **E11, **E12, **E22, **E23, **E33;
double **dE11_da1, **dE11_da2;
double **dE12_da1, **dE12_da2;
double **dE22_da1, **dE22_da2, **dE22_da3;
double **dE23_da1, **dE23_da2, **dE23_da3;
double **dE33_da1, **dE33_da2, **dE33_da3;
double **Ess, **Enn, **Eff, **Esn, **Ens, **Efn, **Enf, **Efs, **Esf;

// Fiber strain
double **ell, **dell_da1, **dell_da2, **dell_da3, **eps_f;
double **eps_fed;
double lambda_mean;

// 2nd Piola-Kirchhoff viscous stress (Not symmetric)
double **Sv11_a1, **Sv12_a1, **Sv13_a1, **Sv21_a1, **Sv22_a1, **Sv23_a1, **Sv31_a1, **Sv32_a1, **Sv33_a1;
double **Sv11_a2, **Sv12_a2, **Sv13_a2, **Sv21_a2, **Sv22_a2, **Sv23_a2, **Sv31_a2, **Sv32_a2, **Sv33_a2;
double **Sv11_a3, **Sv12_a3, **Sv13_a3, **Sv21_a3, **Sv22_a3, **Sv23_a3, **Sv31_a3, **Sv32_a3, **Sv33_a3;

// 2nd Piola-Kirchhoff elastic stress (Symmetric)
double **eW, **Seff, **Sess, **Senn, **Sefn, **Senf, **Sefs, **Sesf, **Sesn, **Sens;
double **Se11, **Se12, **Se13, **Se21, **Se22, **Se23, **Se31, **Se32, **Se33;

// Activation function
double At, At_la;
double Tdlrv, Sdlrv;

// Active fiber stress
double **Sf_ff_const, **Sf_ff_a1, **Sf_ff_a2, **Sf_ff_a3;
double **Sf11_const, **Sf12_const, **Sf13_const, **Sf21_const, **Sf22_const, **Sf23_const, **Sf31_const, **Sf32_const, **Sf33_const;
double **Sf11_a1, **Sf12_a1, **Sf13_a1, **Sf21_a1, **Sf22_a1, **Sf23_a1, **Sf31_a1, **Sf32_a1, **Sf33_a1;
double **Sf11_a2, **Sf12_a2, **Sf13_a2, **Sf21_a2, **Sf22_a2, **Sf23_a2, **Sf31_a2, **Sf32_a2, **Sf33_a2;
double **Sf11_a3, **Sf12_a3, **Sf13_a3, **Sf21_a3, **Sf22_a3, **Sf23_a3, **Sf31_a3, **Sf32_a3, **Sf33_a3;
double Sf_mean, Sf1_mean, Sf2_mean, Sf3_mean, volume;

// Equilibrium integral
double **integrand;

// Volume computation
double *vIntegrand_da1, *vIntegrand_da2, *vol_integrand;

// Valve/vein variables
double B_lv, B_la, B_rv, B_ra, L_lv, L_la, L_rv, L_ra, L_pv, L_sv, VA_sv, VL_sv, VA_pv, VL_pv;
double Zmiv, Zaov, Ztcv, Zpuv;
double Kaovo, Kaovc, Kmivo, Kmivc, Kpuvo, Kpuvc, Ktcvo, Ktcvc;
double VA_aovo, VA_aovc, VA_mivo, VA_mivc, VA_puvo, VA_puvc, VA_tcvo, VA_tcvc, VL_aov, VL_miv, VL_puv, VL_tcv, rho;

//aorta simulation
double Z0start, Z0end, P_echo, Q_trans, P_ao_0, q_lv_0, Vaorta, VELAORTAWAVE;
double *echo, *trans, *q_lv_hist, *p_ao_hist, *tveca;
int Nn;

//External tissue pressure, initial pressures
double Pext, Psp0, Psv0, Ppp0, Ppv0;

//Variables used in optimization
int Nparam, Nopt, Ntarget, NtargetVar, Neval = 0, aortaVar, Niter;
int* ParamVar, * ia;
double tol;
double* Param, * xx, * xxlb, * xxub, * xxbest;
double * outputVector, * outputVector0, * chisq, * deriv_norm, * alamda;
double* Param0, * targetVector, * targetWeight, * targetSigma;
double** derivs;
std::string ParamName[112];	//note: zero-based arrays, include extra element
std::string TargetName[49];
string Optimize, Skip, MCMC, SimpPulm, CalcFiber, SensitivityMatrix;

//variables needed for MCMC
int NstepsMarkov, NburnMarkov, NreportMarkov;
double sigmaMarkov;

int main(int argc, char* argv[]) {
	int iout, Nfit, j;
	int run = 0, runflag = 0;
	char fname[80];
	FILE *ofp;

	do {
		sprintf(fname, "parameters%03i.txt", run);
		//test whether file exists
		ofp = fopen(fname, "r");
		if (ofp) {
			runflag = 1;
			fclose(ofp);
		}
		else {
			runflag = 0;
			printf("*************************\n");
			printf("No more input files found\n");
			printf("*************************\n");
		}
		if(runflag) {
			printf("*******************\n");
			printf("****** %03i ********\n", run);
			printf("*******************\n");	

			readParameters(run);

			if (Skip != "on") {

				updateParams(Param0, Nparam);

				readTargets(run);

				Param = dvector(1, Nparam);		//use for model parameters after adjustment
				outputVector = dvector(1, Ntarget);
				if (Optimize == "on") {
					outputVector0 = dvector(1, Ntarget);
					derivs = dmatrix(1, Ntarget, 1, Nopt);
					chisq = dvector(0, 1); //for first iteration, need to return two values
					alamda = dvector(1, 1);
					deriv_norm = dvector(1, 1);
					derivs = dmatrix(1, Ntarget, 1, Nopt);
				}

				Nta = aorta(Param0, ParamVar, Nparam, xx);
				Nt = 1 + (int)(Tc * Ncycles_show / dtRecord);	// only last cycle if Ncycles_show = 1
				Nstore = 48;
				if (SimpPulm == "on") Nvars = 19;
				else Nvars = 20;
				Nmat = 8;	//this covers all cases
				setuparrays(Nmu, Nnu, Nr, Nmat, Nvars, Nt, Nstore);

				if (Optimize == "on") {
					if (Nopt == 0) printf("*** Error: No optimization variables, cannot optimize\n");
					if (NtargetVar == 0) printf("*** Error: No target variables, cannot optimize\n");
				}
				if (Nopt == 0 || NtargetVar == 0) Optimize = "off";

				if (Optimize == "on") {
					Nfit = 0;
					for (j = 1; j <= Nopt; j++) if (ia[j]) Nfit++;
					printf("Optimizing, %i variable parameters, %i target variables\n", Nfit, NtargetVar);
					for (iout = 1; iout <= Ntarget; iout++) outputVector[iout] = 0.;

					LevMarq(targetVector, targetWeight, outputVector, Ntarget,
						xx, xxlb, xxub, Nopt, chisq, tol, Niter, ia, run);
					makeGraphs(run);
					outputFiles(run);
				}
				else if (MCMC == "on") {

					RunMCMC(run);
				}
				else {
					readTargets(run);
					runModel();
					updateTargets(outputVector, Ntarget);
					makeGraphs(run);
					outputFiles(run);
				}

				freearrays(Nmu, Nnu, Nr, Nmat, Nvars, Nt, Nstore);
				free_dvector(echo, 0, Nta - 1);
				free_dvector(trans, 0, Nta - 1);
				free_dvector(q_lv_hist, 0, Nta - 1);
				free_dvector(p_ao_hist, 0, Nta - 1);
				free_dvector(tveca, 0, Nta - 1);

				free_dvector(Param, 1, Nparam);
				free_dvector(Param0, 1, Nparam);
				free_ivector(ParamVar, 1, Nparam);
				free_dvector(outputVector, 1, Ntarget);

				if (Optimize == "on") {
					free_dvector(alamda, 1, 1);
					free_dvector(deriv_norm, 1, 1);
					free_dvector(xx, 1, Nopt);
					free_dvector(xxlb, 1, Nopt);
					free_dvector(xxub, 1, Nopt);
					free_dvector(outputVector0, 1, Ntarget);
					free_dvector(targetVector, 1, Ntarget);
					free_dvector(targetWeight, 1, Ntarget);
					free_dvector(targetSigma, 1, Ntarget);
					free_dvector(chisq, 0, 1);
					free_ivector(ia, 1, Nopt);
					free_dmatrix(derivs, 1, Ntarget, 1, Nopt);
				}
			}
			else printf("*** File skipped **\n");
		}
		run++;
	} while (runflag);
	return(0);
}
