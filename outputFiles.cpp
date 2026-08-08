// outputFiles.cpp
// TWS, November 2019, July 2025
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>

using namespace std;

void outputFiles(int run) {
	extern int Nt, Nnu, Nr, Nstore, Nparam, Nopt, Ntarget;
	extern int* ParamVar, * ia;
	extern double muin0, muout0, a0, nu_up, Ls0, Lsmax, Lsw, Tc;
	extern double rin0_la, rout0_la, rin0_rv, rout0_rv, rin0_ra, rout0_ra;
	extern double* tvec, * a1vec, * a2vec, * a3vec, * a4vec, * a5vec, * a6vec;
	extern double* targetVector, * targetSigma, * targetWeight, * outputVector, * outputVector0, * chisq;
	extern double* xx, * xxlb, * xxub, * xxbest;
	extern double* Param0, * Param;
	extern double** store, ** Y;
	extern double* deriv_norm, ** derivs;
	extern std::string ParamName[112];	//note: zero-based arrays, include extra element
	extern std::string TargetName[49];
	extern string Optimize;

	int i, k, Nt1 = Nt - 1, iin = 0;
	double result, target, residual, res_sum;
	char fname[80];
	FILE* ofp;
	fstream file1;

	//Write combined text file
	sprintf(fname, "outputFile_All%03i.txt", run);
	fopen_s(&ofp, fname, "w");
	fprintf(ofp, "t a1 a2 a3 a4 a5 a6");
	fprintf(ofp, " Vtot Vlv Plv Vla Pla Vrv Prv Vra Pra null Psa Psp Psv Ppa Ppp Ppv q_la q_lv q_ra q_rv null q_sp q_sv null q_pp q_pv Zmiv Zaov");
	fprintf(ofp, " At     At_la lambda Sfib null Ztcv Zpuv b-a_length dt(length) P_echo Q_trans P_ao_0 Q_lv_0");
	fprintf(ofp, " vel_lv vel_la vel_rv vel_ra LV_diam_int LV_diam_ext Vaorta\n");
	for (i = 0; i < Nt1; i++) {
		fprintf(ofp, "%g %g %g %g %g %g %g", tvec[i], a1vec[i], a2vec[i], a3vec[i], a4vec[i], a5vec[i], a6vec[i]);
		for (k = 0; k < Nstore; k++) fprintf(ofp, " %g", store[k][i]);
		fprintf(ofp, "\n");
	}
	fclose(ofp);

	sprintf(fname, "outputList%03i.txt", run);
	file1.open(fname, ios_base::out);
	if (Optimize == "on") {		//Write list of output parameters, if optimized	
		file1 << "**** Outputs normalized to target values in optimization ****" << endl;
		file1 << "chi2 = " << chisq[0] << endl;
		file1 << "**** Optimized multipliers ****" << endl;
		for (i = 1; i <= Nopt; i++) file1 << " xx[" << i << "] = " << xx[i];
		file1 << endl;
		for (i = 1; i <= Nopt; i++) file1 << " lb[" << i << "] = " << xxlb[i];
		file1 << endl;
		for (i = 1; i <= Nopt; i++) file1 << " ub[" << i << "] = " << xxub[i];
		file1 << endl;
	}
	file1 << "**** Input parameters ****" << endl;
	for (i = 1; i <= Nparam; i++) {
		file1 << ParamName[i] << " = " << Param[i];
		if (ParamVar[i] && Optimize == "on") file1 << " optimized by xx[" << ParamVar[i] << "]";
		file1 << endl;
	}
	file1 << "**** Output variables ****" << endl;

	for (i = 1; i <= Ntarget; i++) {
		target = targetVector[i];
		if (Optimize == "on") outputVector[i] = outputVector0[i];	//use best result, may not be last evaluation
		if (i == 23 || i == 24) result = outputVector[i] / (1. - outputVector[i]);
		else result = outputVector[i];
		file1 << TargetName[i] << " = " << result;
		if (targetSigma[i] > 0. && Optimize == "on") {
			if (i == 23 || i == 24) result = targetVector[i] / (1. - targetVector[i]);
			else result = targetVector[i];
			file1 << ", target = " << result;
			file1 << ", relative sigma = " << targetSigma[i];
		}
		file1 << endl;
	}	
	file1.close();

	sprintf(fname, "residuals%03i.txt", run);
	ofp = fopen(fname, "w");
	res_sum = 0.;
	fprintf(ofp, "     Variable         Target         Fitted       Residual   Target sigma Residual/sigma\n");
	for (i = 1; i <= Ntarget; i++) if (targetWeight[i] > 0.) {
		target = targetVector[i];
		result = outputVector[i];
		residual = result - target;
		res_sum += DSQR(residual) * targetWeight[i];
		residual /= target;
		fprintf(ofp, "%13s  %13g  %13g  %13g  %13g  %13g\n",
			TargetName[i].c_str(), target, result, residual, targetSigma[i], residual / targetSigma[i]);
	}
	fprintf(ofp, "Objective function = %g\n", res_sum);
	printf("Objective function = %g\n", res_sum);
	fclose(ofp);

	if (Optimize == "on") {
		// Normalized sensitivity matrix, correlation matrix
		sprintf(fname, "Sensitivity%03i.txt", run);
		ofp = fopen(fname, "w");
		fprintf(ofp, "Normalized sensitivity matrix\n");
		for (i = 1; i <= Ntarget; i++) if (targetWeight[i] > 0.) {
			for (iin = 1; iin <= Nopt; iin++) if (ia[iin]) fprintf(ofp, "%g ", xx[iin] * derivs[i][iin] / outputVector[i]);
			fprintf(ofp, "\n");
		}
		fclose(ofp);

		//Write list of output parameters, if optimized, using best case values
		for (i = 1; i <= Nparam; i++) {
			Param[i] = Param0[i];
			int pv = ParamVar[i];
			if (pv) Param[i] *= xxbest[pv];
		}
		sprintf(fname, "Newparameters%03i.txt", run);
		file1.open(fname, ios_base::out);
		file1 << "**** Input parameters ****" << endl;
		for (i = 1; i <= Nparam; i++) {
			file1 << ParamName[i] << " = " << std::fixed << std::setprecision(8) << Param[i] << endl;
			if (ParamVar[i]) file1 << "ipar_" << ParamName[i] << " = " << ParamVar[i] << endl;
		}
		file1.close();
	}
}
