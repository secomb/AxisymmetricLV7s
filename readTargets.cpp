// readTargets.cpp
// Reads targets.txt
// If sigma_rel_OUTPUT is not defined, 0 weight is assumed,
// i.e. OUTPUT is not included in optimization
// TWS, July 2025

#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "nrutil.h"

using namespace std;

double findValue(ifstream& fileID, string key);
string findString(ifstream& fileID, string key);

void readTargets(int run)
{
	
	extern double tol;
	extern double* targetVector, * targetSigma, * targetWeight;
	extern int Ntarget, NtargetVar, Niter;
	extern std::string TargetName[49];	//note: zero-based array, include extra element
	extern string Optimize;

	int i;

	Ntarget = 48;
	targetVector = dvector(1, Ntarget);
	targetWeight = dvector(1, Ntarget);
	targetSigma = dvector(1, Ntarget);

	// Volumes (ml)
	TargetName[1] = "LVEDV";
	TargetName[2] = "LVESV";
	TargetName[3] = "LAEDV";
	TargetName[4] = "LAESV";
	TargetName[5] = "RVEDV";
	TargetName[6] = "RVESV";
	TargetName[7] = "RAEDV";
	TargetName[8] = "RAESV";
	TargetName[9] = "SV";
	TargetName[47] = "LVSV";

	// Lengths (cm)
	TargetName[10] = "LVLAEDL";
	TargetName[11] = "LVLAESL";
	TargetName[12] = "LVSAEDD";
	TargetName[13] = "LVSAESD";
	TargetName[14] = "LVEDWT";
	TargetName[15] = "LVESWT";

	// Maximum flows (ml/s)
	TargetName[16] = "PVMAXFLOW";
	TargetName[17] = "LVOTFLOW";
	TargetName[18] = "RVOTFLOW";

	// Velocities (cm/s)
	TargetName[19] = "EWAVEMAXVEL";
	TargetName[20] = "TEWAVEMAXVEL";
	TargetName[21] = "EPRIME";
	TargetName[22] = "VELAORTAWAVE";

	// Ratios (non-dimensional)
	TargetName[23] = "ETOA";
	TargetName[24] = "TETOA";
	TargetName[25] = "ETOEPRIME";

	// Pressures (mmHg)
	TargetName[26] = "SYSTBP";
	TargetName[27] = "DIASTBP";
	TargetName[28] = "PVMEAN";
	TargetName[29] = "PLVMAX";
	TargetName[30] = "PLVED";
	TargetName[31] = "PRVMAX";
	TargetName[32] = "PRVED";
	TargetName[33] = "PPAMEAN";
	TargetName[34] = "PPAED";
	TargetName[35] = "PRAMEAN";
	TargetName[36] = "PAOMEAN";
	TargetName[37] = "PLAMEAN";
	TargetName[38] = "PAOMAX";
	TargetName[46] = "PAOMIN";
	TargetName[48] = "PSVMEAN";

	// Times (s)
	TargetName[39] = "LVOTSTART_END";
	TargetName[40] = "RVOTSTART_END";

	// Strains (non-dimensional)
	TargetName[41] = "LVGLS";
	TargetName[42] = "LVGCS";
	TargetName[43] = "RVSTRAIN";
	TargetName[44] = "LASTRAIN";
	TargetName[45] = "RASTRAIN";
	
	for (i = 1; i <= Ntarget; i++) {
		targetVector[i] = 0.;
		targetWeight[i] = 0.;
		targetSigma[i] = 0.;
	}
	NtargetVar = 0;

	char fname[80];
	sprintf(fname, "targets%03i.txt", run);
	string InputFile = fname;
	ifstream paramFile(InputFile.c_str(), ios::in);

	tol = findValue(paramFile, "tol");
	if (tol == 0.) tol = 0.0001;	//default
	Niter = (int)findValue(paramFile, "Niter");
	if (Niter == 0) Niter = 20;	//default
	for (i = 1; i <= Ntarget; i++) {
		targetVector[i] = findValue(paramFile, TargetName[i]);
		targetSigma[i] = findValue(paramFile, "sigma_rel_" + TargetName[i]);
		if (targetSigma[i] > 0.) NtargetVar++;
		//ETOA and TETOA ratios
		if (i == 23 || i == 24) targetVector[i] = targetVector[i] / (1. + targetVector[i]);
		// Target weight is 1 / (sigma^2), targetSigma is a relative value
		if (targetVector[i] * targetSigma[i] != 0.)
			targetWeight[i] = 1. / DSQR(targetVector[i] * targetSigma[i]);
	}
}