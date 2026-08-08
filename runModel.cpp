// runModel - version for parameter estimation
// TWS, July 2025

#define _CRT_SECURE_NO_WARNINGS
#include <iostream>

void initialgeometry(int Nmu, int Nnu, int Nr);
int ODEsolver();
int aorta(double* Param0, int* ParamVar, int Nparam, double* xx);
void updateParams(double* Param, int Nparam);
void updateTargets(double* outputVector, int Ntarget);

void runModel() {
	extern int Ncycles, Nmu, Nnu, Nr, Neval;

	initialgeometry(Nmu, Nnu, Nr);
	Neval++;
	printf("%i Start ODEsolver with %i cycles ... ", Neval, Ncycles);
	double startTime = std::clock();
	ODEsolver();
	double endTime = std::clock();
	double programTime = (endTime - startTime) / (double)CLOCKS_PER_SEC;
	printf("run time = %f sec\n", programTime);
}

void varyModel(double* Param0, int* ParamVar, int Nparam,
	double* xx, int Nopt, double* outputVector, int Ntarget)
{
	extern int Nta, aortaVar;
	extern double* Param;

	int i, pv;

	// Current parameter vector
	for (i = 1; i <= Nparam; i++) {
		Param[i] = Param0[i];
		pv = ParamVar[i];
		if (pv) Param[i] *= xx[pv];
	}
	updateParams(Param, Nparam);
	//////////////////////////////////////////////////
	if (aortaVar) {
		int Nta1 = aorta(Param0, ParamVar, Nparam, xx);
		if (Nta1 != Nta) printf("*** Error: Nta changed\n");
	}
	runModel();
	//////////////////////////////////////////////////
	updateTargets(outputVector, Ntarget);
}