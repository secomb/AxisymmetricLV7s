// modelDerivs.cpp
// Evaluate output variables, or their derivatives wrt parameters, by finite differences
// according to switch updateDerivs
// TWS, July 2025

#include "nrutil.h"

void varyModel(double* Param0, int* ParamVar, int Nparam,
	double* x, int Nopt, double* outputVector, int Ntarget);

void modelDerivs(double xx[], double* outputVector, double** derivs, int Nopt, int updateDerivs, int *ia)
{
	extern int Nparam, Ntarget, * ParamVar;
	extern double* Param0, * outputVector0;

	int iin, iout;
	double x_diff = 0.01, xin0;

	if (updateDerivs) {//evaluate derivatives using stored outputVector0
		for (iin = 1; iin <= Nopt; iin++) {
			xin0 = xx[iin];
			if (ia[iin]) {
				xx[iin] += x_diff;
				varyModel(Param0, ParamVar, Nparam, xx, Nopt, outputVector, Ntarget);
				xx[iin] -= x_diff;
				for (iout = 1; iout <= Ntarget; iout++) 
					derivs[iout][iin] = (outputVector[iout] - outputVector0[iout]) / x_diff;
			}
			else for (iout = 1; iout <= Ntarget; iout++) derivs[iout][iin] = 0.;
		}
		// restore outputVector values
		for (iout = 1; iout <= Ntarget; iout++) outputVector[iout] = outputVector0[iout];
	}
	else {				//evaluate function, store outputVector0
		varyModel(Param0, ParamVar, Nparam, xx, Nopt, outputVector, Ntarget);
		for (iout = 1; iout <= Ntarget; iout++) outputVector0[iout] = outputVector[iout];
	}
}