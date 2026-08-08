// LevMarq - Levenberg-Marquardt method
//  - based on Numerical Recipes, 2nd edition
//  - modified so derivatives are updated only if a successful step is made
//  - modified to implement bounds on parameters
// TWS, July 2025

#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include <string>
#include "nrutil.h"

using namespace std;

int mrqmin(double y[], double wts[], int ndata, double a[], double alb[], double aub[], int *ia, int ma,
	double** covar, double** alpha, double* atry, double* beta, double* da, double* oneda,
	double* chisq, double* alamda, double* deriv_norm, double* outputVector, int iter, double tol);

void LevMarq(double* targetVector, double* targetWeight, double* outputVector, int Ntarget,
	double* xx, double* xxlb, double* xxub, int Nopt, double *chisq, double tol, int Niter, int *ia, int run)
{
	extern double * deriv_norm, * alamda, ** derivs, * xxbest;
	
	extern string SensitivityMatrix;
	
	int iter, Niter0, updateDerivs, n_fail = 0, i;
	double chisqp = 0., diff, chisqbest = 1.e6;
	double * atry, * beta, *da, *oneda;
	double** covar, ** alpha;
	char fname[80];
	FILE* ofp;

	sprintf(fname, "chisqVals%03i.txt", run);
	ofp = fopen(fname, "w");
	fprintf(ofp, "iter  chisq\n");

	covar = dmatrix(1, Nopt, 1, Nopt);
	alpha = dmatrix(1, Nopt, 1, Nopt);
	atry = dvector(1, Nopt);
	beta = dvector(1, Nopt);
	da = dvector(1, Nopt);
	oneda = dvector(1, Nopt);

	*chisq = 0.;
	*alamda = -1.;
	printf("\n");

	Niter0 = Niter;
	if (SensitivityMatrix == "on") {
		Niter0 = 0;		// This bypasses optimization to give initial sensitivities
		printf("Evaluating initial sensitivity matrix, bypassing optimization\n");
		fprintf(ofp, "Evaluating initial sensitivity matrix, bypassing optimization\n");
	}
	iter = 0;
	do {
		updateDerivs = mrqmin(targetVector, targetWeight, Ntarget, xx, xxlb, xxub, ia, Nopt, covar, alpha,
			atry, beta, da, oneda, chisq, alamda, deriv_norm, outputVector, iter, tol);
		iter++;
		if (iter == 1) {
			diff = fabs(chisq[1] - *chisq);
			fprintf(ofp, "0 %g\n", chisq[1]);
		}
		else diff = fabs(chisqp - *chisq);
		if (updateDerivs == 1) {
			fprintf(ofp, "%i %g\n", iter, *chisq);
			n_fail = 0;
		}
		else n_fail++;		//terminate if stuck on updateDerivs = 0
		chisqp = *chisq;
		if (*chisq < chisqbest) {	// keep track of best result so far
			chisqbest = *chisq;
			for (i = 1; i <= Nopt; i++) xxbest[i] = xx[i];
		}
	} while (iter <= Niter0 && (diff > tol || (updateDerivs == 0 && n_fail < 10)));
	// Note: diff = 0 if algorithm did not reduce chisq 

	if (iter >= Niter) {
		printf("*** Warning: Optimization not converged after %i iterations\n", iter-1);
		fprintf(ofp, "*** Warning: Optimization not converged after %i iterations\n", iter-1);
	}
	fclose(ofp);
	
	free_dmatrix(covar, 1, Nopt, 1, Nopt);
	free_dmatrix(alpha, 1, Nopt, 1, Nopt);
	free_dvector(oneda, 1, Nopt);
	free_dvector(da, 1, Nopt);
	free_dvector(beta, 1, Nopt);
	free_dvector(atry, 1, Nopt);
}
