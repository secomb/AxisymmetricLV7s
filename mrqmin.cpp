// mrqmin.cpp - from Numerical Recipes, 2nd ed.
// modified for function output depend on an index instead of a continous variable
// Also modified to include restriction by upper and lower bounds, following code 
// from MATLAB "projectBox":
/*
trialX = x(:) + dx(:);
activeLB = isfinite(lb) & trialX < lb;
activeUB = isfinite(ub) & trialX > ub;
dx(activeLB) = lb(activeLB) - x(activeLB);
dx(activeUB) = ub(activeUB) - x(activeUB);
*/
#include <stdio.h>
#include <math.h>
#include "nrutil.h"

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

void gaussj1(double** a, int n, double* b);
void mrqcof(double y[], double sig[], int ndata, double a[], int *ia,
	int ma, double** alpha, double beta[], double* chisq, double* deriv_norm, double* outputVector, int updateDerivs);

void covsrt(double** covar, int ma, int ia[], int mfit)
// Expand in storage the covariance matrix covar, so as to take into account parameters that are
// being held fixed. (For the latter, return zero covariances.)
{
	int i, j, k;
	double swap;
	for (i = mfit + 1; i <= ma; i++)
		for (j = 1; j <= i; j++) covar[i][j] = covar[j][i] = 0.;
	k = mfit;
	for (j = ma; j >= 1; j--) {
		if (ia[j]) {
			for (i = 1; i <= ma; i++) SWAP(covar[i][k], covar[i][j])
				for (i = 1; i <= ma; i++) SWAP(covar[k][i], covar[j][i])
					k--;
		}
	}
}

int mrqmin(double y[], double wts[], int ndata, double a[], double alb[], double aub[], int *ia, int ma,
	double** covar, double** alpha, double * atry, double * beta, double * da, double * oneda,
	double* chisq, double* alamda, double* deriv_norm, double* outputVector, int iter, double tol)
{
	extern int Ntarget;
	extern double* outputVector0;

	int j, k, l, updateDerivs = 0;
	double diff, lam2;
	static int mfit;
	static double ochisq;


	if (*alamda < 0.) {				//Initialization		
		for (mfit = 0, j = 1; j <= ma; j++) if (ia[j]) mfit++;
		*alamda = 0.001;
		mrqcof(y, wts, ndata, a, ia, ma, alpha, beta, chisq, deriv_norm, outputVector, 0);	//*****
		ochisq = (*chisq);
		printf("iter = 0, chi2 = %g\n", *chisq);
		//for (j = 1; j <= ma; j++) printf(" %8.4f", a[j]);
		//printf("\n");
		mrqcof(y, wts, ndata, a, ia, ma, alpha, beta, chisq, deriv_norm, outputVector, 1);	//*****
		printf("deriv_norm = %g\n", *deriv_norm);

		chisq[1] = ochisq;
		for (j = 1; j <= ma; j++) atry[j] = a[j];
	}
	for (j = 1; j <= mfit; j++) {		//Alter linearized fitting matrix, by augmenting diagonal elements.
		for (k = 1; k <= mfit; k++) covar[j][k] = alpha[j][k];
		covar[j][j] = alpha[j][j] * (1. + (*alamda));
		oneda[j] = beta[j];
	}
	gaussj1(covar, mfit, oneda);		//Matrix solution - 1-based solver
	for (j = 1; j <= mfit; j++) da[j] = oneda[j];
	if (*alamda == 0.) {				//Once converged, evaluate covariance matrix.
		covsrt(covar, ma, ia, mfit);
		covsrt(alpha, ma, ia, mfit);	//Spread out alpha to its full size too.
		updateDerivs = 1;
		return updateDerivs;
	}
	if (*alamda <= 1.) lam2 = 1.;	//modified step size for steepest descent type steps 
	else lam2 = sqrt(*alamda);
	for (j = 0, l = 1; l <= ma; l++)	//Did the trial succeed?
		if (ia[l]) atry[l] = a[l] + da[++j];
	for (l = 1; l <= ma; l++) {			//Impose bounds on parameters
		if (atry[l] < alb[l]) atry[l] = alb[l];
		if (atry[l] > aub[l]) atry[l] = aub[l];
	}

	mrqcof(y, wts, ndata, atry, ia, ma, covar, da, chisq, deriv_norm, outputVector, 0);  //*****
	if (*chisq < ochisq) { 
		diff = ochisq - *chisq;		
		updateDerivs = 1;
		for (j = 1; j <= Ntarget; j++) outputVector0[j] = outputVector[j];
	}
	else {
		diff = 0.;
		updateDerivs = 0;
	}
	printf("iter = %i, updateDerivs = %i, chi2 = %g, lamda = %g, diff =  %10.6f\nxx =",
		iter, updateDerivs, *chisq, *alamda, diff);
	for (j = 1; j <= ma; j++) printf(" %8.4f", atry[j]);
	printf("\n\n");
	if (*chisq < ochisq) {//Success, accept the new solution
		if (diff > tol) {		// otherwise no need to calculate derivatives
			//*alamda *= 0.1;		//more weight to Gauss-Newton 
			*alamda /= 3.;		//more weight to Gauss-Newton 
			ochisq = (*chisq);
			mrqcof(y, wts, ndata, atry, ia, ma, covar, da, chisq, deriv_norm, outputVector, 1);  //*****
			printf("deriv_norm = %g\n", *deriv_norm);
			for (j = 1; j <= mfit; j++) {
				for (k = 1; k <= mfit; k++) alpha[j][k] = covar[j][k];
				beta[j] = da[j];
			}
			for (l = 1; l <= ma; l++) a[l] = atry[l];
		}
		printf("\n");
	}
	else {							//Failure, increase alamda and return.
		//*alamda *= 10.;				//more weight to steepest descent
		*alamda *= 3.;				//more weight to steepest descent
		*chisq = ochisq;
		updateDerivs = 0;
	}
	return updateDerivs;
}