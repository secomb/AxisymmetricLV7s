#include "nrutil.h"

void modelDerivs(double xx[], double* outputVector, double** derivs, int Nopt, int updateDerivs, int * ia);

void mrqcof(double y[], double wts[], int ndata, double a[], int *ia, int ma,
	double** alpha, double beta[], double* chisq, double* deriv_norm, double* outputVector, int updateDerivs)
	// Used by mrqmin to evaluate the linearized fitting matrix alpha,
	// and vector beta as in (15.5.8), and calculate χ2.
{
	int i, j, k, l, m, mfit = 0;
	double wt, dy, deriv_j;
	extern double ** derivs;

	if (updateDerivs) {
		for (j = 1; j <= ma; j++) if (ia[j]) mfit++;
		for (j = 1; j <= mfit; j++) {			//Initialize (symmetric) alpha, beta.
			for (k = 1; k <= j; k++) alpha[j][k] = 0.;
			beta[j] = 0.;
		}
		*chisq = 0.;
		modelDerivs(a, outputVector, derivs, ma, updateDerivs, ia);
		for (i = 1; i <= ndata; i++) {				//Summation loop over all data.
			dy = y[i] - outputVector[i];			//deviation from target
			for (j = 0, l = 1; l <= ma; l++) {
				if (ia[l]) {
					wt = derivs[i][l] * wts[i];
					for (j++, k = 0, m = 1; m <= l; m++)
						if (ia[m]) alpha[j][++k] += wt * derivs[i][m];
					beta[j] += dy * wt;
				}
			}
			*chisq += dy * dy * wts[i];		//And find χ2.
		}
		*deriv_norm = 0.;
		for (j = 1; j <= ma; j++) if (ia[j]) {			//Summation over parameters
			deriv_j = 0.;
			for (i = 1; i <= ndata; i++) {				//Summation over targets
				dy = y[i] - outputVector[i];			//deviation from target
				deriv_j += 2. * dy * derivs[i][j] * wts[i];		//derivative of objective function wrt parameter j
			}
			*deriv_norm += DSQR(deriv_j);
		}
		for (j = 2; j <= mfit; j++)			//Fill in the symmetric side.
			for (k = 1; k < j; k++) alpha[k][j] = alpha[j][k];
	}
	else {
		*chisq = 0.;
		modelDerivs(a, outputVector, derivs, ma, updateDerivs, ia);
		for (i = 1; i <= ndata; i++) {				//Summation loop over all data.
			dy = y[i] - outputVector[i];			//deviation from target
			*chisq += dy * dy * wts[i];		//And find χ2.
		}
	}
}