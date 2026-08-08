// RunMCMC
// Sample parameter space using Markov Chain Monte Carlo method
// Use for estimation of uncertainty in parameters
// Based on code by H. Paul Keeler, 2024.
// TWS, March 2026

#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include "nrutil.h"

void updateParams(double* Param, int Nparam);
void updateTargets(double* outputVector, int Ntarget);
void runModel();
void histogram(double** p_tHist, double* xx_lb, double* xx_ub,
    int* ia, int Nopt, int NstepsMarkov, int NburnMarkov, const char* fname1);

double residual(double* x_input)
{
    // returns the residual (sum of squared deviations) at a single point inside a simulation window
    extern int* ParamVar, * ia, Nopt, Ntarget, Nparam;
    extern double* Param, * Param0, * xx, * xxlb, * xxub, * outputVector;
    extern double* targetVector, * targetWeight;

    int pv, i, j, k, ifit;
    double Error, Resid;
    
    ifit = 0;
    for (k = 1; k <= Nopt; k++) if (ia[k]) {
        xx[k] = xxlb[k] + (xxub[k] - xxlb[k]) * x_input[ifit];
        ifit++;
    }
    for (i = 1; i <= Nparam; i++) {
        Param[i] = Param0[i];
        pv = ParamVar[i];
        if (pv) Param[i] *= xx[pv];
    }
    // Compute residual with current parameters
    updateParams(Param, Nparam);
    runModel();
    updateTargets(outputVector, Ntarget);
    Resid = 0.;
    for (j = 1; j <= Ntarget; j++) if (targetWeight[j] > 0.) {  //Summation loop over data
        Error = targetVector[j] - outputVector[j];			//deviation from target
        Resid += DSQR(Error) * targetWeight[j];
    }    
    return Resid;
}

double* normrand(double* randValues, int numbRand, double mu, double sigma)
{
    // simulate pairs of iid normal variables using Box-Muller transform
    int i = 0;
    double x, y, r, rho;
    while (i < numbRand) {
        do {
            x = 2. * rand() / float(RAND_MAX) - 1.;
            y = 2. * rand() / float(RAND_MAX) - 1.;
            r = x * x + y * y;
        } while (r >= 1.);
        rho = sqrt(-2. * log(r) / r);
        randValues[i] = sigma * rho * x + mu; // assign first of random variable pair
        i++;
        if (i < numbRand) { // if more variables are needed, generate second value of random pair
            randValues[i] = sigma * rho * y + mu; // assign second of random variable pair
            i++;
        }
        else break;
    }
    return randValues;
}

double* unirand(double* randValues, int numbRand)
{
    // simulate numbRand uniform random variables on the unit interval
    // storing them in randValues which must be allocated by the caller
    int i;
    for (i = 0; i < numbRand; i++) randValues[i] = (double)rand() / RAND_MAX;
    return randValues;
}

double mean_var(double* set_sample, int ns, double* varX)
{
    // mean and variance of set_sample
    int i;
    double meanX = 0, meanXSquared = 0, tempX;  // initialize statistics variables (for testing results)

    for (i = 1; i <= ns; i++) {
        tempX = set_sample[i];
        meanX += tempX / ((double)ns);
        meanXSquared += tempX * tempX / ((double)ns);
    }
    *varX = meanXSquared - meanX * meanX;
    return meanX;
}

void RunMCMC(int run)
{
    extern int* ia, Nopt, NstepsMarkov, NburnMarkov, NreportMarkov;
    extern double* xx, * xxlb, * xxub, * targetWeight, sigmaMarkov;

    int i, j, j1, k, Nfit = 0, ifit, count, count1, trial = 0, accept = 0;
    double* tRandProposal, * p_numbNormT, ** p_tHist, * p_tPrev, * xx_vals;
    double pdfProposal;      // density for proposed position
    double pdfCurrent;       // density for current position
    double ratioAccept;      // ratio of densities (i.e. acceptance probability)
    double uRand;            // uniform variable for Bernoulli trial
    double meanTemp = 0., varTemp = 0., stdTemp = 0.;
    double Resid0, Resid;     // Initial residual, current residual

    bool booleInsideWindow;
    char fname[80], fname1[80];
    FILE* ofp;

    // initialize (pseudo)-random number generator
    //time_t timeCPU; // use CPU time for seed
    //srand((int)time(&timeCPU));
    srand(42); //to reproduce results

    // Nopt = max value of ipar_ parameters
    // Nfit = actual number of optimized parameters
    for (j = 1; j <= Nopt; j++) if (ia[j]) Nfit++;

    tRandProposal = dvector(0, Nfit - 1);    // proposal for a new position in the random walk
    p_numbNormT = dvector(0, 0);    // random step (normally distributed)
    p_tHist = dmatrix(0, Nfit - 1, 0, NstepsMarkov - NburnMarkov - 1);    // history, 0-based
    p_tPrev = dvector(0, Nfit - 1);    // previous, 0-based
    xx_vals = dvector(1, NstepsMarkov - NburnMarkov);       // for statistics, 1-based

    ///// initial state /////
    ifit = 0;
    for (k = 1; k <= Nopt; k++) if (ia[k]) {    // not all Nopt ipar_ parameters are varied
        tRandProposal[ifit] = 0.5;              // initial reference state: middle of interval
        p_tPrev[ifit] = tRandProposal[ifit]; // store in array
        ifit++;
    }
    Resid0 = residual(tRandProposal);
    pdfCurrent = 1.;    // Normalize pdf by initial pdf

    sprintf(fname, "MCMC_Summary%03i.txt", run);
    ofp = fopen(fname, "w");
    fprintf(ofp, "MCMC estimation of parameter distributions\n");
    fprintf(ofp, "Nsteps = %i, Nburn = %i, proposal sigma = %g\n", NstepsMarkov, NburnMarkov, sigmaMarkov);
    fprintf(ofp, "***** Statistics exclude first Nburn values *****\n");
    fprintf(ofp, "*************************************************\n");
    fclose(ofp);

    ///// start Markov chain /////
    for (j = 0; j < NstepsMarkov; j++) {
        count1 = 0;
        do {                // make sure proposal is acceptable 
            printf("%i ", count1);
            count = 0;
            do {                            // make sure proposal is in bounds
                ifit = 0;
                booleInsideWindow = true;
                for (k = 1; k <= Nopt; k++) if (ia[k]) {
                    (void)normrand(p_numbNormT, 1, 0, sigmaMarkov);   // normally distributed random step
                    tRandProposal[ifit] = p_tPrev[ifit] + p_numbNormT[0];
                    booleInsideWindow = booleInsideWindow && (tRandProposal[ifit] >= 0.) && (tRandProposal[ifit] <= 1.);
                    ifit++;
                }
                count++;
            } while (booleInsideWindow == false && count < 100);
            if (count >= 100) printf("*** Error: no proposal in bounds after 100 trials\n");
            Resid = residual(tRandProposal);
            pdfProposal = exp(-(Resid - Resid0) / 2.);  // Normalize pdf by initial pdf 
            (void)unirand(&uRand, 1);
            ratioAccept = pdfProposal / pdfCurrent;
            count1++;
            trial++;
        } while (uRand >= ratioAccept && count1 < 100);
        if (count1 >= 100) printf("*** Error: no acceptable proposal after 100 trials\n");
        accept++;
        ifit = 0;
        for (k = 1; k <= Nopt; k++) if (ia[k]) {
            p_tPrev[ifit] = tRandProposal[ifit];
            ifit++;
        }
        if (j >= NburnMarkov) { // Save for statistical analysis
            j1 = j - NburnMarkov;
            ifit = 0;
            for (k = 1; k <= Nopt; k++) if (ia[k]) {
                p_tHist[ifit][j1] = tRandProposal[ifit]; // store in 0-based array
                ifit++;
            }
        }
        pdfCurrent = pdfProposal; 
        if (j + 1 > NburnMarkov && ((j + 1) % NreportMarkov == 0 || j + 1 == NstepsMarkov)) {
            // generate statistics report to this point (excluding first NburnMarkov values)
            j1 = j - NburnMarkov;
            ofp = fopen(fname, "a");
            fprintf(ofp, "Steps = %i, Acceptance ratio = %f\n", j + 1, (float)accept / (float)trial);
            ifit = 0; 
            for (k = 1; k <= Nopt; k++) if (ia[k]) {
                for (i = 0; i <= j1; i++) {
                    xx_vals[i + 1] = xxlb[k] + (xxub[k] - xxlb[k]) * p_tHist[ifit][i];
                }
                meanTemp = mean_var(xx_vals, j1, &varTemp);
                stdTemp = sqrt(varTemp);
                fprintf(ofp, "variable %i, mean = %f, s.d. = %f\n", k, meanTemp, stdTemp);
                ifit++;
            }
            fclose(ofp);
        }
    }
    ///// end markov chain /////
    // generate histograms
    sprintf(fname1, "Histogram_xx%03i.txt", run);
    histogram(p_tHist, xxlb, xxub, ia, Nopt, NstepsMarkov, NburnMarkov, fname1);

    free_dvector(p_numbNormT, 0, 0);
    free_dvector(tRandProposal, 0, Nfit - 1);
    free_dmatrix(p_tHist, 0, Nfit - 1, 0, NstepsMarkov - NburnMarkov - 1);
    free_dvector(p_tPrev, 0, Nfit - 1);
    free_dvector(xx_vals, 1, NstepsMarkov - NburnMarkov);
}
