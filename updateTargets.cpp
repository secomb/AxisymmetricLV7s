// updateTargets
// TWS, July 2025

#define _CRT_SECURE_NO_WARNINGS
# include <math.h>
#include "nrutil.h"

void updateTargets(double* outputVector, int Ntarget)
{
	extern int Nt;
	extern double dtRecord, Ta, Tc, Tca, Tca_shift, rin0_rv, rin0_la, rin0_ra;
	extern double muin0, nu_up, VELAORTAWAVE;

	extern double** store, *a4vec, *a5vec, *a6vec;

	///////////////////////////////////////////////////////////

	int j;
	
	// Target variables
	// Time index values in store 
	int NtEC = Nt - 1;										// end cycle
	int NtED = Nt - 2 - (int)(Ta / dtRecord);				// end diastole, -2 because of rounding
	int NtEC1 = Nt - 1 - (int)(Tc / dtRecord);				// end previous cycle
	int NtSAS = Nt - 1 - (int)((Ta + Tca_shift + Tca / 2.) / dtRecord);	// start atrial systole
	int NtESLV = 0;			//end LV systole, to be computed
	int NtESRV = 0;			//end RV systole, to be computed
	
	// Times (s)
	double LVOTSTART_END = 0., RVOTSTART_END = 0., Q0 = 0.1;
	//for (j = NtEC1; j < NtEC; j++) {		//just last cycle
	//	if (store[17][j] > Q0) {
	//		LVOTSTART_END += dtRecord;
	//		NtESLV = j;
	//	}
	//	if (store[19][j] > Q0) {
	//		RVOTSTART_END += dtRecord;
	//		NtESRV = j;
	//	}
	//}
	for (j = NtEC1 + 1; j < NtEC; j++) {		//just last cycle, time estimated using linear interpolation
		if (store[17][j] >= Q0) {
			if (store[17][j - 1] >= Q0) LVOTSTART_END += dtRecord;
			else LVOTSTART_END += dtRecord * (store[17][j] - Q0) / (store[17][j] - store[17][j - 1]);
		}
		else if (store[17][j - 1] >= Q0) {
			LVOTSTART_END += dtRecord * (store[17][j - 1] - Q0) / (store[17][j - 1] - store[17][j]);
			NtESLV = j;
		}
		if (store[19][j] >= Q0) {
			if (store[19][j - 1] >= Q0) RVOTSTART_END += dtRecord;
			else RVOTSTART_END += dtRecord * (store[19][j] - Q0) / (store[19][j] - store[19][j - 1]);
		}
		else if (store[19][j - 1] >= Q0) {
			RVOTSTART_END += dtRecord * (store[19][j - 1] - Q0) / (store[19][j - 1] - store[19][j]);
			NtESRV = j;
		}
	}

	// Volumes (ml)
	double LVEDV = 0., LVESV = 1000., RVEDV = 0., RVESV = 1000.;
	double LAEDV = 0., LAESV = 1000., RAEDV = 0., RAESV = 1000.;
	for (j = NtEC1; j < NtEC; j++) {		//just last cycle
		LVEDV = DMAX(LVEDV, store[1][j]);
		LVESV = DMIN(LVESV, store[1][j]);
		RVEDV = DMAX(RVEDV, store[5][j]);
		RVESV = DMIN(RVESV, store[5][j]);
		LAEDV = DMAX(LAEDV, store[3][j]);
		LAESV = DMIN(LAESV, store[3][j]);
		RAEDV = DMAX(RAEDV, store[7][j]);
		RAESV = DMIN(RAESV, store[7][j]);
	}
	double SV = (LVEDV + RVEDV - LVESV - RVESV) / 2;
	double LVSV = LVEDV - LVESV;

	// Lengths (cm) (internal)
	double LVLAEDL = store[35][NtED];
	double LVLAESL = store[35][NtESLV];
	double LVSAEDD = store[45][NtED];
	double LVSAESD = store[45][NtESLV];
	double LVEDWT = (store[46][NtED] - store[45][NtED]) / 2.;
	double LVESWT = (store[46][NtESLV] - store[45][NtESLV]) / 2.;

	// Maximum flows (ml/s)
	double PVMAXFLOW = 0., LVOTFLOW = 0., RVOTFLOW = 0.;
	for (j = NtEC1; j < NtEC; j++) {		//just last cycle
		PVMAXFLOW = DMAX(PVMAXFLOW, store[25][j]);
		LVOTFLOW = DMAX(LVOTFLOW, store[17][j]);
		RVOTFLOW = DMAX(RVOTFLOW, store[19][j]);
	}

	// Maximum velocities (cm/s)
	double EWAVEMAXVEL = 0., TEWAVEMAXVEL = 0., EPRIME = 0.;
	for (j = NtEC1; j < NtSAS; j++) {		// just until start of atrial systole
		EWAVEMAXVEL = DMAX(EWAVEMAXVEL, store[42][j]);
		TEWAVEMAXVEL = DMAX(TEWAVEMAXVEL, store[44][j]);
		EPRIME = DMAX(EPRIME, -store[36][j]);
	}
	double AWAVEMAXVEL = 0., TAWAVEMAXVEL = 0.;
	for (j = NtSAS; j < NtED; j++) {		// just after start of atrial systole
		AWAVEMAXVEL = DMAX(AWAVEMAXVEL, store[42][j]);
		TAWAVEMAXVEL = DMAX(TAWAVEMAXVEL, store[44][j]);
	}

	// Ratios (non-dimensional)
	double ETOA = 0.;
	if (AWAVEMAXVEL > 0.) ETOA = EWAVEMAXVEL / AWAVEMAXVEL;
	double TETOA = 0.;
	if (TAWAVEMAXVEL > 0.) TETOA = TEWAVEMAXVEL / TAWAVEMAXVEL;
	double ETOEPRIME = 0.;
	if (EPRIME > 0.) ETOEPRIME = EWAVEMAXVEL / EPRIME;

	// Pressures (mmHg)
	double SYSTBP = 0., PLVMAX = 0., PRVMAX = 0., DIASTBP = 1000., PAOMAX = 0., PAOMIN = 1000.;
	for (j = NtEC1; j < NtEC; j++) {		//just last cycle
		SYSTBP = DMAX(SYSTBP, store[11][j]);	// this is Psp 		
		DIASTBP = DMIN(DIASTBP, store[11][j]);
		PLVMAX = DMAX(PLVMAX, store[2][j]);
		PRVMAX = DMAX(PRVMAX, store[6][j]);
		PAOMAX = DMAX(PAOMAX, store[10][j]);	// this is Psa
		PAOMIN = DMIN(PAOMIN, store[10][j]);	// this is Psa
	}
	SYSTBP *= 7.5;
	PLVMAX *= 7.5;
	PRVMAX *= 7.5;
	DIASTBP *= 7.5;
	PAOMAX *= 7.5;
	PAOMIN *= 7.5;
	double PVMEAN = 0., PPAMEAN = 0., PLAMEAN = 0., PRAMEAN = 0., PAOMEAN = 0., PSVMEAN = 0.;
	for (j = NtEC1; j < NtEC; j++) {		//just last cycle
		PVMEAN += store[15][j];
		PPAMEAN += store[13][j];
		PLAMEAN += store[4][j];
		PRAMEAN += store[8][j];
		PAOMEAN += store[39][j];
		PSVMEAN += store[12][j];
	}
	PVMEAN *= 7.5 / (NtEC - NtEC1);
	PPAMEAN *= 7.5 / (NtEC - NtEC1);
	PLAMEAN *= 7.5 / (NtEC - NtEC1);
	PRAMEAN *= 7.5 / (NtEC - NtEC1);
	PAOMEAN *= 7.5 / (NtEC - NtEC1);
	PSVMEAN *= 7.5 / (NtEC - NtEC1);
	double PLVED = store[2][NtED] * 7.5;
	double PRVED = store[6][NtED] * 7.5;
	double PPAED = store[13][NtED] * 7.5;

	// Strains (non-dimensional)
	double LVGLS = (DSQR(LVLAESL / LVLAEDL) - 1.) / 2.;
	double LVGCS = (DSQR(LVSAESD / LVSAEDD) - 1.) / 2.;
	double a4Min = 1000., a5Min = 1000., a6Min = 1000.;
	double a4Max = -1000., a5Max = -1000., a6Max = -1000.;
	for (j = NtEC1; j < NtEC; j++) {		//just last cycle
		a4Min = DMIN(a4Min, a4vec[j]);
		a5Min = DMIN(a5Min, a5vec[j]);
		a6Min = DMIN(a6Min, a6vec[j]);
		a4Max = DMAX(a4Max, a4vec[j]);
		a5Max = DMAX(a5Max, a5vec[j]);
		a6Max = DMAX(a6Max, a6vec[j]);
	}
	double third = 1. / 3.;
	double LASTRAIN = (pow(((DCUB(rin0_la) - a4Max) / (DCUB(rin0_la) - a4Min)), 2. * third) - 1.) / 2.;
	double RVSTRAIN = (pow(((DCUB(rin0_rv) - a5Max) / (DCUB(rin0_rv) - a5Min)), 2. * third) - 1.) / 2.;
	double RASTRAIN = (pow(((DCUB(rin0_ra) - a6Max) / (DCUB(rin0_ra) - a6Min)), 2. * third) - 1.) / 2.;

	// Volumes (ml)
	outputVector[1] = LVEDV;
	outputVector[2] = LVESV;
	outputVector[3] = LAEDV;
	outputVector[4] = LAESV;
	outputVector[5] = RVEDV;
	outputVector[6] = RVESV;
	outputVector[7] = RAEDV;
	outputVector[8] = RAESV;
	outputVector[9] = LVSV;
	
	// Lengths (cm)
	outputVector[10] = LVLAEDL;
	outputVector[11] = LVLAESL;
	outputVector[12] = LVSAEDD;
	outputVector[13] = LVSAESD;
	outputVector[14] = LVEDWT;
	outputVector[15] = LVESWT;
	
	// Maximum flows (ml/s)
	outputVector[16] = PVMAXFLOW;
	outputVector[17] = LVOTFLOW;
	outputVector[18] = RVOTFLOW;
	
	// Maximum velocities (cm/s)
	outputVector[19] = EWAVEMAXVEL;
	outputVector[20] = TEWAVEMAXVEL;
	outputVector[21] = EPRIME;
	outputVector[22] = VELAORTAWAVE;
	
	// Ratios (non-dimensional)
	outputVector[23] = ETOA / (1. + ETOA);
	outputVector[24] = TETOA / (1. + TETOA);
	outputVector[25] = ETOEPRIME;
	
	// Pressures (mmHg)
	outputVector[26] = SYSTBP;
	outputVector[27] = DIASTBP;
	outputVector[28] = PVMEAN;
	outputVector[29] = PLVMAX;
	outputVector[30] = PLVED;
	outputVector[31] = PRVMAX;
	outputVector[32] = PRVED;
	outputVector[33] = PPAMEAN;
	outputVector[34] = PPAED;
	outputVector[35] = PRAMEAN;
	outputVector[36] = PAOMEAN;
	outputVector[37] = PLAMEAN;
	outputVector[38] = PAOMAX;
	outputVector[48] = PSVMEAN;

	// Times (s)
	outputVector[39] = LVOTSTART_END;
	outputVector[40] = RVOTSTART_END;
	
	// Strains (non-dimensional)
	outputVector[41] = LVGLS;
	outputVector[42] = LVGCS;
	outputVector[43] = RVSTRAIN;
	outputVector[44] = LASTRAIN;
	outputVector[45] = RASTRAIN;
	outputVector[46] = PAOMIN;
	outputVector[47] = LVSV;
}
