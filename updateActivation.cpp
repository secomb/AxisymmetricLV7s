// updateActivation.cpp
// edited by TWS, November 2019
#include <math.h>
#include "nrutil.h"
#include <string>
using namespace std;

void updateActivation(double t) {
	extern int Nmu, Nnu;
	extern double At, At_la, Tc, Ta, kdfac, Tca_shift, Tca;
	extern double Tdlrv, Sdlrv;
	
	double truncT = t - Tc + Ta - (int)((t - Tc + Ta )/ Tc) * Tc ;	//truncT = 0 at the start of systole
	double T05;

	// Left and right ventricle activation function
	At = 0.;
	if (truncT < Ta + Tdlrv && truncT >= 0) { // eliminated kd, no need to do multiple calculations. TWS 6/2025
		At = sin(PI * pow(truncT / Ta, kdfac));	// Left ventricle activation function, modified 6/25
		if (Tdlrv > 0.) {	//delayed inactivation
			T05 = Ta * pow((1. - asin(Sdlrv) / PI), 1./kdfac);
			if (truncT > T05) At = Sdlrv * DSQR(1. - (truncT - T05) / (Ta + Tdlrv - T05));
		}
	}

	// Left and right atrium activation function
	At_la = 0.;
	if (fabs(truncT - Tc + Tca_shift) <= Tca / 2.) At_la = cos(PI / Tca * (truncT - Tc + Tca_shift));
	else if (fabs(truncT + Tca_shift) <= Tca / 2.) At_la = cos(PI / Tca * (truncT + Tca_shift));

}
