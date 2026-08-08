// sphericalChamber.cpp
// spherical chamber model for LA, RV and RA 
// TWS, December 2019
// Updated February 2022
#include <math.h>
#include <stdio.h>
#include "nrutil.h"

double simpsons(double* f, double dx, int Npoints);

void sphericalChamber(double aj, double Atj, double kmj, double kavj, double kvj, double Lswj,
	double Ls0j, double Lsmaxj, double brj, double claj, double bperpj,
	double *r0j, double *etaj, double *kappaj, double *chij, double *vol)
{
	extern int Nr;
	extern double *r, *chij_integrand, *kappaj_integrand;

	double Sv1j_rr, Sv1j_perp, Sf1j_perp, Sf0_perp, Se_perp, Se_rr, sigperp, dErr_daj, dEperp_daj, Eperp, Err, eW;
	double LS, lambda, G, r03, twoLsw_sq = 2. * DSQR(Lswj);

	for (int k = 0; k < Nr; k++) {
		// Calculate new wall locations
		r03 = pow(r0j[k], 3);
		if (aj > r03) {
			printf("***Error: aj too large\n");
			aj = 0.9 * r03;
		}
		r[k] = pow(r03 - aj, 1./3.);
		lambda = r[k] / r0j[k];
		// Calculate green strain and derivatives
		Err = 0.5 * (1. / DSQR(DSQR(lambda)) - 1.);
		Eperp = 0.5 * (DSQR(lambda) - 1);
		dErr_daj = 2. / (3. * r03 * DSQR(DSQR(lambda)) * DCUB(lambda));
		dEperp_daj = -1. / (3. * r03 * lambda);
		// Calculate PK2 passive elastic stress
		eW = exp(brj * DSQR(Err) + 2. * bperpj * DSQR(Eperp));
		Se_rr = claj * brj * eW * Err;
		Se_perp = claj * bperpj * eW * Eperp;
		// Calculate PK2 passive viscous stress
		Sv1j_rr = kvj * 4. / (3. * r03) * lambda;
		Sv1j_perp = -kvj * 2. / (3. * r03 * DSQR(lambda) * DCUB(lambda));
		// Calculate active Cauchy stress
		// Note - does not include strain-dependent force generation terms
		LS = Ls0j * lambda;
		G = exp(-DSQR(LS - Lsmaxj) / twoLsw_sq);
		sigperp = kmj * Atj * G;
		// Calculate active PK2 stress
		Sf0_perp = sigperp / DSQR(lambda) / 2.;		//factor 1/2 here due to 2D isotropic fiber stress
		// Calculate PK2 active viscous stress
		Sf1j_perp = -kavj * Atj * G / (6. * r03 * DCUB(lambda));
		// Calculate integrands
		kappaj_integrand[k] = (Se_rr * dErr_daj + 2. * (Se_perp + Sf0_perp) * dEperp_daj) * DSQR(r0j[k]) * 4. * PI;
		chij_integrand[k] = (Sv1j_rr * dErr_daj + 2. * (Sv1j_perp + Sf1j_perp) * dEperp_daj) * DSQR(r0j[k]) * 4. * PI;
	}

	double dr0 = r0j[1] - r0j[0];
	r03 = pow(r0j[0], 3);				//corrected July 2025
	*kappaj = simpsons(kappaj_integrand, dr0, Nr);
	*chij = simpsons(chij_integrand, dr0, Nr);
	*etaj = -4./3. * PI;
	*vol = 4./3. * PI * (r03 - aj);
}
