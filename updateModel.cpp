// updateModel.cpp
// edited by TWS, November 2019
#include <math.h>
#include <stdio.h>
#include <string>
#include "nrutil.h"

#include <chrono>

using namespace std;

void updateMapping();
void updateDefGradTensor();
void updateCauchyTensor();
void updateGreenStrain();
void updateFiberStrain(double t);
void computeViscous();
void computeElastic();
void updateActivation(double t);
void computeActive();
void updateIntegrals(double t);
void sphericalChamber(double aj, double Atj, double kmj, double kavj, double kvj, double Lswj,
	double Ls0j, double Lsmaxj, double brj, double claj, double bperpj,
	double *r0j, double *etaj, double *kappaj, double *chij, double *vol);
int updateDerivs(double t);

int updateModel(double t) {
	extern double Vrv, Vra, Vla, Plv, Pla, Ppa, Psa, Psp, Psv, Pra, Ppv, Ppp, Prv;
	extern double a1, a2, a3, a4, a5, a6;
	extern double At_la, At_rv, f_rv, Zmiv, Zaov, Ztcv, Zpuv, Ls0, Lsmax, Lsw, At;
	extern double q_la, q_lv, q_ra, q_rv, q_pv, q_sv;
	extern double km_la, kav_la, kv_la, br_la, c1_la, bperp_la;
	extern double km_rv, kav_rv, kv_rv, br_rv, c1_rv, bperp_rv;
	extern double km_ra, kav_ra, kv_ra, br_ra, c1_ra, bperp_ra;
	extern double eta_la, kappa_la, chi_la;
	extern double eta_rv, kappa_rv, chi_rv;
	extern double eta_ra, kappa_ra, chi_ra;
	extern double *r0_la, *r0_rv, *r0_ra;
	extern double *y;
	extern string SimpPulm;

	a1 = y[0];
	a2 = y[1];
	a3 = y[2];
	a4 = y[3];
	a5 = y[4];
	a6 = y[5];
	q_pv = y[6];
	q_sv = y[7];
	q_lv = y[8];
	q_la = y[9];
	q_rv = y[10];
	q_ra = y[11];
	Psp = y[12];
	Psv = y[13];
	Ppp = y[14];
	if (SimpPulm == "on") {
		Zaov = y[15];
		Zmiv = y[16];
		Zpuv = y[17];
		Ztcv = y[18];
	}
	else {
		Ppv = y[15];
		Zaov = y[16];
		Zmiv = y[17];
		Zpuv = y[18];
		Ztcv = y[19];
	}

	// left ventricle
	int timer = 0;
	//if (fabs(t - 0.5) < 0.0005) {
	if (timer){
		auto begin = std::chrono::high_resolution_clock::now();
		auto end = std::chrono::high_resolution_clock::now();
		auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
		printf("Time for null: %.3f ms\n", elapsed.count() * 1e-6);
		begin = end;

		updateMapping();
		end = std::chrono::high_resolution_clock::now();
		elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
		printf("Time for updateMapping(): %.3f ms\n", elapsed.count() * 1e-6);
		begin = end;

		updateDefGradTensor();
		end = std::chrono::high_resolution_clock::now();
		elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
		printf("Time for updateDefGradTensor(): %.3f ms\n", elapsed.count() * 1e-6);
		begin = end;

		updateCauchyTensor();
		end = std::chrono::high_resolution_clock::now();
		elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
		printf("Time for updateCauchyTensor(): %.3f ms\n", elapsed.count() * 1e-6);
		begin = end;

		updateGreenStrain();
		end = std::chrono::high_resolution_clock::now();
		elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
		printf("Time for updateGreenStrain(): %.3f ms\n", elapsed.count() * 1e-6);
		begin = end;

		updateFiberStrain(t);
		end = std::chrono::high_resolution_clock::now();
		elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
		printf("Time for updateFiberStrain(t): %.3f ms\n", elapsed.count() * 1e-6);
		begin = end;

		computeViscous();
		end = std::chrono::high_resolution_clock::now();
		elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
		printf("Time for computeViscous(): %.3f ms\n", elapsed.count() * 1e-6);
		begin = end;

		computeElastic();
		end = std::chrono::high_resolution_clock::now();
		elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
		printf("Time for computeElastic(): %.3f ms\n", elapsed.count() * 1e-6);
		begin = end;

		updateActivation(t);
		end = std::chrono::high_resolution_clock::now();
		elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
		printf("Time for updateActivation(t): %.3f ms\n", elapsed.count() * 1e-6);
		begin = end;

		computeActive();
		end = std::chrono::high_resolution_clock::now();
		elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
		printf("Time for computeActive(): %.3f ms\n", elapsed.count() * 1e-6);
		begin = end;

		updateIntegrals(t);
		end = std::chrono::high_resolution_clock::now();
		elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
		printf("Time for updateIntegrals(t): %.3f ms\n", elapsed.count() * 1e-6);
		begin = end;
		// left atrium
		sphericalChamber(a4, At_la, km_la, kav_la, kv_la, Lsw, Ls0, Lsmax, br_la, c1_la, bperp_la,
			r0_la, &eta_la, &kappa_la, &chi_la, &Vla);
		// right ventricle
		sphericalChamber(a5, At, km_rv, kav_rv, kv_rv, Lsw, Ls0, Lsmax, br_rv, c1_rv, bperp_rv,
			r0_rv, &eta_rv, &kappa_rv, &chi_rv, &Vrv);
		Vrv = Vrv * f_rv;
		// right atrium
		sphericalChamber(a6, At_la, km_ra, kav_ra, kv_ra, Lsw, Ls0, Lsmax, br_ra, c1_ra, bperp_ra,
			r0_ra, &eta_ra, &kappa_ra, &chi_ra, &Vra);
		elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
		printf("Time for sphericalChamber(): %.3f ms\n", elapsed.count() * 1e-6);
		begin = end;

		end = std::chrono::high_resolution_clock::now();
		elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
		printf("Time for null: %.3f ms\n", elapsed.count() * 1e-6);
		begin = end;
	}
	else {
		updateMapping();
		updateDefGradTensor();
		updateCauchyTensor();
		updateGreenStrain();
		updateFiberStrain(t);
		computeViscous();
		computeElastic();
		updateActivation(t);
		computeActive();
		updateIntegrals(t);

		// left atrium
		sphericalChamber(a4, At_la, km_la, kav_la, kv_la, Lsw, Ls0, Lsmax, br_la, c1_la, bperp_la,
			r0_la, &eta_la, &kappa_la, &chi_la, &Vla);
		// right ventricle
		sphericalChamber(a5, At, km_rv, kav_rv, kv_rv, Lsw, Ls0, Lsmax, br_rv, c1_rv, bperp_rv,
			r0_rv, &eta_rv, &kappa_rv, &chi_rv, &Vrv);
		Vrv = Vrv * f_rv;
		// right atrium
		sphericalChamber(a6, At_la, km_ra, kav_ra, kv_ra, Lsw, Ls0, Lsmax, br_ra, c1_ra, bperp_ra,
			r0_ra, &eta_ra, &kappa_ra, &chi_ra, &Vra);
	}
	int error = updateDerivs(t);
	if (error) {
		printf("*** Error: updateDerivs failed\n");
		return 1;
	}
	return 0;
}

