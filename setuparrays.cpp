// setuparrays.cpp
// TWS, December 2019

#include <stdio.h>
#include <math.h>
#include "nrutil.h"

void setuparrays(int Nmu, int Nnu, int Nr, int Nmat, int Nvars, int Nt, int Nstore) {
	extern double *mu0Vec, *nuVec, **mu0, **nu;
	extern double *psi_eq, *psiStar_eq, **omega;
	extern double **cps, **sps, **cps_sq, **sps_sq; // Cos(psi^\star), Sin(psi^\star)
	extern double **snu, **cnu, **shmu0, **chmu0;
	extern double **snu_sq, **cnu_sq, **sum_sq0, **root_sum_sq0;
	extern double **ell0, **gmu0, **gnu0, **gphi0, **integratingFactor;
	extern double *r0_la, *r0_rv, *r0_ra;

	extern double *X, *f;
	extern double **J, **J1;
	extern double *y, *ytest;

	//storage of results
	extern double *tvec, *a1vec, *a2vec, *a3vec, *a4vec, *a5vec, *a6vec, **Y, **store;

	extern double **shmu, **chmu;
	extern double **sh_sq, **ch_sq, **sum_sq, **root_sum_sq;
	extern double **dmu_dmu0, **dmu_dnu, **dmu_da1, **dmu_da2;
	extern double **ddmu_dnuda1, **ddmu_dnuda2;
	extern double **ddmu_dmu0da1, **ddmu_dmu0da2;

	// Deformation gradient tensor (Symmetric)
	extern double **F11, **F12, **F22, **F23, **F33;
	extern double **dF11_da1, **dF11_da2;
	extern double **dF12_da1, **dF12_da2;
	extern double **dF22_da1, **dF22_da2;
	extern double **dF23_da1, **dF23_da2, **dF23_da3;
	extern double **dF33_da1, **dF33_da2;

	// Cauchy deformation tensor (Symmetric)
	extern double **C11, **C12, **C22, **C23, **C33;
	extern double **Cinv11, **Cinv12, **Cinv13, **Cinv22, **Cinv23, **Cinv33;

	// Green strain declarations (Symmetric)
	extern double **E11, **E12, **E22, **E23, **E33;
	extern double **dE11_da1, **dE11_da2;
	extern double **dE12_da1, **dE12_da2;
	extern double **dE22_da1, **dE22_da2, **dE22_da3;
	extern double **dE23_da1, **dE23_da2, **dE23_da3;
	extern double **dE33_da1, **dE33_da2, **dE33_da3;
	extern double **Ess, **Enn, **Eff, **Esn, **Ens, **Efn, **Enf, **Efs, **Esf;

	// Fiber strain declaration
	extern double **ell, **dell_da1, **dell_da2, **dell_da3, **eps_f, **eps_fed;

	// 2nd Piola-Kirchhoff viscous stress (Not symmetric)
	extern double **Sv11_a1, **Sv12_a1, **Sv13_a1, **Sv21_a1, **Sv22_a1, **Sv23_a1, **Sv31_a1, **Sv32_a1, **Sv33_a1;
	extern double **Sv11_a2, **Sv12_a2, **Sv13_a2, **Sv21_a2, **Sv22_a2, **Sv23_a2, **Sv31_a2, **Sv32_a2, **Sv33_a2;
	extern double **Sv11_a3, **Sv12_a3, **Sv13_a3, **Sv21_a3, **Sv22_a3, **Sv23_a3, **Sv31_a3, **Sv32_a3, **Sv33_a3;

	// 2nd Piola-Kirchhoff elastic stress (Symmetric)
	extern double **eW, **Seff, **Sess, **Senn, **Sefn, **Senf, **Sefs, **Sesf, **Sesn, **Sens;
	extern double **Se11, **Se12, **Se13, **Se21, **Se22, **Se23, **Se31, **Se32, **Se33;

	// Active fiber stress 
	extern double **Sf_ff_const, **Sf_ff_a1, **Sf_ff_a2, **Sf_ff_a3;
	extern double **Sf11_const, **Sf12_const, **Sf13_const, **Sf21_const, **Sf22_const, **Sf23_const, **Sf31_const, **Sf32_const, **Sf33_const;
	extern double **Sf11_a1, **Sf12_a1, **Sf13_a1, **Sf21_a1, **Sf22_a1, **Sf23_a1, **Sf31_a1, **Sf32_a1, **Sf33_a1;
	extern double **Sf11_a2, **Sf12_a2, **Sf13_a2, **Sf21_a2, **Sf22_a2, **Sf23_a2, **Sf31_a2, **Sf32_a2, **Sf33_a2;
	extern double **Sf11_a3, **Sf12_a3, **Sf13_a3, **Sf21_a3, **Sf22_a3, **Sf23_a3, **Sf31_a3, **Sf32_a3, **Sf33_a3;

	// Spherical chamber variables
	extern double *r, *kappaj_integrand, *chij_integrand;

	// System variables
	extern double *dy;

	// Equilibrium integral declarations
	extern double **integrand;

	// Volume computation declarations
	extern double *vol_integrand, *vIntegrand_da1, *vIntegrand_da2;
/////////////////////////////////////////////////////////////////////////////////////////

	// Initialize system variables
	dy = dvector(0, Nvars - 1);
	// Set up solution vectors
	y = dvector(0, Nvars - 1);
	ytest = dvector(0, Nvars - 1);

	// Initialize iteration vectors
	X = dvector(0,Nmat-1);
	f = dvector(0,Nmat-1);

	// Initialize Jacobian matrix
	J = dmatrix(0,Nmat-1,0,Nmat-1);
	J1 = dmatrix(0,Nmat-1,0,Nmat-1);

	mu0Vec = dvector(0,Nmu-1);
	nuVec = dvector(0,Nnu-1);
	psi_eq = dvector(0,Nmu-1);
	psiStar_eq = dvector(0,Nmu-1);

	omega = dmatrix(0,Nmu-1,0,Nnu-1);
	cps = dmatrix(0,Nmu-1,0,Nnu-1);
	sps = dmatrix(0,Nmu-1,0,Nnu-1);
	cps_sq = dmatrix(0,Nmu-1,0,Nnu-1);
	sps_sq = dmatrix(0,Nmu-1,0,Nnu-1);
	mu0 = dmatrix(0,Nmu-1,0,Nnu-1);
	nu = dmatrix(0,Nmu-1,0,Nnu-1);
	snu = dmatrix(0,Nmu-1,0,Nnu-1);
	cnu = dmatrix(0,Nmu-1,0,Nnu-1);
	shmu0 = dmatrix(0,Nmu-1,0,Nnu-1);
	chmu0 = dmatrix(0,Nmu-1,0,Nnu-1);
	snu_sq = dmatrix(0,Nmu-1,0,Nnu-1);
	cnu_sq = dmatrix(0,Nmu-1,0,Nnu-1);
	sum_sq0 = dmatrix(0,Nmu-1,0,Nnu-1);
	root_sum_sq0 = dmatrix(0,Nmu-1,0,Nnu-1);
	ell0 = dmatrix(0,Nmu-1,0,Nnu-1);
	gmu0 = dmatrix(0,Nmu-1,0,Nnu-1);
	gnu0 = dmatrix(0,Nmu-1,0,Nnu-1);
	gphi0 = dmatrix(0,Nmu-1,0,Nnu-1);
	integratingFactor = dmatrix(0,Nmu-1,0,Nnu-1);

	// Initialize subelements
	shmu = dmatrix(0,Nmu-1,0,Nnu-1);
	chmu = dmatrix(0,Nmu-1,0,Nnu-1);
	for (int i = 0; i < Nmu; i++) for (int j = 0; j < Nnu; j++) chmu[i][j] = 0.;
	sh_sq = dmatrix(0,Nmu-1,0,Nnu-1);
	ch_sq = dmatrix(0,Nmu-1,0,Nnu-1);
	sum_sq = dmatrix(0,Nmu-1,0,Nnu-1);
	root_sum_sq = dmatrix(0,Nmu-1,0,Nnu-1);
	dmu_dmu0 = dmatrix(0,Nmu-1,0,Nnu-1);
	dmu_dnu = dmatrix(0,Nmu-1,0,Nnu-1);
	dmu_da1 = dmatrix(0,Nmu-1,0,Nnu-1);
	dmu_da2 = dmatrix(0,Nmu-1,0,Nnu-1);
	ddmu_dnuda1 = dmatrix(0,Nmu-1,0,Nnu-1);
	ddmu_dnuda2 = dmatrix(0,Nmu-1,0,Nnu-1);
	ddmu_dmu0da1 = dmatrix(0,Nmu-1,0,Nnu-1);
	ddmu_dmu0da2 = dmatrix(0,Nmu-1,0,Nnu-1);

	// deformation gradient tensor
	F11 = dmatrix(0,Nmu-1,0,Nnu-1);
	F12 = dmatrix(0,Nmu-1,0,Nnu-1);
	F22 = dmatrix(0,Nmu-1,0,Nnu-1);
	F23 = dmatrix(0,Nmu-1,0,Nnu-1);
	F33 = dmatrix(0,Nmu-1,0,Nnu-1);
	dF11_da1 = dmatrix(0,Nmu-1,0,Nnu-1);
	dF11_da2 = dmatrix(0,Nmu-1,0,Nnu-1);
	dF12_da1 = dmatrix(0,Nmu-1,0,Nnu-1);
	dF12_da2 = dmatrix(0,Nmu-1,0,Nnu-1);
	dF22_da1 = dmatrix(0,Nmu-1,0,Nnu-1);
	dF22_da2 = dmatrix(0,Nmu-1,0,Nnu-1);
	dF23_da1 = dmatrix(0,Nmu-1,0,Nnu-1);
	dF23_da2 = dmatrix(0,Nmu-1,0,Nnu-1);
	dF23_da3 = dmatrix(0,Nmu-1,0,Nnu-1);
	dF33_da1 = dmatrix(0,Nmu-1,0,Nnu-1);
	dF33_da2 = dmatrix(0,Nmu-1,0,Nnu-1);

	// Cauchy deformation tensor
	C11 = dmatrix(0,Nmu-1,0,Nnu-1);
	C12 = dmatrix(0,Nmu-1,0,Nnu-1);
	C22 = dmatrix(0,Nmu-1,0,Nnu-1);
	C23 = dmatrix(0,Nmu-1,0,Nnu-1);
	C33 = dmatrix(0,Nmu-1,0,Nnu-1);
	Cinv11 = dmatrix(0,Nmu-1,0,Nnu-1);
	Cinv12 = dmatrix(0,Nmu-1,0,Nnu-1);
	Cinv13 = dmatrix(0,Nmu-1,0,Nnu-1);
	Cinv22 = dmatrix(0,Nmu-1,0,Nnu-1);
	Cinv23 = dmatrix(0,Nmu-1,0,Nnu-1);
	Cinv33 = dmatrix(0,Nmu-1,0,Nnu-1);
	
	// Green strain tensor
	E11 = dmatrix(0,Nmu-1,0,Nnu-1);
	E12 = dmatrix(0,Nmu-1,0,Nnu-1);
	E22 = dmatrix(0,Nmu-1,0,Nnu-1);
	E23 = dmatrix(0,Nmu-1,0,Nnu-1);
	E33 = dmatrix(0,Nmu-1,0,Nnu-1);
	dE11_da1 = dmatrix(0,Nmu-1,0,Nnu-1);
	dE11_da2 = dmatrix(0,Nmu-1,0,Nnu-1);
	dE12_da1 = dmatrix(0,Nmu-1,0,Nnu-1);
	dE12_da2 = dmatrix(0,Nmu-1,0,Nnu-1);
	dE22_da1 = dmatrix(0,Nmu-1,0,Nnu-1);
	dE22_da2 = dmatrix(0,Nmu-1,0,Nnu-1);
	dE22_da3 = dmatrix(0,Nmu-1,0,Nnu-1);
	dE23_da1 = dmatrix(0,Nmu-1,0,Nnu-1);
	dE23_da2 = dmatrix(0,Nmu-1,0,Nnu-1);
	dE23_da3 = dmatrix(0,Nmu-1,0,Nnu-1);
	dE33_da1 = dmatrix(0,Nmu-1,0,Nnu-1);
	dE33_da2 = dmatrix(0,Nmu-1,0,Nnu-1);
	dE33_da3 = dmatrix(0,Nmu-1,0,Nnu-1);
	Ess = dmatrix(0,Nmu-1,0,Nnu-1);
	Enn = dmatrix(0,Nmu-1,0,Nnu-1);
	Eff = dmatrix(0,Nmu-1,0,Nnu-1);
	Esn = dmatrix(0,Nmu-1,0,Nnu-1);
	Ens = dmatrix(0,Nmu-1,0,Nnu-1);
	Efn = dmatrix(0,Nmu-1,0,Nnu-1);
	Enf = dmatrix(0,Nmu-1,0,Nnu-1);
	Efs = dmatrix(0,Nmu-1,0,Nnu-1);
	Esf = dmatrix(0,Nmu-1,0,Nnu-1);	

	// Fiber strain
	ell = dmatrix(0,Nmu-1,0,Nnu-1);
	dell_da1 = dmatrix(0,Nmu-1,0,Nnu-1);
	dell_da2 = dmatrix(0,Nmu-1,0,Nnu-1);
	dell_da3 = dmatrix(0,Nmu-1,0,Nnu-1);
	eps_f = dmatrix(0,Nmu-1,0,Nnu-1);
	eps_fed = dmatrix(0,Nmu-1,0,Nnu-1);

	// 2nd Piola-Kirchhoff viscous stress
	Sv11_a1 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sv12_a1 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sv13_a1 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sv21_a1 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sv22_a1 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sv23_a1 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sv31_a1 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sv32_a1 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sv33_a1 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sv11_a2 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sv12_a2 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sv13_a2 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sv21_a2 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sv22_a2 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sv23_a2 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sv31_a2 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sv32_a2 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sv33_a2 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sv11_a3 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sv12_a3 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sv13_a3 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sv21_a3 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sv22_a3 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sv23_a3 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sv31_a3 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sv32_a3 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sv33_a3 = dmatrix(0,Nmu-1,0,Nnu-1);

	// 2nd Piola-Kirchhoff elastic stress 
	eW = dmatrix(0,Nmu-1,0,Nnu-1);
	Se11 = dmatrix(0,Nmu-1,0,Nnu-1);
	Se12 = dmatrix(0,Nmu-1,0,Nnu-1);
	Se13 = dmatrix(0,Nmu-1,0,Nnu-1);
	Se21 = dmatrix(0,Nmu-1,0,Nnu-1);
	Se22 = dmatrix(0,Nmu-1,0,Nnu-1);
	Se23 = dmatrix(0,Nmu-1,0,Nnu-1);
	Se31 = dmatrix(0,Nmu-1,0,Nnu-1);
	Se32 = dmatrix(0,Nmu-1,0,Nnu-1);
	Se33 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sess = dmatrix(0,Nmu-1,0,Nnu-1);
	Senn = dmatrix(0,Nmu-1,0,Nnu-1);
	Seff = dmatrix(0,Nmu-1,0,Nnu-1);
	Sesn = dmatrix(0,Nmu-1,0,Nnu-1);
	Sens = dmatrix(0,Nmu-1,0,Nnu-1);
	Sefn = dmatrix(0,Nmu-1,0,Nnu-1);
	Senf = dmatrix(0,Nmu-1,0,Nnu-1);
	Sefs = dmatrix(0,Nmu-1,0,Nnu-1);
	Sesf = dmatrix(0,Nmu-1,0,Nnu-1);	

	// active fiber strain vectors
	Sf_ff_const = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf_ff_a1 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf_ff_a2 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf_ff_a3 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf11_const = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf12_const = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf13_const = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf21_const = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf22_const = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf23_const = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf31_const = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf32_const = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf33_const = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf11_a1 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf12_a1 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf13_a1 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf21_a1 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf22_a1 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf23_a1 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf31_a1 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf32_a1 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf33_a1 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf11_a2 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf12_a2 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf13_a2 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf21_a2 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf22_a2 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf23_a2 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf31_a2 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf32_a2 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf33_a2 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf11_a3 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf12_a3 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf13_a3 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf21_a3 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf22_a3 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf23_a3 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf31_a3 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf32_a3 = dmatrix(0,Nmu-1,0,Nnu-1);
	Sf33_a3 = dmatrix(0,Nmu-1,0,Nnu-1);

	// Initialize equilibrium integrands
	integrand = dmatrix(0,Nmu-1,0,Nnu-1);

	// Initialize volume derivative
	vIntegrand_da1 = dvector(0,Nnu-1);
	vIntegrand_da2 = dvector(0,Nnu-1);
	vol_integrand = dvector(0,Nnu-1);

	// Initialize spherical chamber variables
	r = dvector(0, Nr - 1);
	r0_la = dvector(0, Nr - 1);
	r0_rv = dvector(0, Nr - 1);
	r0_ra = dvector(0, Nr - 1);
	kappaj_integrand = dvector(0,Nr-1);
	chij_integrand = dvector(0,Nr-1);

	// Set up solution storage vector
	Y = dmatrix(0, Nvars - 1, 0, Nt - 1);
	// Set up storage for other variables
	store = dmatrix(0,Nstore-1,0,Nt-1);

	// Set up storage for a1,a2,a3,a4,a5,a6
	a1vec = dvector(0, Nt - 1);
	a2vec = dvector(0, Nt - 1);
	a3vec = dvector(0, Nt - 1);
	a4vec = dvector(0, Nt - 1);
	a5vec = dvector(0, Nt - 1);
	a6vec = dvector(0, Nt - 1);

	// Set up storage for the time variable
	tvec = dvector(0, Nt - 1);
}
