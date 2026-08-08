// freearrays.cpp
// TWS, February 2020

#include <stdio.h>
#include <math.h>
#include "nrutil.h"

void freearrays(int Nmu, int Nnu, int Nr, int Nmat, int Nvars, int Nt, int Nstore) {
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

	free_dvector(dy, 0, Nvars - 1);
	free_dvector(y, 0, Nvars - 1);
	free_dvector(ytest, 0, Nvars - 1);
	
	free_dvector(X, 0,Nmat-1);
	free_dvector(f, 0,Nmat-1);
	free_dmatrix(J, 0,Nmat-1,0,Nmat-1);
	free_dmatrix(J1, 0,Nmat-1,0,Nmat-1);

	free_dvector(mu0Vec, 0,Nmu-1);
	free_dvector(nuVec, 0,Nnu-1);
	free_dvector(psi_eq, 0,Nmu-1);
	free_dvector(psiStar_eq, 0,Nmu-1);

	free_dmatrix(omega, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(cps, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(sps, 0, Nmu-1,0,Nnu-1);
	free_dmatrix(cps_sq, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(sps_sq, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(mu0, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(nu, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(snu, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(cnu, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(shmu0, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(chmu0, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(snu_sq, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(cnu_sq, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(sum_sq0, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(root_sum_sq0, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(ell0, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(gmu0, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(gnu0, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(gphi0, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(integratingFactor, 0,Nmu-1,0,Nnu-1);

	// Initialize subelements
	free_dmatrix(shmu, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(chmu, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(sh_sq, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(ch_sq, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(sum_sq, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(root_sum_sq, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dmu_dmu0, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dmu_dnu, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dmu_da1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dmu_da2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(ddmu_dnuda1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(ddmu_dnuda2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(ddmu_dmu0da1, 0, Nmu-1,0,Nnu-1);
	free_dmatrix(ddmu_dmu0da2, 0,Nmu-1,0,Nnu-1);

	// deformation gradient tensor
	free_dmatrix(F11, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(F12, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(F22, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(F23, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(F33, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dF11_da1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dF11_da2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dF12_da1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dF12_da2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dF22_da1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dF22_da2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dF23_da1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dF23_da2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dF23_da3, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dF33_da1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dF33_da2, 0,Nmu-1,0,Nnu-1);

	// Cauchy deformation tensor
	free_dmatrix(C11, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(C12, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(C22, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(C23, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(C33, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Cinv11, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Cinv12, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Cinv13, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Cinv22, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Cinv23, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Cinv33, 0,Nmu-1,0,Nnu-1);
	
	// Green strain tensor
	free_dmatrix(E11, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(E12, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(E22, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(E23, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(E33, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dE11_da1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dE11_da2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dE12_da1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dE12_da2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dE22_da1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dE22_da2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dE22_da3, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dE23_da1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dE23_da2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dE23_da3, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dE33_da1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dE33_da2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dE33_da3, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Ess, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Enn, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Eff, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Esn, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Ens, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Efn, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Enf, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Efs, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Esf, 0,Nmu-1,0,Nnu-1);

	// Fiber strain
	free_dmatrix(ell, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dell_da1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dell_da2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(dell_da3, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(eps_f, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(eps_fed, 0,Nmu-1,0,Nnu-1);

	// 2nd Piola-Kirchhoff viscous stress
	free_dmatrix(Sv11_a1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sv12_a1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sv13_a1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sv21_a1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sv22_a1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sv23_a1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sv31_a1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sv32_a1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sv33_a1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sv11_a2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sv12_a2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sv13_a2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sv21_a2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sv22_a2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sv23_a2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sv31_a2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sv32_a2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sv33_a2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sv11_a3, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sv12_a3, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sv13_a3, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sv21_a3, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sv22_a3, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sv23_a3, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sv31_a3, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sv32_a3, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sv33_a3, 0,Nmu-1,0,Nnu-1);

	// 2nd Piola-Kirchhoff elastic stress 
	free_dmatrix(eW, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Se11, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Se12, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Se13, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Se21, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Se22, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Se23, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Se31, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Se32, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Se33, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sess, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Senn, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Seff, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sesn, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sens, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sefn, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Senf, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sefs, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sesf, 0,Nmu-1,0,Nnu-1);	

	// active fiber strain vectors
	free_dmatrix(Sf_ff_const, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf_ff_a1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf_ff_a2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf_ff_a3, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf11_const, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf12_const, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf13_const, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf21_const, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf22_const, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf23_const, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf31_const, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf32_const, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf33_const, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf11_a1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf12_a1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf13_a1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf21_a1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf22_a1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf23_a1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf31_a1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf32_a1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf33_a1, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf11_a2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf12_a2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf13_a2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf21_a2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf22_a2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf23_a2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf31_a2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf32_a2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf33_a2, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf11_a3, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf12_a3, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf13_a3, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf21_a3, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf22_a3, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf23_a3, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf31_a3, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf32_a3, 0,Nmu-1,0,Nnu-1);
	free_dmatrix(Sf33_a3, 0,Nmu-1,0,Nnu-1);
	
	// Initialize equilibrium integrands
	free_dmatrix(integrand, 0,Nmu-1,0,Nnu-1);

	// Initialize volume derivative
	free_dvector(vIntegrand_da1, 0,Nnu-1);
	free_dvector(vIntegrand_da2, 0,Nnu-1);
	free_dvector(vol_integrand, 0,Nnu-1);

	// Initialize spherical chamber variables
	free_dvector(r, 0, Nr - 1);
	free_dvector(r0_la, 0, Nr - 1);
	free_dvector(r0_rv, 0, Nr - 1);
	free_dvector(r0_ra, 0, Nr - 1);
	free_dvector(kappaj_integrand, 0,Nr-1);
	free_dvector(chij_integrand, 0,Nr-1);

	// Set up solution storage vector
	free_dmatrix(Y, 0, Nvars - 1, 0, Nt - 1);
	// Set up storage for other variables
	free_dmatrix(store, 0,Nstore-1,0,Nt-1);	

	// Set up storage for a1,a2,a3,a4,a5,a6
	free_dvector(a1vec, 0, Nt - 1);
	free_dvector(a2vec, 0, Nt - 1);
	free_dvector(a3vec, 0, Nt - 1);
	free_dvector(a4vec, 0, Nt - 1);
	free_dvector(a5vec, 0, Nt - 1);
	free_dvector(a6vec, 0, Nt - 1);

	// Set up storage for the time variable
	free_dvector(tvec, 0, Nt - 1);
}
