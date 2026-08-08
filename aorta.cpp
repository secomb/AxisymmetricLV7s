// aorta.cpp
// TWS, June, 2020
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include <string>
#include "nrutil.h"

using namespace std;

/**********************************************************
Model for tapered aorta using finite differences for x-derivative
Put pressures on the nodes, flows on the half-nodes
Equations:
dt(p) + a0/(A0 gamma) dx(q) = 0
dt(q) + A0/rho dx(p) = 0
gamma = 2 da0/dp
compliance G = (1/A0) dA0/dp = gamma/a0
so dt(p) + 1/(A0 G) dx(q) = 0
wave speed = 1/sqrt(G rho)
Assume constant G and rho, variable A0
Units here are cgs
TWS, May 2020
*************************************************************/

void derivs(double *p, double *q, double *pt, double *qt, double *A0, double *G, double t, int Nn, 
	double dx, double tpulse, double impedancefac)
{
	extern double Z0start, Z0end, rho;
	
	int i;
	double qin, qint, qintt;

	//input is short pulse of flow. Peak is at t = tpulse/2.
	if (t < tpulse) {
		qin = (1 - cos(2. * PI * t / tpulse)) / tpulse;
		qint = 2. * PI / DSQR(tpulse) * sin(2. * PI * t / tpulse);
		qintt = 4. * DSQR(PI / tpulse) / tpulse * cos(2. * PI * t / tpulse);
	}
	else {
		qin = 0.;
		qint = 0.;
		qintt = 0.;
	}
	// Equations
	// dt(p) = - 1 / (A0 G) dx(q)
	// dt(q) = - A0 / rho dx(p)
	for (i = 1; i < Nn; i++) { 		// pressure and flow time derivatives
		if (i == 1) { // flow boundary condition
			pt[i] = -1. / (A0[i] * G[i] * dx) * (q[2] - qin) + rho * dx / A0[i] * qintt;
			qt[i] = qint;
		}
		else {
			if (i == 2) pt[i] = -1. / (A0[i] * G[i] * dx) * (q[i] - qin);
			else pt[i] = -1. / (A0[i] * G[i] * dx) * (q[i] - q[i - 1]);
			qt[i] = -A0[i] / (rho * dx) * (p[i + 1] - p[i]);
		}
	}
	// impedance boundary condition to avoid reflection
	pt[Nn] = -impedancefac * Z0end * A0[Nn] / rho / dx * (p[Nn] - p[Nn - 1]);
}

int aorta(double* Param0, int* ParamVar, int Nparam, double* xx){
	extern double Z0start, Z0end, Vaorta, rho, VELAORTAWAVE;
	extern double *echo, *trans, *q_lv_hist, *p_ao_hist, *tveca, *Param;
	extern string Optimize;

	int i, ii, j, k;
	int Ncycles = 3;				//number of cardiac cycles to include
	int Nn = 501;					// number of nodes
	int Nta, Npulse, Ntdiv = 10;	// very small time steps are needed for stability with explicit method		
	double *x, *p, *q, *A0, *pold, *qold, *c, *G;
	double *pt1, *pt2, *pt3, *pt4, *qt1, *qt2, *qt3, *qt4;
	double t, dt, dtMax, dx, echo_av, Tc, Laorta, A00, A01, Gaorta0, Gaorta1, impedancefac, tpulse, t_transit;

	// Current parameter vector
	for (i = 1; i <= Nparam; i++) {
		Param[i] = Param0[i];
		if (ParamVar[i] && Optimize == "on") Param[i] *= xx[ParamVar[i]];
	}
	// Time parameters (s) 
	dtMax = Param[6];
	Tc = Param[9];	
	//aorta parameters
	Laorta = Param[89];
	A00 = Param[90];
	A01 = Param[91];
	Gaorta0 = Param[92];
	Gaorta1 = Param[93];
	tpulse = Param[94];
	impedancefac = Param[95];

	dx = Laorta / (Nn - 1);
	dt = dtMax / Ntdiv;					// subdivision of time steps for integration
	Nta = (int)(Tc / dtMax + 0.1);		// number of time points output, used for convolution
	Npulse = (int)(tpulse / 2. / dtMax);// number of time points before peak of pulse is reached

	x = dvector(1, Nn);
	p = dvector(1, Nn);
	pold = dvector(1, Nn);
	pt1 = dvector(1, Nn);
	pt2 = dvector(1, Nn);
	pt3 = dvector(1, Nn);
	pt4 = dvector(1, Nn);
	q = dvector(1, Nn - 1);
	qold = dvector(1, Nn - 1);
	qt1 = dvector(1, Nn - 1);
	qt2 = dvector(1, Nn - 1);
	qt3 = dvector(1, Nn - 1);
	qt4 = dvector(1, Nn - 1);
	A0 = dvector(1, Nn);
	G = dvector(1, Nn);
	c = dvector(1, Nn);
	echo = dvector(0, Nta - 1);
	trans = dvector(0, Nta - 1);
	tveca = dvector(0, Nta - 1);
	q_lv_hist = dvector(0, Nta - 1);
	p_ao_hist = dvector(0, Nta - 1);

	Vaorta = 0.;
	t_transit = 0.;
	for (i = 1; i <= Nn; i++) {
		x[i] = (i - 1) * Laorta / (Nn - 1);
		A0[i] = A00 * exp(-A01 * x[i]);	// A0 values
		G[i] = Gaorta0 * exp(Gaorta1 * x[i]);	// G values
		c[i] = 1. / sqrt(G[i] * rho); // wave speed 
		//A0[i] = A00 * (1 - A01 * x[i]); // linear taper
		if (i == 1 || i == Nn) {
			Vaorta += A0[i] * Laorta / (Nn - 1) / 2.;
			t_transit += 1. / c[i] / 2.;
		}
		else {
			Vaorta += A0[i] * Laorta / (Nn - 1);	//volume of aorta by trapezium rule (cm3)
			t_transit += 1. / c[i];	//wave transit time of aorta by trapezium rule
		}
		p[i] = 0.;
		if (i < Nn) q[i] = 0.;		
	}
	t_transit = t_transit * Laorta / (Nn - 1);
	VELAORTAWAVE = Laorta / t_transit;
	Z0end = rho * c[Nn] / A0[Nn]; // impedance at far end
	Z0start = rho * c[1] / A0[1]; // impedance at near end

	t = 0.;
	for (i = 0; i < Nta; i++) {
		tveca[i] = (i + 1) * dtMax;
		echo[i] = 0.;
		trans[i] = 0.;
	}
	//solve system using RK4 scheme. TWS, May 2022.
	for (i = 1; i <= Nta * Ncycles; i++) {
		for (j = 1; j <= Ntdiv; j++) {
			derivs(p, q, pt1, qt1, A0, G, t, Nn,
				dx, tpulse, impedancefac);	//m1=f(x0,y0);
			for (k = 1; k <= Nn; k++) {
				pold[k] = p[k];
				p[k] += pt1[k] * dt / 2.;
				if (k < Nn) {
					qold[k] = q[k];
					q[k] += qt1[k] * dt / 2.;
				}
			}
			t += dt / 2.;
			derivs(p, q, pt2, qt2, A0, G, t, Nn,
				dx, tpulse, impedancefac);	//m2 = f((x0 + h / 2.0), (y0 + m1 * h / 2.0)); 
			for (k = 1; k <= Nn; k++) {
				p[k] = pold[k] + pt2[k] * dt / 2.;
				if (k < Nn) q[k] = qold[k] + qt2[k] * dt / 2.;
			}
			derivs(p, q, pt3, qt3, A0, G, t, Nn,
				dx, tpulse, impedancefac);	//m3 = f((x0 + h / 2.0), (y0 + m2 * h / 2.0));
			for (k = 1; k <= Nn; k++) {
				p[k] = pold[k] + pt3[k] * dt;
				if (k < Nn) q[k] = qold[k] + qt3[k] * dt;
			}
			t += dt / 2.;
			derivs(p, q, pt4, qt4, A0, G, t, Nn,
				dx, tpulse, impedancefac);	//m4 = f((x0 + h), (y0 + m3 * h));				
			for (k = 1; k <= Nn; k++) {				//y1 = y0 + (m1 + 2 * m2 + 2 * m3 + m4) * h / 6;
				p[k] = pold[k] + (pt1[k] + 2. * pt2[k] + 2. * pt3[k] + pt4[k]) * dt / 6.;
				if (k < Nn) q[k] = qold[k] + (qt1[k] + 2. * qt2[k] + 2. * qt3[k] + qt4[k]) * dt / 6.;
			}
		}
		ii = (i + Nta - Npulse) % Nta;	//map into interval length Nta, starting at peak of pulse
		echo[ii] += p[1];				//accumulate values over Ncycle cycles
		trans[ii] += q[Nn - 1];			//accumulate values over Ncycle cycles
		if (t < tpulse) echo[ii] -= Z0start * (1 - cos(2 * PI * t / tpulse)) / tpulse; //subtract initial pressure signal
	}
	echo_av = 0.;
	for (i = 0; i < Nta; i++) echo_av += echo[i];
	echo_av = echo_av / Nta;
	for (i = 0; i < Nta; i++) echo[i] = (echo[i] - echo_av) / 1.e4;		//convert to kPa

	free_dvector(x, 1, Nn);
	free_dvector(p, 1, Nn);
	free_dvector(pold, 1, Nn);
	free_dvector(pt1, 1, Nn);
	free_dvector(pt2, 1, Nn);
	free_dvector(pt3, 1, Nn);
	free_dvector(pt4, 1, Nn);		
	free_dvector(q, 1, Nn - 1);
	free_dvector(qold, 1, Nn - 1);
	free_dvector(qt1, 1, Nn - 1);
	free_dvector(qt2, 1, Nn - 1);
	free_dvector(qt3, 1, Nn - 1);
	free_dvector(qt4, 1, Nn - 1);		
	free_dvector(G, 1, Nn);
	free_dvector(c, 1, Nn);
	free_dvector(A0, 1, Nn);
	//don't free tvec, q_lv_hist, p_ao_hist, trans and echo - needed for model
	return Nta;
}
