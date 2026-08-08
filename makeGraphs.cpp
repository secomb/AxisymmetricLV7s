// MakeGraphs.cpp
// TWS, November 2019
// Version for "Model2020C" only
#define _CRT_SECURE_NO_WARNINGS	//TWS
#include <stdio.h>
#include <math.h>
#include <string>
#include "nrutil.h"

using namespace std;

void addAxes(FILE *ofp, double yoffsetplot, double xmin, double xmax, double ymin, double ymax) {
	int i;
	double xtick, ytick;
	fprintf(ofp, "0 0 0 setrgbcolor\n");
	double scalefacx = 400. / (xmax - xmin);
	double scalefacy = 250. / (ymax - ymin);
	fprintf(ofp, "/mx {%g mul 50 add} def\n", scalefacx);
	fprintf(ofp, "/my {%g mul 50 add} def\n", scalefacy);
	fprintf(ofp, "0 0 0 setrgbcolor\n");
	fprintf(ofp, "newpath\n");
	fprintf(ofp, "%g mx %g my %g add m\n", 0., 0., yoffsetplot);
	fprintf(ofp, "%g mx %g my %g add l\n", 0., ymax - ymin, yoffsetplot);
	fprintf(ofp, "stroke\n");
	fprintf(ofp, "%g mx %g my %g add m\n", 0., DMAX(-ymin,0.), yoffsetplot);
	fprintf(ofp, "%g mx %g my %g add l\n", xmax - xmin, DMAX(-ymin, 0.), yoffsetplot);
	fprintf(ofp, "stroke\n");
	fprintf(ofp, "/Times-Roman findfont 10 scalefont setfont\n");
	for (i = 0; i <= 10; i++) {
		ytick = i / 10.*(ymax - ymin);
		fprintf(ofp, "30 %g my %g add m (%g) show\n", ytick, yoffsetplot, ymin + ytick);	//y-axis labels
		fprintf(ofp, "47 %g my %g add m\n", ytick, yoffsetplot);			//tick marks
		fprintf(ofp, "50 %g my %g add l\n", ytick, yoffsetplot);
		fprintf(ofp, "stroke\n");
	}
	for (i = 0; i <= 10; i++) {
		xtick = i / 10.*(xmax - xmin);
		fprintf(ofp, "%g mx %g my -10 add %g add m (%g) show\n", xtick, DMAX(-ymin, 0.), yoffsetplot, xmin + xtick);	//x-axis labels
		fprintf(ofp, "%g mx %g my %g add m\n", xtick, DMAX(-ymin, 0.), yoffsetplot);			//tick marks
		fprintf(ofp, "%g mx %g my -3 add %g add l\n", xtick, DMAX(-ymin, 0.), yoffsetplot);
		fprintf(ofp, "stroke\n");
	}
}

void addPlot(FILE *ofp, const char *name, double yoffsetplot, double yoffsetname,
	double *xvals, double *yvals, double xmin, double xmax, double ymin, double ymax, int nvals) {
	int i;
	fprintf(ofp, "250 %g m (%s) show\n", yoffsetname, name);
	for (i = 0; i < nvals; i++) {
		if (i == 0) fprintf(ofp, "%g mx %g my %g add m\n", xvals[i] - xmin, yvals[i] - ymin, yoffsetplot);
		else fprintf(ofp, "%g mx %g my %g add l\n", xvals[i] - xmin, yvals[i] - ymin, yoffsetplot);
	}
	fprintf(ofp, "stroke\n");
}

void makeGraphs(int run) {		//TWS November 2019
	extern int Nt, Nta;
	extern double Rsp, Rpp, f_rv, Csa, Csv, Cpv, Cpa, Cpp;
	extern double *tvec, *a1vec, *a2vec, *a3vec, *a4vec, *a5vec, *a6vec;
	extern double **store, **Y;
	extern double *echo, *trans, *tveca;
	extern string SimpPulm;

	int iRec, pagenum, i;
	double xmax, xmin, ymax, ymin;
	double yoffset = 350., ylabel, ylabel0 = 375., dylabel = 15.;
	double *xvals, *yvals;
	FILE *ofp;
	char fname[80];

	sprintf(fname, "graphs%03i.ps",run);

	int Nt1 = Nt - 1;
	xvals = dvector(0, Nt1);
	yvals = dvector(0, Nt1);

	fopen_s(&ofp, fname, "w");
	fprintf(ofp, "%%!PS-Adobe-2.0\n");
	pagenum = 8;
	fprintf(ofp, "%%%%Pages: %i\n", pagenum);	//number of pages must be specified
	fprintf(ofp, "%%%%EndComments\n");
	pagenum = 0;

	///////////// Page 1: PV loops /////////////////
	pagenum++;
	fprintf(ofp, "%%%%Page: %i %i\n", pagenum, pagenum);
	fprintf(ofp, "/m {moveto} def\n");
	fprintf(ofp, "/l {lineto} def\n");
	ymin = 0;
	ymax = 150.;
	xmin = 0.;
	xmax = 200.;
	ylabel = yoffset + ylabel0;
	addAxes(ofp, yoffset, xmin, xmax, ymin, ymax);
	fprintf(ofp, "1 0 0 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) {
		xvals[iRec] = store[1][iRec];
		yvals[iRec] = store[2][iRec] / 0.1333;
	}
	addPlot(ofp, "P-V loop - LV", yoffset, ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	fprintf(ofp, "0 0 1 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) {
		xvals[iRec] = store[5][iRec];
		yvals[iRec] = store[6][iRec] / 0.1333;
	}
	ylabel -= dylabel;
	addPlot(ofp, "P-V loop - RV", yoffset, ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	////////////////////////
	ymin = 0;
	ymax = 25.;
	xmin = 0.;
	xmax = 100.;
	ylabel = ylabel0;
	addAxes(ofp, 0., xmin, xmax, ymin, ymax);
	fprintf(ofp, "1 0 0 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) {
		xvals[iRec] = store[3][iRec];
		yvals[iRec] = store[4][iRec] / 0.1333;
	}
	addPlot(ofp, "P-V loop - LA", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	fprintf(ofp, "0 0 1 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) {
		xvals[iRec] = store[7][iRec];
		yvals[iRec] = store[8][iRec] / 0.1333;
	}
	ylabel -= dylabel;
	addPlot(ofp, "P-V loop - RA", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	fprintf(ofp, "showpage\n");

	//////////////// Page 2: a1(t), a2(t), a3(t), A(t), valve status	////////////////	
	pagenum++;
	fprintf(ofp, "%%%%Page: %i %i\n", pagenum, pagenum);
	fprintf(ofp, "/m {moveto} def\n");
	fprintf(ofp, "/l {lineto} def\n");
	xmin = floor(tvec[0]);
	xmax = ceil(tvec[Nt1 - 1]);
	ymin = -1.5;
	ymax = 0.5;
	ylabel = yoffset + ylabel0;
	addAxes(ofp, yoffset, xmin, xmax, ymin, ymax);
	for (iRec = 0; iRec < Nt1; iRec++) xvals[iRec] = tvec[iRec];
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = a1vec[iRec];
	fprintf(ofp, "0 0 0 setrgbcolor\n");
	addPlot(ofp, "a1(t)", yoffset, ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = a2vec[iRec] / 10.;
	fprintf(ofp, "1 0 0 setrgbcolor\n");
	addPlot(ofp, "a2(t)/10", yoffset, ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = a3vec[iRec];
	fprintf(ofp, "0 1 0 setrgbcolor\n");
	addPlot(ofp, "a3(t)", yoffset, ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = a4vec[iRec] / 10.;
	fprintf(ofp, "0 0 1 setrgbcolor\n");
	addPlot(ofp, "a4(t)/10", yoffset, ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = f_rv * a5vec[iRec] / 20.;
	fprintf(ofp, "0 1 1 setrgbcolor\n");
	addPlot(ofp, "frv*a5(t)/20", yoffset, ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = a6vec[iRec] / 10.;
	fprintf(ofp, "1 0 1 setrgbcolor\n");
	addPlot(ofp, "a6(t)/10", yoffset, ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	///////////////////////
	ymin = 0.;
	ymax = 9.;
	ylabel = ylabel0;
	addAxes(ofp, 0., xmin, xmax, ymin, ymax);
	fprintf(ofp, "0 0 0 setrgbcolor\n");				//black
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[28][iRec];
	addPlot(ofp, "A(t)", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "1 0 0 setrgbcolor\n");				//red	
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[29][iRec] + 1.2;
	addPlot(ofp, "A_la(t)", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "0 0 1 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[26][iRec] + 2.4;
	addPlot(ofp, "Mitral valve open", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "1 0 1 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[27][iRec] + 3.6;
	addPlot(ofp, "Aortic valve open", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "0.5 0.5 0.5 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = 0.2 + store[33][iRec] + 4.8;
	addPlot(ofp, "Tricuspid valve open", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "0 1 1 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = 0.2 + store[34][iRec] + 6;
	addPlot(ofp, "Pulmonary valve open", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	fprintf(ofp, "showpage\n");

	/////////////////// Page 3: Pressures ////////////////		
	pagenum++;
	fprintf(ofp, "%%%%Page: %i %i\n", pagenum, pagenum);
	fprintf(ofp, "/m {moveto} def\n");
	fprintf(ofp, "/l {lineto} def\n");
	ymin = 0;
	ymax = 150.;
	ylabel = yoffset + ylabel0;
	addAxes(ofp, yoffset, xmin, xmax, ymin, ymax);
	fprintf(ofp, "0 0 0 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[2][iRec] / 0.1333;
	addPlot(ofp, "Plv(t)", yoffset, ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "1 0 0 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[4][iRec] / 0.1333;
	addPlot(ofp, "Pla(t)", yoffset, ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "0 1 0 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = Y[12][iRec] / 0.1333;
	addPlot(ofp, "Psp(t)", yoffset, ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "0 0 1 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = Y[15][iRec] / 0.1333;
	addPlot(ofp, "Ppv(t)", yoffset, ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "1 0 1 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[10][iRec] / 0.1333;
	addPlot(ofp, "Psa(t)", yoffset, ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	////////////////////
	ylabel = ylabel0;
	ymax = 40.;
	addAxes(ofp, 0., xmin, xmax, ymin, ymax);
	fprintf(ofp, "0 0 0 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[6][iRec] / 0.1333;
	addPlot(ofp, "Prv(t)", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "1 0 0 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[12][iRec] / 0.1333;
	addPlot(ofp, "Psv(t)", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "0 1 0 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[8][iRec] / 0.1333;
	addPlot(ofp, "Pra(t)", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "0 0 1 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[14][iRec] / 0.1333;
	addPlot(ofp, "Ppp(t)", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "1 0 1 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[13][iRec] / 0.1333;
	addPlot(ofp, "Ppa(t)", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	fprintf(ofp, "showpage\n");

	/////////////////// Page 4: Volumes ////////////////		
	pagenum++;
	fprintf(ofp, "%%%%Page: %i %i\n", pagenum, pagenum);
	fprintf(ofp, "/m {moveto} def\n");
	fprintf(ofp, "/l {lineto} def\n");
	ymin = 0;
	ymax = 200.;
	ylabel = yoffset + ylabel0;
	addAxes(ofp, yoffset, xmin, xmax, ymin, ymax);
	fprintf(ofp, "0 0 0 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[1][iRec];
	addPlot(ofp, "Vlv(t)", yoffset, ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "1 0 0 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[3][iRec];
	addPlot(ofp, "Vla(t)", yoffset, ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "0 1 0 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[5][iRec];
	addPlot(ofp, "Vrv(t)", yoffset, ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "0 0 1 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[7][iRec];
	addPlot(ofp, "Vra(t)", yoffset, ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	////////////////////
	ymax = 2000.;
	ylabel = ylabel0;
	addAxes(ofp, 0., xmin, xmax, ymin, ymax);
	ylabel -= dylabel;
	fprintf(ofp, "0 0 0 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = Y[12][iRec] * Csa;
	addPlot(ofp, "Vsa(t)", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "1 0 0 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = Y[13][iRec] * Csv;
	addPlot(ofp, "Vsv(t)", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	if (SimpPulm == "on") {
		fprintf(ofp, "0 1 0 setrgbcolor\n");
		for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = Y[14][iRec] * Cpp;
		addPlot(ofp, "Vpp(t)", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
		ylabel -= dylabel;
	}
	else {
		fprintf(ofp, "0 1 0 setrgbcolor\n");
		for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = Y[14][iRec] * Cpa;
		addPlot(ofp, "Vpp(t)", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
		ylabel -= dylabel;
		fprintf(ofp, "0 0 1 setrgbcolor\n");
		for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = Y[15][iRec] * Cpv;
		ylabel -= dylabel;
	}
	fprintf(ofp, "0 1 1 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[47][iRec];
	addPlot(ofp, "Vaorta(t)", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "1 0 1 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[0][iRec];
	addPlot(ofp, "Vtot(t)", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	fprintf(ofp, "showpage\n");

	/////////////////// Page 5: Flows ////////////////		
	pagenum++;
	fprintf(ofp, "%%%%Page: %i %i\n", pagenum, pagenum);
	fprintf(ofp, "/m {moveto} def\n");
	fprintf(ofp, "/l {lineto} def\n");
	ymin = 0;
	ymax = 600.;
	ylabel = yoffset + ylabel0;
	addAxes(ofp, yoffset, xmin, xmax, ymin, ymax);
	fprintf(ofp, "0 0 0 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[17][iRec];
	addPlot(ofp, "q_lv(t)", yoffset, ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "1 0 0 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[16][iRec];
	addPlot(ofp, "q_la(t)", yoffset, ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "0 1 0 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[21][iRec];
	addPlot(ofp, "q_sp(t)", yoffset, ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "0 0 1 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[24][iRec];
	addPlot(ofp, "q_pp(t)", yoffset, ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	//////////////////////
	ylabel = ylabel0;
	ymax = 600.;
	addAxes(ofp, 0., xmin, xmax, ymin, ymax);
	fprintf(ofp, "0 0 0 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[19][iRec];
	addPlot(ofp, "q_rv(t)", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "1 0 0 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[18][iRec];
	addPlot(ofp, "q_ra(t)", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "0 1 0 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[22][iRec];
	addPlot(ofp, "q_sv(t)", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "0 0 1 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[25][iRec];
	addPlot(ofp, "q_pv(t)", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	fprintf(ofp, "showpage\n");

	/////////////////// Page 6: Geometry ////////////////
	pagenum++;
	fprintf(ofp, "%%%%Page: %i %i\n", pagenum, pagenum);
	fprintf(ofp, "/m {moveto} def\n");
	fprintf(ofp, "/l {lineto} def\n");
	ymin = 3.;
	ymax = 9.;
	ylabel = yoffset + ylabel0 - 4. * dylabel;
	fprintf(ofp, "0 0 0 setrgbcolor\n");
	addAxes(ofp, yoffset, xmin, xmax, ymin, ymax);
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[35][iRec];
	addPlot(ofp, "Internal base-apex length", yoffset, ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "1 0 0 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[45][iRec];
	addPlot(ofp, "Internal diameter", yoffset, ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	//////////////////////////		
	ylabel = ylabel0 - 4. * dylabel;
	ymin = -20.;
	ymax = 10.;
	addAxes(ofp, 0., xmin, xmax, ymin, ymax);
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[36][iRec];
	addPlot(ofp, "-d/dt(base-apex length)", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	fprintf(ofp, "showpage\n");

	/////////////////// Page 7: Flow velocities in valves ////////////////
	pagenum++;
	fprintf(ofp, "%%%%Page: %i %i\n", pagenum, pagenum);
	fprintf(ofp, "/m {moveto} def\n");
	fprintf(ofp, "/l {lineto} def\n");
	ymin = 0.;
	ymax = 200.;
	ylabel = yoffset + ylabel0 - 4. * dylabel;
	fprintf(ofp, "0 0 0 setrgbcolor\n");
	addAxes(ofp, yoffset, xmin, xmax, ymin, ymax);
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[41][iRec];
	addPlot(ofp, "Aortic valve velocity (based on open or closed area >= Amin)", yoffset, ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "1 0 0 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[42][iRec];
	addPlot(ofp, "Mitral valve velocity", yoffset, ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	//////////////////////////		
	ylabel = ylabel0 - 4. * dylabel;
	ymin = 0.;
	ymax = 200.;
	addAxes(ofp, 0., xmin, xmax, ymin, ymax);
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[43][iRec];
	addPlot(ofp, "Pulmonary valve velocity", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "1 0 0 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[44][iRec];
	addPlot(ofp, "Tricuspid valve velocity", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	fprintf(ofp, "showpage\n");

	/////////////////// Page 8: Aorta model ////////////////		
	pagenum++;
	fprintf(ofp, "%%%%Page: %i %i\n", pagenum, pagenum);
	fprintf(ofp, "/m {moveto} def\n");
	fprintf(ofp, "/l {lineto} def\n");
	ymin = -0.05;
	ymax = 0.1;
	xmax = tveca[Nta - 1];
	xmin = 0.;
	for (i = 0; i < Nta; i++) trans[i] /= 50.;
	ylabel = yoffset + ylabel0 - 4. * dylabel;
	fprintf(ofp, "0 0 0 setrgbcolor\n");
	addAxes(ofp, yoffset, xmin, xmax, ymin, ymax);
	addPlot(ofp, "echo(t)", yoffset, ylabel, tveca, echo, xmin, xmax, ymin, ymax, Nta);
	ylabel -= dylabel;
	fprintf(ofp, "0 0 1 setrgbcolor\n");
	addPlot(ofp, "trans(t)/50", yoffset, ylabel, tveca, trans, xmin, xmax, ymin, ymax, Nta);
	//////////////////////////		
	ylabel = ylabel0 - 4. * dylabel;
	ymin = -10.;
	ymax = 140.;
	xmax = ceil(tvec[Nt1 - 1]);
	xmin = tvec[0];
	addAxes(ofp, 0., xmin, xmax, ymin, ymax);
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[37][iRec] / 0.1333;
	addPlot(ofp, "P_echo(t)", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "1 0 0 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[38][iRec] / 5.;
	addPlot(ofp, "Q_trans(t)/5", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "0 0 1 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[39][iRec] / 0.1333;
	addPlot(ofp, "P_ao_0(t)", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	ylabel -= dylabel;
	fprintf(ofp, "0 1 0 setrgbcolor\n");
	for (iRec = 0; iRec < Nt1; iRec++) yvals[iRec] = store[40][iRec];
	addPlot(ofp, "q_lv_0(t)", 0., ylabel, xvals, yvals, xmin, xmax, ymin, ymax, Nt1);
	fprintf(ofp, "showpage\n");
	///////////////////////////////////////////////////////////////////////////
	free_dvector(xvals, 0, Nt1);
	free_dvector(yvals, 0, Nt1);
	fclose(ofp);
}
