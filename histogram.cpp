/*****************************************************
Evaluate histograms of solute levels.  TWS December 07.
Version 2.0, May 1, 2010.
Version 3.0, May 17,2011.
Revised Bohan Li, 2018.
******************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <math.h>
#include "nrutil.h"
#include <stdio.h>
#include <string.h>

void histogram(double** p_tHist, double* xxlb, double* xxub,
	int* ia, int Nopt, int NstepsMarkov, int NburnMarkov, const char* fname1)
{
	int i, j, itp, nctop, ncbot, numbins, imaxstat;
	int ifit, k, n = NstepsMarkov - NburnMarkov;
	float step, xmax, ymax, scalefacx, scalefacy, xshift, yshift;
	float* stepmax, * stat, * cumu, * mstat;
	float min, max, maxstat;
	char histogramimagename[100];
	FILE* ofp, * ofp1;
	double* var;

	var = dvector(1, NstepsMarkov - NburnMarkov);

	strcpy(histogramimagename, fname1);
	strcat(histogramimagename, ".ps");
	ofp = fopen(histogramimagename, "w");
	fprintf(ofp, "%%!PS\n");
	fprintf(ofp, "/m {moveto} def\n");
	fprintf(ofp, "/l {lineto} def\n");
	fprintf(ofp, "0 0 0 setrgbcolor\n");
	fprintf(ofp, "/Times-Roman findfont 7 scalefont setfont\n");
	fprintf(ofp, "0.5 setlinewidth\n");

	ofp1 = fopen(fname1, "w");
	
	ifit = 0;
	for (k = 1; k <= Nopt; k++) if (ia[k]) {
		fprintf(ofp1, "Histogram data for xx[%i]\n", k);
		min = xxlb[k];
		max = xxub[k];
		for (j = 0; j < NstepsMarkov - NburnMarkov; j++) var[j + 1] = min + (max - min) * p_tHist[ifit][j];
		step = 1e6;
		for (i = 0; i <= 36; i++) {
			numbins = (int(floor(max / step)) - int(floor(min / step)) + 1);
			if (numbins >= 20) goto done;	//found number of bins
			if (i % 3 != 1) step = step * 0.5;
			else step = step * 0.4;
		}
		printf("Error: bin size not found for %s\n", fname1);
		return;
	done:;
		nctop = floor(max / step) + 1.;
		ncbot = floor(min / step);
		if (ncbot <= 2 && ncbot > 0) ncbot = 0;
		if (nctop >= -2 && nctop < 0) nctop = 0;
		numbins = nctop - ncbot;
		stepmax = vector(1, numbins);
		stat = vector(1, numbins);
		cumu = vector(1, numbins);
		mstat = vector(1, numbins);
		for (i = 1; i <= numbins; i++) {
			stepmax[i] = step * (i + ncbot);
			mstat[i] = 0;
		}
		for (itp = 1; itp <= n; itp++) {
			for (j = 1; j <= numbins; j++) if (var[itp] < stepmax[j]) {
				mstat[j]++;
				goto binned;
			}
		binned:;
		}
		maxstat = 0;
		for (i = 1; i <= numbins; i++) {
			stat[i] = mstat[i] * 100. / n;
			maxstat = FMAX(maxstat, stat[i]);;
		}
		imaxstat = (maxstat + 5.) / 5.;
		cumu[1] = stat[1];
		for (i = 2; i <= numbins; i++) cumu[i] = cumu[i - 1] + stat[i];

		fprintf(ofp1, "value  %% cumul. %%\n");
		fprintf(ofp1, "%g %7.2f %7.2f\n", step * ncbot, 0., 0.);
		for (i = 1; i <= numbins; i++) fprintf(ofp1, "%g %7.2f %7.2f\n", stepmax[i], stat[i], cumu[i]);

		xmax = max;
		scalefacx = 120. / (numbins * step);
		ymax = imaxstat * 5.;
		scalefacy = 120. / ymax;
		xshift = (ncbot * step) * scalefacx - 140. * float(ifit%4);
		yshift = 150. * float(ifit / 4);

		fprintf(ofp, "/mx {%g mul %f add} def\n", scalefacx, 30. - xshift);
		fprintf(ofp, "/my {%g mul %f add} def\n", scalefacy, 650. - yshift);
		fprintf(ofp, "/mx1 {%f add} def\n", 30. - xshift);
		fprintf(ofp, "/my1 {%f add} def\n", 650. - yshift);
		fprintf(ofp, "newpath\n");
		fprintf(ofp, "%g mx %g my m\n", ncbot * step, ymax);
		fprintf(ofp, "%g mx %g my l\n", ncbot * step, 0.);
		fprintf(ofp, "%g mx %g my l\n", nctop * step, 0.);
		fprintf(ofp, "stroke\n");
		fprintf(ofp, "100 mx1 -20 my1 m (xx[%i]) show\n", k);	//graph label
		for (i = 0; i <= imaxstat; i++)
			fprintf(ofp, "45 mx1 %g my m (%g) show\n", i * 5., i * 5.);	//y-axis labels
		for (i = 0; i <= numbins; i++)	if (i % 2 == 0)
			fprintf(ofp, "%g mx -10 my1 m (%3.1f) show\n", step * (ncbot + i), step * (ncbot + i));	//x-axis labels
		fprintf(ofp, "stroke\n");	//needed to avoid artifact in plot
		fprintf(ofp, "newpath\n");
		for (i = 1; i <= numbins; i++) {
			fprintf(ofp, "%g mx %g my m\n", (i - 1 + ncbot) * step, 0.);
			fprintf(ofp, "%g mx %g my l\n", (i - 1 + ncbot) * step, stat[i]);
			fprintf(ofp, "%g mx %g my l\n", (i + ncbot) * step, stat[i]);
			fprintf(ofp, "%g mx %g my l\n", (i + ncbot) * step, 0.);
		}
		fprintf(ofp, "stroke\n");
		free_vector(stepmax, 1, numbins);
		free_vector(stat, 1, numbins);
		free_vector(cumu, 1, numbins);
		free_vector(mstat, 1, numbins);
		ifit++;
	}
	fprintf(ofp, "showpage\n");
	fclose(ofp);
	fclose(ofp1);
	free_dvector(var, 1, NstepsMarkov - NburnMarkov);
}
