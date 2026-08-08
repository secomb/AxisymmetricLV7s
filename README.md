# AxisymmetricLV7s
Fast computational simulation of cardiovascular system, with axisymmetric left ventricle (Version 7). This version of the program includes parameter estimation, based on target values of output variables, and uncertainty quantification. For details, including definitions of parameters and output variables, see: A reference model for adult human cardiovascular mechanics, by Timothy W. Secomb and Michael J. Moulton (in review).

Multiple sets of input files, labeled 000, 001, 002 etc. can be processed. In these files, lines starting // are ignored.

A number of switches are specified in parameters000.txt. If switches are off, a single forward simulation is performed. The switches are as follows:  
Optimize = on: Run parameter optimization using the Levenberg-Marquardt algorithm  
Skip = on: Skip this file, useful when processing multiple files  
MCMC = on: Use Markov Chain Monte Carlo to find parameter distributions. Requires values for
number of steps for the Markov process, number of steps for the "burn-in", number of steps between reports of statistics, standard deviation for normal random proposal steps.  
SimpPulm = on: Use simplified pulmonary model, as in the above referenced work. SimpPulm = off: Use the model of Moulton and Secomb (2023).  
CalcFiber = on: Calculate fiber averages of stretch and active stress as functions of time.  
SensitivityMatrix = on: Calculate sensitivity matrix of output variables with respect to parameters  
Note: Only one of the switches Optimize, MCMC and SensitivityMatrix should be set "on."

Initial parameter values are specified in parameters000.txt. If parameter "param" is to be optimized, then the line  
ipar_param = n  
should be included in the file, where the nth multiplicative factor xx[n] is optimized. This allows for more than one input parameter to be varied by a single factor. The factors are initially set to 1 and their upper and lower bounds are specified in the last section of the file.

Target values of output properties are specified in targets000.txt. Target "tar" is included in the objective function if the line  
sigma_rel_tar = s  
is included in the file, where s is a relative standard deviation associated with the determination of tar, for example 0.1. The weights of each target in the objective function are given by 1/s^2.

Multiple output files are generated, as follows.  
graphs000.ps: A set of graphs of relevant quantities in PostScript  
Histogram_xx000.txt.ps: If MCMC is run, another PostScript file gives probability distributions for the multiplicative factors xx[n].   
outputFile_All000.txt: All variables as functions of time, can be copied into a spreadsheet   
Sensitivity000.txt: Sensitivity matrix  
residuals000.txt: Fitted and target values of output variables, showing deviations  
Newparameters000.txt: Updated parameter values after optimization. This can be copied into parameters000.txt to carry out further optimization.   
MCMC_Summary000.txt: Summary of MCMC statistics, including relative standard deviations of estimates   
Histogram_xx000.txt: Data used to generate MCMC histograms  
chisqVals000.txt: Variation of objective function during optimization 

This code runs under Visual Studio 22 C++. Please contact Timothy W. Secomb secomb@arizona.edu if futher informaton is needed.
