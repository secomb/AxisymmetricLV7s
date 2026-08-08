// readParameters.cpp
// Reads parameters.txt
// TWS, July 2025

#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "nrutil.h"

using namespace std;

// This finds a single value from a line that contains the key "string"
// Doesn't read lines that contain "//"
double findValue(ifstream& fileID, string key) {
	double value;
	bool success = 0;
	string str;
	size_t foundKey;
	stringstream numStr;
	fileID.clear();
	fileID.seekg(0, ios::beg);

	while (std::getline(fileID, str)) {				// look through file for the key
		foundKey = str.find(key);
		char* checkChar = new char[str.length() + 1];
		std::strcpy(checkChar, str.c_str());
		size_t foundEscape = str.find("//");

		if (foundKey != string::npos
			&& (foundKey == 0 || isspace(checkChar[foundKey - 1]))
			&& (isspace(checkChar[foundKey + key.length()]) || (checkChar[foundKey + key.length()] == '='))
			&& foundEscape == string::npos) {

			// erase the key, as it might contain a digit
			str.erase(foundKey, key.length());

			// once the key is found, search for the number
			char* strChar = new char[str.length() + 1];
			std::strcpy(strChar, str.c_str());
			int i = 0;
			char c;
			char cprev = 'e';
			while (strChar[i]) {
				c = strChar[i];
				if (isdigit(c) || c == '.' || (c == 'e' && isdigit(cprev)) || c == '+' || c == '-') numStr << c;
				else if (numStr.str().length() > 0) break;
				cprev = c;
				i++;
			}
			numStr >> value;
			success = 1;
			delete[] strChar;
		}
		delete[] checkChar;
	}
	if (success == 1) printf(""); //printf("Set %s = %6f\n", key.c_str(), value);
	else {
		//cout << "*** Parameter " << key << " not found. Set to 0. ***" << endl;
		value = 0;
	}
	return value;
}

// This finds a vector of values from a line that contains the key "string"
void findArray(ifstream& fileID, string key, int N, double* values) {
	double success = 0;
	string str;
	size_t foundKey;
	stringstream numStr;
	fileID.clear();
	fileID.seekg(0, ios::beg);
	int k = 0;
	while (std::getline(fileID, str)) {
		// look through file for the key
		foundKey = str.find(key);
		char* checkChar = new char[str.length() + 1];
		std::strcpy(checkChar, str.c_str());
		size_t foundEscape = str.find("//");

		if (foundKey != string::npos
			&& (foundKey == 0 || isspace(checkChar[foundKey - 1]))
			&& (isspace(checkChar[foundKey + key.length()]) || (checkChar[foundKey + key.length()] == '='))
			&& foundEscape == string::npos) {

			// erase the key, as it might contain a digit
			str.erase(foundKey, key.length());

			// add something on the end so it doesnt miss the last number
			str.append(" end line ");

			// once the key is found, search for the number
			char* strChar = new char[str.length() + 1];
			std::strcpy(strChar, str.c_str());
			success = 0.5;

			// once the key is found, search for the values
			char c;
			char cprev = 'e';
			int i = 0;

			while (strChar[i]) {
				c = strChar[i];
				if (isdigit(c) || c == '.' || (c == 'e' && isdigit(cprev)) || c == '+' || c == '-') numStr << c;
				else if (numStr.str().length() > 0) {
					numStr >> values[k];
					k = k + 1;
					if (k > N) {
						success = -1;
						break;
					}
					numStr.str(std::string());
					numStr.clear();
				}
				else {
					// if its not a number or a . and no number is stored then just skip
				}
				cprev = c;
				i++;
			}
			delete[] strChar;
		}
		delete[] checkChar;
	}
	// success if k got to n
	if (k == N && success != -1) {
		cout << "Vector " << key << " successfully copied with " << N << " elements." << endl;
		success = 1;
	}
	if (success == 0) cout << "Vector " << key << " not found. Vector arbitrarily set to 0." << endl;
	else if (success == 0.5)
		cout << "Vector " << key << " did not have enough elements. Remaining elements arbitrarily set to 0." << endl;
	else if (success == -1)
		cout << "Vector " << key << " had too many elements. Returning vector with only " << N << " elements." << endl;
}

// This finds a string that follows the key
string findString(ifstream& fileID, string key) {
	string outString;
	bool success = 0;
	string str;
	size_t foundKey;
	stringstream valueStr;
	fileID.clear();
	fileID.seekg(0, ios::beg);

	while (std::getline(fileID, str)) {
		// look through file for the key
		foundKey = str.find(key);
		char* checkChar = new char[str.length() + 1];
		std::strcpy(checkChar, str.c_str());
		size_t foundEscape = str.find("//");
		if (foundKey != string::npos
			&& (foundKey == 0 || isspace(checkChar[foundKey - 1]))
			&& (isspace(checkChar[foundKey + key.length()]) || (checkChar[foundKey + key.length()] == '='))
			&& foundEscape == string::npos) {

			// erase the key, as it might contain a digit
			str.erase(foundKey, key.length());

			// once the key is found, search for the number
			char* strChar = new char[str.length() + 1];
			std::strcpy(strChar, str.c_str());
			int i = 0;
			char c;
			char cprev = 'e';
			while (strChar[i]) {
				c = strChar[i];
				// skip '=' or spaces
				if (c != ' ' && c != '=' && c != '\r') valueStr << c;
				else if (valueStr.str().length() > 0) 	break;
				cprev = c;
				i++;
			}
			outString = valueStr.str();
			success = 1;
			delete[] strChar;
		}
		delete[] checkChar;
	}
	if (success == 1) printf(""); //cout << "Set string " << key << " to \"" << outString << "\"" << endl;
	else {
		cout << "*** String " << key << " not found. String is left empty. ***" << endl;
		outString = "";
	}
	return outString;
}

void readParameters(int run)
{
	extern int Nparam, Nopt, aortaVar;		//test if aorta parameters are being varied
	extern int NstepsMarkov, NburnMarkov, NreportMarkov;
	extern int* ParamVar, * ia;
	extern double sigmaMarkov;
	extern double* Param0, *xx, *xxlb, *xxub, *xxbest;
	extern std::string ParamName[112];	//note: zero-based array, include extra element
	extern string Optimize, Skip, MCMC, SimpPulm, CalcFiber, SensitivityMatrix;

	int i, j;
	double* xx00;

	Nparam = 111;
	Param0 = dvector(1, Nparam);
	ParamVar = ivector(1, Nparam);

	char fname[80];
	sprintf(fname, "parameters%03i.txt", run);
	string InputFile = fname;
	ifstream paramFile(InputFile.c_str(), ios::in);
	SimpPulm = findString(paramFile, "SimpPulm");			// Simplified pulmonary model?

	// Number of cycles
	ParamName[1] = "Ncycles";
	ParamName[2] = "Ncycles_show";

	// Domain discretization parameters. LV: Nmu, Nnu other chambers: Nr
	ParamName[3] = "Nmu";
	ParamName[4] = "Nnu";
	ParamName[5] = "Nr";

	// Time steps
	ParamName[6] = "dtMax";
	ParamName[7] = "dtRecord";

	// density
	ParamName[8] = "rho";

	// Time parameters (s) 
	ParamName[9] = "Tc";

	// *** The following parameters can be used as optimization variables ***
	// Time parameters (s) 
	ParamName[10] = "Ta";
	ParamName[11] = "Tca";
	ParamName[12] = "Tca_shift";

	// Reference configuration parameters (LV)
	ParamName[13] = "a0";
	ParamName[14] = "muin0";
	ParamName[15] = "muout0";
	ParamName[16] = "nu_up";

	// Alternative geometric parameters
	ParamName[17] = "RefVolLV";
	ParamName[18] = "WallVolLV";
	ParamName[19] = "LAtoSA";
	
	// Muscle fiber parameterization angles
	ParamName[20] = "psi_in_b0";
	ParamName[21] = "psi_out_b0";

	// Muscle fiber sarcomere parameters
	ParamName[22] = "Ls0";
	ParamName[23] = "Lsmax";
	ParamName[24] = "Lsw";

	// LV active force generation, activation-dependent viscous resistance
	ParamName[25] = "km";
	ParamName[26] = "kav";
	// LV activation with sarcomere length
	ParamName[27] = "kp";
	// LV, RV steepen activation onset (if kdfac < 1)
	ParamName[28] = "kdfac";

	// LV passive stress parameters: viscous, elastic
	ParamName[29] = "kv";
	ParamName[30] = "c1";
	ParamName[31] = "bff";
	ParamName[32] = "bxx";
	ParamName[33] = "bfx";

	// Right ventricle parameters
	ParamName[34] = "rin0_rv";
	ParamName[35] = "rout0_rv";
	ParamName[36] = "f_rv";
	ParamName[37] = "RefVolRV";
	ParamName[38] = "WallVolRV";
	ParamName[39] = "km_rv";
	ParamName[40] = "kav_rv";
	ParamName[41] = "kv_rv";
	ParamName[42] = "c1_rv";
	ParamName[43] = "br_rv";
	ParamName[44] = "bperp_rv";

	// Left atrium parameters
	ParamName[45] = "rin0_la";
	ParamName[46] = "rout0_la";
	ParamName[47] = "RefVolLA";
	ParamName[48] = "WallVolLA";
	ParamName[49] = "km_la";
	ParamName[50] = "kav_la";
	ParamName[51] = "kv_la";
	ParamName[52] = "c1_la";
	ParamName[53] = "br_la";
	ParamName[54] = "bperp_la";

	// Right atrium parameters
	ParamName[55] = "rin0_ra";
	ParamName[56] = "rout0_ra";
	ParamName[57] = "RefVolRA";
	ParamName[58] = "WallVolRA";
	ParamName[59] = "km_ra";
	ParamName[60] = "kav_ra";
	ParamName[61] = "kv_ra";
	ParamName[62] = "c1_ra";
	ParamName[63] = "br_ra";
	ParamName[64] = "bperp_ra";

	// Valve parameters (Mynard dissertation) area VA in cm2, length VL in cm
	ParamName[65] = "VA_aovo";
	ParamName[66] = "VL_aov";
	ParamName[67] = "VA_mivo";
	ParamName[68] = "VL_miv";
	ParamName[69] = "VA_puvo";
	ParamName[70] = "VL_puv";
	ParamName[71] = "VA_tcvo";
	ParamName[72] = "VL_tcv";
	ParamName[73] = "VA_aovc";
	ParamName[74] = "VA_mivc";
	ParamName[75] = "VA_puvc";
	ParamName[76] = "VA_tcvc";

	//Valve rate constants,1/(kPa s) from Mynard in cm2/dyn/s
	ParamName[77] = "Kaovo";
	ParamName[78] = "Kmivo";
	ParamName[79] = "Kpuvo";
	ParamName[80] = "Ktcvo";
	ParamName[81] = "Kaovc";
	ParamName[82] = "Kmivc";
	ParamName[83] = "Kpuvc";
	ParamName[84] = "Ktcvc";

	// Vein inertia parameters
	ParamName[85] = "VA_pv";
	ParamName[86] = "VL_pv";
	ParamName[87] = "VA_sv";
	ParamName[88] = "VL_sv";

	//aorta parameters: length cm, area at root cm2, exponential taper per cm, compliance (dyn/cm2)^-1, length of flow pulse (s)
	//distal impedance factor, matched = 1, Gaorta0, compliance at entrance Gaorta1, increase of compliance with distance
	ParamName[89] = "Laorta";
	ParamName[90] = "A00";
	ParamName[91] = "A01";
	ParamName[92] = "Gaorta0"; // entrance compliance in (dyn/cm2)^-1
	ParamName[93] = "Gaorta1"; // decrease in compliance with distance
	ParamName[94] = "tpulse"; //length of pulse
	ParamName[95] = "impedancefac";	//distal impedance factor, matched = 1

	// Systemic resistance
	ParamName[96] = "Rsp";
	ParamName[97] = "Rsv";

	// Pulmonary resistance
	ParamName[98] = "Rpa";
	if (SimpPulm == "on") ParamName[99] = "Null";
	else ParamName[99] = "Rpp";
	ParamName[100] = "Rpv";

	// Compliances
	ParamName[101] = "Csa";
	ParamName[102] = "Csv";
	if (SimpPulm == "on") {
		ParamName[103] = "Cpp";
		ParamName[104] = "Null";
	}
	else {
		ParamName[103] = "Cpa";
		ParamName[104] = "Cpv";
	}

	// External tissue pressure
	ParamName[105] = "Pext";

	// Initial pressures (kPa) - determine fluid loading
	ParamName[106] = "Psp0";
	ParamName[107] = "Psv0";
	ParamName[108] = "Ppp0";
	ParamName[109] = "Ppv0";

	// delayed LV and RV inactivation: delay time and shoulder level
	ParamName[110] = "Tdlrv";
	ParamName[111] = "Sdlrv";


	Skip = findString(paramFile, "Skip");		// Skip?
	 
	if (Skip != "on") {
		for (i = 1; i <= Nparam; i++) {
			Param0[i] = findValue(paramFile, ParamName[i]);
			ParamVar[i] = (int)findValue(paramFile, "ipar_" + ParamName[i]);
		}

		aortaVar = 0;
		for (i = 89; i <= 95; i++) if (ParamVar[i]) aortaVar = 1;	//set if any aorta parameter is being varied

		SensitivityMatrix = findString(paramFile, "SensitivityMatrix");			// Optimize?
		Optimize = findString(paramFile, "Optimize");			// Optimize?
		MCMC = findString(paramFile, "MCMC");					// MCMC estimation of parameter distributions?
		if (SensitivityMatrix == "on") {						// Override other settings
			Optimize = "on";
			MCMC = "off";
		}
		CalcFiber = findString(paramFile, "CalcFiber");			// Calculate fiber stretch and stress averages?
		if (MCMC == "on") {
			NstepsMarkov = findValue(paramFile, "NstepsMarkov");	// number of steps for the Markov process
			NburnMarkov = findValue(paramFile, "NburnMarkov");		// number of steps for the Markov process "burn-in"
			NreportMarkov = findValue(paramFile, "NreportMarkov");	// number of steps between reports of results
			sigmaMarkov = findValue(paramFile, "sigmaMarkov");		// standard deviation for normal random proposal steps
		}
		if (Optimize == "on" || MCMC == "on") {			// Determine Nopt
			Nopt = 0;
			for (i = 1; i <= Nparam; i++) if (ParamVar[i] > Nopt) Nopt = ParamVar[i];
			printf("Nopt = %i\n", Nopt);
			if (Nopt > 0) {		//check which optimization parameters are used
				printf("If parameter is used, ia = 1 ... ");
				ia = ivector(1, Nopt);
				for (j = 1; j <= Nopt; j++) {
					ia[j] = 0;
					for (i = 1; i <= Nparam; i++) if (ParamVar[i] == j) ia[j] = 1;
					printf("ia[%i] = %i ", j, ia[j]);
				}
				printf("\n");
				// Read initial values and bounds for multiplicative factors
				xx00 = dvector(0, Nopt - 1);
				xx = dvector(1, Nopt);				//optimization parameters
				xxlb = dvector(1, Nopt);			//lower bounds
				xxub = dvector(1, Nopt);			//upper bounds
				xxbest = dvector(1, Nopt);			//best result to this point
				findArray(paramFile, "xx", Nopt, xx00);
				for (i = 1; i <= Nopt; i++) xx[i] = xx00[i - 1];
				findArray(paramFile, "xxlb", Nopt, xx00);
				for (i = 1; i <= Nopt; i++) xxlb[i] = xx00[i - 1];
				findArray(paramFile, "xxub", Nopt, xx00);
				for (i = 1; i <= Nopt; i++) xxub[i] = xx00[i - 1];
				free_dvector(xx00, 0, Nopt - 1);
			}
		}
		if (Optimize == "on" && MCMC == "on") printf("*** Error: Optimize = on and MCMC = on\n");
	}
	paramFile.close();
}
