#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <TROOT.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>
#include <TLorentzVector.h>
#include <math.h>


using namespace std;

class Analyzer{

public:
	void citanje_iz_datoteke(string filename);
	void anti_kt_histo(string filename);
	bool _skipFirstLine;

	double particle1_px;
	double particle1_py;
	double particle1_pz;
	double particle1_en;
	double particle1_m;

	double particle2_px;
	double particle2_py;
	double particle2_pz;
	double particle2_en;
	double particle2_m;

	double Higgs_px;
	double Higgs_py;
	double Higgs_pz;
	double Higgs_en;
	double Higgs_m;

	double temp;

	TLorentzVector b,bbar,higgs_rekonstr,higgs_original;
	TH1F *higgs_pt_rekonstr, *higgs_m_rekonstr, *higgs_pt_original, *higgs_m_original, *histo_akt;
	
};
