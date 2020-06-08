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
#include <TLegend.h>


using namespace std;

class Analyzer{

public:
	void citanje_iz_datoteke(string filename);
	void anti_kt_histo(string filename);
	void dipole_kt_histo(string filename);
	void kt_histo(string filename);
	void histogram_rek_jetova(string filename1, string filename2, string filename3);
	void Crtanje();
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
	
	double particle1_px_true;
	double particle1_py_true;
	double particle1_pz_true;
	double particle1_en_true;

	double particle2_px_true;
	double particle2_py_true;
	double particle2_pz_true;
	double particle2_en_true;

	double Higgs_px;
	double Higgs_py;
	double Higgs_pz;
	double Higgs_en;
	double Higgs_m;

	double temp;

	TLorentzVector b,bbar,b_true, bbar_true, higgs_true,higgs_rekonstr,higgs_pyt_outpt;
	TH1F *higgs_pt_rekonstr, *higgs_m_rekonstr, *higgs_rap_rekonstr, *higgs_pt_pyt_outpt, *higgs_m_pyt_outpt, *higgs_rap_pyt_outpt, *histo_akt_pt, *histo_akt_rap, *histo_dkt_pt, *histo_dkt_rap, *histo_dipole, *histo_anti_kt,*histo_kt, *histo_kt_pt, *histo_kt_rap, *higgs_m_nakon_kt, *higgs_m_nakon_akt, *higgs_m_nakon_dkt;
	
};
