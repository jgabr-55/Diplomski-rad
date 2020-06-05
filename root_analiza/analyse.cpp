#include "Analyzer.h"

using namespace std;

int main(){

	Analyzer *a = new Analyzer();
	a->citanje_iz_datoteke("true_data_higgs_bbar.txt");
	a->dipole_kt_histo("dipol_kt_alg_rekonstr.txt");
	a->anti_kt_histo("anti_kt_alg_rekonstr.txt");
	a->kt_histo("kt_alg_rekonstr.txt");
	
	a->histogram_rek_jetova("dipole_kt_broj_rek_jetova.txt", "anti_kt_broj_rek_jetova.txt","kt_broj_rek_jetova.txt");
	a->Crtanje();

	
}
