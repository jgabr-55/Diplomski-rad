#include "Analyzer.h"

using namespace std;

int main(){
	
	bool signal = false;
	Analyzer *a = new Analyzer();
	a->citanje_iz_datoteke("true_data_higgs_bbar.txt", signal);
	a->dipole_kt_histo("dipol_kt_alg_rekonstr.txt",signal);
	a->anti_kt_histo("anti_kt_alg_rekonstr.txt",signal);
	a->kt_histo("kt_alg_rekonstr.txt",signal);
	
	a->histogram_rek_jetova("dipole_kt_broj_rek_jetova.txt", "anti_kt_broj_rek_jetova.txt","kt_broj_rek_jetova.txt");
	if(signal == false)
		a->Crtanje_pozadina();
	else
		a->Crtanje_signal();

	
}
