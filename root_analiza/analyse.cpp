#include "Analyzer.h"

using namespace std;

int main(){

	Analyzer *a = new Analyzer();
	//a->citanje_iz_datoteke("higgs_bb_raspad.txt");
	a->dipole_kt_histo("dipol_kt_alg.txt");
	a->anti_kt_histo("anti_kt_rekonstr_bbbar.txt");
	a->kt_histo("kt_alg_rekonstr.txt");
	
	a->histogram_rek_jetova("broj_rek_jetova_po_dog_dipolekt.txt", "kt_broj_rek_jetova.txt");
	a->Crtanje();

	
}
