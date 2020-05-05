#include "Analyzer.h"

using namespace std;

int main(){

	Analyzer *a = new Analyzer();
	a->citanje_iz_datoteke("higgs_bb_raspad.txt");
	a->anti_kt_histo("broj_jetova_po_dog.txt");
	
}
