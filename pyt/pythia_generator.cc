// main01.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk.
// It studies the charged multiplicity distribution at the LHC.


#include "Pythia8/Pythia.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include <math.h>
#include <fstream>




using namespace Pythia8;

int main() {
  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  pythia.readString("Beams:eCM = 13000"); //energija 
  //pythia.readString("HardQCD:all = on");  // podprocesi zeljeni
  pythia.readString("HiggsSM:gg2H = on");   // podprocesi zeljeni

  pythia.readString("PartonLevel:ISR = off"); 
  pythia.readString("PartonLevel:FSR = off");  
  pythia.readString("PartonLevel:MPI = off");
  



  //pythia.readString("HiggsSM:all = on");
//pythia.readString("HadronLevel:all= on");


  //pythia.readString("25:m0 = 200.0");
  pythia.readString("25:onMode = off");
  
  pythia.readString("25:onIfMatch = 5 -5");

  pythia.readString("HadronLevel:Hadronize = off");

  //pythia.readString("23:onMode = off");

  /*pythia.readString("5:onMode = off");
  pythia.readString("-5:onMode = off");*/

  /*pythia.readString("23:onMode = off");

  pythia.readString("23:onIfAll = 13 -13");*/



  
  //pythia.readString("PartonLevel:all = off");

  //pythia.readString("25:onMode = off");
  //pythia.readString("25:onIfMatch = 23 23");
  //pythia.readString("25:onIfMatch = 5 -5");
  ///pythia.readString("25:onIfMatch = -5 5");

  //pythia.readString("HiggsSM:all = on");
  //pythia.readString("Next:numberShowEvent = 4");   // display n događaja
  //pythia.readString("PhaseSpace:pTHatMin = 20.");
  pythia.init();
  Hist mult("charged multiplicity", 100, -0.5, 799.5);

  ofstream dat;
  dat.open("higgs_bb_raspad.txt");


  // Fastjet analysis - select algorithm and parameters
  /*double Rparam = 0.4;
  fastjet::Strategy               strategy = fastjet::Best;
  fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;
  fastjet::JetDefinition         *jetDef = NULL;
  jetDef = new fastjet::JetDefinition(fastjet::kt_algorithm, Rparam,
                                      recombScheme, strategy);
  // Fastjet input
  std::vector <fastjet::PseudoJet> fjInputs;*/


  string ime, daughter1, daughter2;
  int d1_id,d2_id, brojac = 0;
  string ime_nep;
  double m_higgs;
	
  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 1; ++iEvent) {
    if (!pythia.next()) continue;
    // Find number of all final charged particles and fill histogram.
    
    for (int i = 0; i < pythia.event.size(); ++i){ //petlja po svakoj čestici (tu možemo gledati njihova svojstva
      if(pythia.event[i].id() == 25)
	{

		ime = pythia.event[i].name();
		m_higgs = pythia.event[i].m();
		d1_id = pythia.event[i].daughter1();
		d2_id = pythia.event[i].daughter2();
		
		/*if(d1_id==5)
			if(d2_id==-5)
				brojac++;

		if(d1_id==-5)
			if(d2_id==5)
				brojac++;*/

		if(d1_id == 23 || d2_id == 23)
			brojac++;



		dat<<ime<<"\t"<<d1_id<<"\t"<<d2_id<<"\t"<<m_higgs<<endl;
		  
	}

	
	
  }
	
  // End of event loop. Statistics. Histogram. Done.
  }

  
  dat<<brojac<<endl<<endl;
  dat.close();
 
	
  pythia.stat();
  cout << mult;
  return 0;
}
