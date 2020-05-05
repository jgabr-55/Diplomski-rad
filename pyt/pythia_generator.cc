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
  //pythia.readString("HiggsSM:gg2H = on");   // podprocesi zeljeni
  pythia.readString("HiggsSM:ffbar2HZ = on");

  pythia.readString("PartonLevel:ISR = off"); 
  pythia.readString("PartonLevel:FSR = off");  
  pythia.readString("PartonLevel:MPI = off");
  
  //pythia.readString("HiggsSM:all = on");
  pythia.readString("25:m0 = 125.0");
   
  pythia.readString("25:onMode = off");  
  pythia.readString("25:onIfMatch = 5 -5");

  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfMatch = 11 -11");

  pythia.readString("HadronLevel:Hadronize = on");
 
  pythia.readString("PhaseSpace:pTHatMin=100");
  pythia.readString("PhaseSpace:pTHatMax=120");
  
  pythia.init();
  //Hist mult("charged multiplicity", 100, -0.5, 799.5);

  ofstream dat,jet_anti_kt;
  dat.open("higgs_bb_raspad.txt");
  jet_anti_kt.open("broj_jetova_po_dog.txt");


  // Fastjet analysis - select algorithm and parameters
  double Rparam = 0.3;
  fastjet::Strategy               strategy = fastjet::Best;
  fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;
  fastjet::JetDefinition         *jetDef = NULL;
  jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, Rparam,
                                      recombScheme, strategy);
  // Fastjet input
  


 //u histogram idu masa i pT higgsa///////////////////////////


  
  int d1_no = -1,d2_no = -1, h_no=-1, anti_kt_brojac = 0, pom_akt_brojac=0;
   
  
 
  dat<<"px(Higgs)  py(Higgs)  pz(Higgs)  energ(Higgs)  masa(Higgs)  px(cest_1)  py(cest_1)  pz(cest_1)  energ(cest_1)  masa(cest_1)  px(cest_2)  py(cest_2)  pz(cest_2)  energ(cest_2)  masa(cest_2)"<<endl;
	
  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 100; ++iEvent) {
    if (!pythia.next()) continue;
    // Find number of all final charged particles and fill histogram.

    d1_no=-1;
    d2_no=-1;

    std::vector <fastjet::PseudoJet> finalParticles;

    for (int i = 0; i < pythia.event.size(); ++i){ //petlja po svakoj čestici (tu možemo gledati njihova svojstva
      if(pythia.event[i].id() == 25 && pythia.event[i].daughter1()!=pythia.event[i].daughter2())
	{
		d1_no = pythia.event[i].daughter1();
		d2_no = pythia.event[i].daughter2();
		h_no = i;
   		//dat<<ime<<"\t"<<d1_id<<"\t"<<d2_id<<"\t"<<m_higgs<<endl;
		  
	}

      if(pythia.event[i].status()>0)
	{
		finalParticles.push_back(fastjet::PseudoJet(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e()));
	}
	
  }
	
  //jet rekonstrukcija
	fastjet::ClusterSequence clustSeq(finalParticles, *jetDef);
        vector <fastjet::PseudoJet> sortedJets;      
	sortedJets = clustSeq.inclusive_jets(20);
	jet_anti_kt<<sortedJets.size()<<endl;
	

  if(d1_no != -1 && d2_no != -1)
	  {
		//dat<<iEvent<<"\t"<<pythia.event[d1_no].id()<<"\t"<<pythia.event[d2_no].id()<<endl;
		dat<<pythia.event[h_no].px()<<"     "<<pythia.event[h_no].py()<<"     "<<pythia.event[h_no].pz()<<"     "<<pythia.event[h_no].e()<<"     "<<pythia.event[h_no].m()<<"     "<<pythia.event[d1_no].px()<<"     "<<pythia.event[d1_no].py()<<"     "<<pythia.event[d1_no].pz()<<"     "<<pythia.event[d1_no].e()<<"     "<<pythia.event[d1_no].m()<<"     "<<pythia.event[d2_no].px()<<"     "<<pythia.event[d2_no].py()<<"     "<<pythia.event[d2_no].pz()<<"     "<<pythia.event[d2_no].e()<<"     "<<pythia.event[d2_no].m()<<endl;
		
          }
	
  // End of event loop. Statistics. Histogram. Done.
  }
  
  //stvarnje objekta za clustring i spremanje clustera u vektor

  

  
/*cout<<"alg: "<<jetDef->description()<<endl;
  for(unsigned int i = 0; i<sortedJets.size();i++)
  {
	cout<<sortedJets[i].pt()<<endl;

  }*/

  
  dat.close();
  jet_anti_kt.close();
	
  pythia.stat();
  //cout << mult;
  return 0;
}
