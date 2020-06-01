// main01.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk.
// It studies the charged multiplicity distribution at the LHC.

//kompajler vidi promjenu u .h i .cc samo ako napravimo i promjenu u ovom kodu

#include "Pythia8/Pythia.h"
#include "InterKTPythia8.h"
#include "InterKT.h"

#include <fenv.h>
#include <map>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>




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
  pythia.readString("25:m0 = 100");
   
  pythia.readString("25:onMode = off");  
  pythia.readString("25:onIfMatch = 5 -5");

  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfMatch = 11 -11");

  pythia.readString("HadronLevel:Hadronize = on");
 
  pythia.readString("PhaseSpace:pTHatMin=100");
  pythia.readString("PhaseSpace:pTHatMax=120");
  
  pythia.init();
  //Hist mult("charged multiplicity", 100, -0.5, 799.5);

  std::ofstream dat,jet_anti_kt,dipole_kt_dat, test, dipole_hist,antikt_hist;
  dat.open("higgs_bb_raspad.txt");
  jet_anti_kt.open("anti_kt_rekonstr_bbbar.txt");
  dipole_kt_dat.open("dipol_kt_alg.txt");

  dipole_hist.open("broj_rek_jetova_po_dog_dipolekt.txt");
  antikt_hist.open("antikt_broj_rek_jetova.txt");

  test.open("test.txt");



  // Fastjet analysis - select algorithm and parameters
  double Rparam = 0.3;
  fastjet::Strategy               strategy = fastjet::Best;
  fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;
  fastjet::JetDefinition         *jetDef = NULL;
  jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, Rparam,
                                      recombScheme, strategy);
 
  //priprema za dipole_kt

  InterKT::Clustering<Vec4> clus;
  
  clus.cutkt = 10.0;
  //clus.nmin = 2;
  
  int d1_no = -1,d2_no = -1, h_no=-1;
   
  dat<<"px(Higgs)  py(Higgs)  pz(Higgs)  energ(Higgs)  masa(Higgs)  px(cest_1)  py(cest_1)  pz(cest_1)  energ(cest_1)  masa(cest_1)  px(cest_2)  py(cest_2)  pz(cest_2)  energ(cest_2)  masa(cest_2)"<<endl;
  jet_anti_kt<<"px(1)   py(1)   pz(1)   E(1)  px(2)   py(2)   pz(2)   E(2)"<<endl;
  dipole_kt_dat<<"px(1)   py(1)   pz(1)   E(1)  px(2)   py(2)   pz(2)   E(2)"<<endl;

  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 1000; ++iEvent) {
    if (!pythia.next()) continue;

    d1_no=-1;
    d2_no=-1;

    std::vector <fastjet::PseudoJet> finalParticles;
    std::vector<Vec4> tracks;
    multimap<double,Vec4> sorted;

    for (int i = 0; i < pythia.event.size(); ++i){ //petlja po svakoj čestici (tu možemo gledati njihova svojstva
      if(pythia.event[i].id() == 25 && pythia.event[i].daughter1()!=pythia.event[i].daughter2())
	{
		//pronalazimo zadnjeg higgsa koji se raspao i spremamo njegove kceri
		d1_no = pythia.event[i].daughter1();
		d2_no = pythia.event[i].daughter2();
		h_no = i;
		  
	}

	//sve finalne cestice osim el/poz ubacujemo u alg za rekonstr jetova i dipole kt alg

      if(pythia.event[i].status()>0 && pythia.event[i].id()!=11 && pythia.event[i].id()!=-11)
	{
		finalParticles.push_back(fastjet::PseudoJet(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e()));

	  /// donja linija za dipoleKT
		tracks.push_back(pythia.event[i].p()); //// vektor za funkciju cluster!!!!!!! ---> radi
	}
	
  }
	
  //jet rekonstrukcija anti_kt
	fastjet::ClusterSequence clustSeq(finalParticles, *jetDef);
        vector <fastjet::PseudoJet> sortedJets;      
	sortedJets = clustSeq.inclusive_jets(10); //broj oznacava min pt

	antikt_hist<<sortedJets.size()<<endl;

  //jet rekonstrukcija dipoleKT

	
	clus.cluster(tracks,test);
        std::vector<Vec4> dipolekt_jets = clus.getJets();

	dipole_hist<< dipolekt_jets.size()<<endl;

  //gledamo samo slucajeve sa 2 jeta (nakon dipole kt alg)
	/*if(dipolekt_jets.size()==2)
	{
		//cout<<dipolekt_jets[0]<<endl;
		dipole_kt_dat<<dipolekt_jets[0].px()<<"     "<<dipolekt_jets[0].py()<<"     "<<dipolekt_jets[0].pz()<<"     "<<dipolekt_jets[0].e()<<"     "<<
			       dipolekt_jets[1].px()<<"     "<<dipolekt_jets[1].py()<<"     "<<dipolekt_jets[1].pz()<<"     "<<dipolekt_jets[1].e()<<endl;
		
	}*/





  //gledamo samo slucajeve sa 2 jeta  (nakon anti kt alg)
	/*if(sortedJets.size()==2)
	{
		jet_anti_kt<<sortedJets[0].px()<<"     "<<sortedJets[0].py()<<"     "<<sortedJets[0].pz()<<"     "<<sortedJets[0].E()<<"     "<<
			     sortedJets[1].px()<<"     "<<sortedJets[1].py()<<"     "<<sortedJets[1].pz()<<"     "<<sortedJets[1].E()<<endl;
			
	}*/
	
 //true values od higgsa i podaci za b i bbar direktno iz pythia outputa
 //za navedene uvjete upisujemo u sve 3 datoteke tako da u svima bude isti broj događaja koji odgovaraju jedan drugom

  if(d1_no != -1 && d2_no != -1 && sortedJets.size()==2 && dipolekt_jets.size()==2)
	  {
		//dat<<iEvent<<"\t"<<pythia.event[d1_no].id()<<"\t"<<pythia.event[d2_no].id()<<endl;
		dat<<pythia.event[h_no].px()<<"     "<<pythia.event[h_no].py()<<"     "<<pythia.event[h_no].pz()<<"     "<<pythia.event[h_no].e()<<"     "<<pythia.event[h_no].m()<<"     "<<pythia.event[d1_no].px()<<"     "<<pythia.event[d1_no].py()<<"     "<<pythia.event[d1_no].pz()<<"     "<<pythia.event[d1_no].e()<<"     "<<pythia.event[d1_no].m()<<"     "<<pythia.event[d2_no].px()<<"     "<<pythia.event[d2_no].py()<<"     "<<pythia.event[d2_no].pz()<<"     "<<pythia.event[d2_no].e()<<"     "<<pythia.event[d2_no].m()<<endl;


	//upis jetova iz dipol_kt alg
		dipole_kt_dat<<dipolekt_jets[0].px()<<"     "<<dipolekt_jets[0].py()<<"     "<<dipolekt_jets[0].pz()<<"     "<<dipolekt_jets[0].e()<<"     "<<
			       dipolekt_jets[1].px()<<"     "<<dipolekt_jets[1].py()<<"     "<<dipolekt_jets[1].pz()<<"     "<<dipolekt_jets[1].e()<<endl;


	//upis jetova iz anti_kt alg
		jet_anti_kt<<sortedJets[0].px()<<"     "<<sortedJets[0].py()<<"     "<<sortedJets[0].pz()<<"     "<<sortedJets[0].E()<<"     "<<
			     sortedJets[1].px()<<"     "<<sortedJets[1].py()<<"     "<<sortedJets[1].pz()<<"     "<<sortedJets[1].E()<<endl;
		
		
          }
	
  // End of event loop. Statistics. Histogram. Done.
  }
  

  dat.close();
  jet_anti_kt.close();
  dipole_kt_dat.close();
  test.close();
  dipole_hist.close();
  antikt_hist.close();
  pythia.stat();
  //cout << mult;
  return 0;
}
