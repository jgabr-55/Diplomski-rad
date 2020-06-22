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


  //pythia.readString("Beams:eCM = 13000"); //energija 
  //pythia.readString("HardQCD:all = on");  // podprocesi zeljeni
  //pythia.readString("HiggsSM:gg2H = on");   // podprocesi zeljeni
 //pythia.readString("HiggsSM:ffbar2HZ = on"); // signal
  pythia.readString("Beams:idB = -2212");
  pythia.readString("Beams:eCM = 14000.");

  pythia.readString("PhaseSpace:mHatMin = 40.");
  pythia.readString("PhaseSpace:pTHatMin = 100."); ////////////////////to variramo

 /* pythia.readString("WeakBosonAndParton:qqbar2Wg = on");  ////// mozda tu
  pythia.particleData.mayDecay(-24, false);
  pythia.particleData.mayDecay(24, false); */

  pythia.readString("HardQCD:qqbar2gg = on"); ///2 jeta
  bool signal = true;

 // pythia.readString("HardQCD:gg2bbbar = on"); // pozadina

  pythia.readString("PartonLevel:ISR = on"); 
  pythia.readString("PartonLevel:FSR = on");  
  pythia.readString("PartonLevel:MPI = on");
  
  pythia.readString("Next:numberShowEvent = 4");
  
  //pythia.readString("HiggsSM:all = on");
  /*pythia.readString("25:m0 = 125.0");
   
  pythia.readString("25:onMode = off");  
  pythia.readString("25:onIfMatch = 5 -5");

  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfMatch = 11 -11");*/

  pythia.readString("HadronLevel:Hadronize = on");
  //pythia.readString("HadronLevel:all = on");
  
 
  //pythia.readString("PhaseSpace:pTHatMax=520");
  
  


  pythia.init();
  

  std::ofstream dat,jet_anti_kt,dipole_kt_dat, test, dipole_hist,antikt_hist, jet_kt, kt_hist, validacija;
  dat.open("true_data_higgs_bbar.txt");

  jet_anti_kt.open("anti_kt_alg_rekonstr.txt");
  dipole_kt_dat.open("dipol_kt_alg_rekonstr.txt");
  jet_kt.open("kt_alg_rekonstr.txt");

  dipole_hist.open("dipole_kt_broj_rek_jetova.txt");
  antikt_hist.open("anti_kt_broj_rek_jetova.txt");
  kt_hist.open("kt_broj_rek_jetova.txt");

  validacija.open("validacija.txt");

  

  
  
  // Fastjet analysis - select algorithm and parameters
  double Rparam = 0.7;
  fastjet::Strategy               strategy = fastjet::Best;
  fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;
  fastjet::JetDefinition         *jetDef = NULL, *jetDef_kt = NULL;
  jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, Rparam,
                                      recombScheme, strategy);

  jetDef_kt = new fastjet::JetDefinition(fastjet::kt_algorithm, Rparam, recombScheme, strategy);
 
  //priprema za dipole_kt

  InterKT::Clustering<Vec4> clus;
  //clus.maxMIpt2 = 1.0;
  //clus.cutkt = 0.0; ne utjece nista
  //InterKT::Calorimeter<Vec4> calo(60, 100, -5.0, 5.0);

  //clus.cutkt = 10.0;
  clus.nmin = 2; //  2
  
  int d1_no = -1,d2_no = -1, h_no=0, testna = -1;
   
  dat<<"px(Higgs)  py(Higgs)  pz(Higgs)  energ(Higgs)  masa(Higgs)  px(cest_1)  py(cest_1)  pz(cest_1)  energ(cest_1)  masa(cest_1)  px(cest_2)  py(cest_2)  pz(cest_2)  energ(cest_2)  masa(cest_2)"<<endl;
  jet_anti_kt<<"px(1, true)   py(1, true)   pz(1, true)   E(1, true)  px(2, true)   py(2, true)   pz(2, true)   E(2, true)   px(1)   py(1)   pz(1)   E(1)  px(2)   py(2)   pz(2)   E(2)"<<endl;
  dipole_kt_dat<<"px(1, true)   py(1, true)   pz(1, true)   E(1, true)  px(2, true)   py(2, true)   pz(2, true)   E(2, true)   px(1)   py(1)   pz(1)   E(1)  px(2)   py(2)   pz(2)   E(2)"<<endl;
  jet_kt<<"px(1, true)   py(1, true)   pz(1, true)   E(1, true)  px(2, true)   py(2, true)   pz(2, true)   E(2, true)   px(1)   py(1)   pz(1)   E(1)  px(2)   py(2)   pz(2)   E(2)"<<endl;
 
  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 1000; ++iEvent) {
    if (!pythia.next()) continue;

    d1_no=-1;
    d2_no=-1;
    testna = -1;

    std::vector <fastjet::PseudoJet> finalParticles;
    std::vector<fastjet::PseudoJet> finalParticles_kt;
    std::vector <Vec4> tracks;
    std::vector <Vec4> others;

    test.open("test.txt");

    
    //multimap<double,Vec4> sorted;
     
if(signal)
{

    for (int i = 0; i < pythia.event.size(); ++i){ //petlja po svakoj čestici (tu možemo gledati njihova svojstva


      if(pythia.event[i].id() == 21 && pythia.event[i].mother1()==3 && pythia.event[i].mother2()==4 && testna == -1) // 25
	{
		//pronalazimo zadnjeg higgsa koji se raspao i spremamo njegove kceri
		/*d1_no = pythia.event[i].daughter1();
		d2_no = pythia.event[i].daughter2();
		h_no = i;*/
		d1_no = i;
		testna = testna + 1;
		//d2_no = i;
		h_no++;
		
	}
      else if(pythia.event[i].id() == 21 && pythia.event[i].mother1()==3 && pythia.event[i].mother2()==4 && testna == 0)
		{d2_no = i;
		 h_no++;}		
		
	//sve finalne cestice osim el/poz ubacujemo u anti_kt i kt alg za rekonstr jetova i dipole kt alg

      //if(pythia.event[i].status()>0 && pythia.event[i].id()!=24 && pythia.event[i].id()!=-24) // 11 , -11
      if(pythia.event[i].status()>0)
	{
		finalParticles.push_back(fastjet::PseudoJet(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e()));

		finalParticles_kt.push_back(fastjet::PseudoJet(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e()));

	  /// donja linija za dipoleKT
		tracks.push_back(pythia.event[i].p()); //// vektor za funkciju cluster!!!!!!! ---> radi
	}

	
       //ubacujemo leptone za dipole
      /*if((pythia.event[i].id()==24 || pythia.event[i].id()==-24) && pythia.event[i].status()>0) // 11 , -11
	{
		others.push_back(pythia.event[i].p());
		//h_no++;
	}*/
      
      //else{}
	
  }
}

////// pronalazimo b i bbar pozadine i ubacujemo sve finalne u algoritme
else
{
 for (int i = 0; i < pythia.event.size(); ++i){ 
      if(pythia.event[i].id() == 5 && pythia.event[i].daughter1()!=pythia.event[i].daughter2())
	{		
		d1_no = i;
 		//cout<<i<<"\t"<<pythia.event[i].id()<<endl;		  
	}
      
      else if(pythia.event[i].id() == -5 && pythia.event[i].daughter1()!=pythia.event[i].daughter2())
	{
		d2_no = i;
		//cout<<i<<endl;
	}

	//sve finalne cestice osim el/poz ubacujemo u anti_kt i kt alg za rekonstr jetova i dipole kt alg

      if(pythia.event[i].status()>0)
	{
		finalParticles.push_back(fastjet::PseudoJet(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e()));

		finalParticles_kt.push_back(fastjet::PseudoJet(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e()));

	  /// donja linija za dipoleKT
		tracks.push_back(pythia.event[i].p()); //// vektor za funkciju cluster!!!!!!! ---> radi
	}
	
  }

}
	
  //jet rekonstrukcija anti_kt
	fastjet::ClusterSequence clustSeq(finalParticles, *jetDef);
        vector <fastjet::PseudoJet> sortedJets;      
	sortedJets = clustSeq.inclusive_jets(10); //broj oznacava min pt
	//sortedJets = clustSeq.exclusive_jets(2); //moze dat error

  //provjeravamo rekonstruirane jetove za svaki dog. --> (anti_kt)
	antikt_hist<<sortedJets.size()<<endl;

  //jet rekonstrukcija kt
	fastjet::ClusterSequence clustSeq_kt(finalParticles_kt, *jetDef_kt);
	vector <fastjet::PseudoJet> jets_kt;
	jets_kt = clustSeq_kt.exclusive_jets(2); // 2

	kt_hist<<jets_kt.size()<<endl;

  //jet rekonstrukcija dipoleKT
	if(signal)
		clus.cluster(tracks,test,others); //dodat others
	else
		clus.cluster(tracks,test);

        std::vector<Vec4> dipolekt_jets = clus.getJets();

  //provjeravamo rekonstruirane jetove za svaki dogadjaj (dipole_kt)
	dipole_hist<< dipolekt_jets.size()<<endl;

  //gledamo samo slucajeve sa 2 jeta (nakon dipole kt alg)
	if(dipolekt_jets.size()==2) // 2
	{
		

		//cout<<dipolekt_jets[0]<<endl;
		dipole_kt_dat<<pythia.event[d1_no].px()<<" "<<pythia.event[d1_no].py()<<" "<<pythia.event[d1_no].pz()<<" "<<pythia.event[d1_no].e()<<" "<<
			       pythia.event[d2_no].px()<<" "<<pythia.event[d2_no].py()<<" "<<pythia.event[d2_no].pz()<<" "<<pythia.event[d2_no].e()<<" "<<
			       dipolekt_jets[0].px()<<" "<<dipolekt_jets[0].py()<<" "<<dipolekt_jets[0].pz()<<" "<<dipolekt_jets[0].e()<<" "<<
			       dipolekt_jets[1].px()<<" "<<dipolekt_jets[1].py()<<" "<<dipolekt_jets[1].pz()<<" "<<dipolekt_jets[1].e()<<endl;

		/*dipole_kt_dat<<pythia.event[d1_no].px()<<" "<<pythia.event[d1_no].py()<<" "<<pythia.event[d1_no].pz()<<" "<<pythia.event[d1_no].e()<<" "<<
			       dipolekt_jets[0].px()<<" "<<dipolekt_jets[0].py()<<" "<<dipolekt_jets[0].pz()<<" "<<dipolekt_jets[0].e()<<endl;*/
		
	}


	if(jets_kt.size() == 2) // 2
	{
		       jet_kt<<pythia.event[d1_no].px()<<" "<<pythia.event[d1_no].py()<<" "<<pythia.event[d1_no].pz()<<" "<<pythia.event[d1_no].e()<<" "<<
			       pythia.event[d2_no].px()<<" "<<pythia.event[d2_no].py()<<" "<<pythia.event[d2_no].pz()<<" "<<pythia.event[d2_no].e()<<" "<<
			       jets_kt[0].px()<<" "<<jets_kt[0].py()<<" "<<jets_kt[0].pz()<<" "<<jets_kt[0].e()<<" "<<
			       jets_kt[1].px()<<" "<<jets_kt[1].py()<<" "<<jets_kt[1].pz()<<" "<<jets_kt[1].e()<<endl;

			/*jet_kt<<pythia.event[d1_no].px()<<" "<<pythia.event[d1_no].py()<<" "<<pythia.event[d1_no].pz()<<" "<<pythia.event[d1_no].e()<<" "<<
			       jets_kt[0].px()<<" "<<jets_kt[0].py()<<" "<<jets_kt[0].pz()<<" "<<jets_kt[0].e()<<endl;*/
	}





  //gledamo samo slucajeve sa 2 jeta  (nakon anti kt alg)
	if(sortedJets.size()==2)
	{

   		///////
		jet_anti_kt<<pythia.event[d1_no].px()<<"     "<<pythia.event[d1_no].py()<<"     "<<pythia.event[d1_no].pz()<<"     "<<pythia.event[d1_no].e()<<"     "<<
			     pythia.event[d2_no].px()<<"     "<<pythia.event[d2_no].py()<<"     "<<pythia.event[d2_no].pz()<<"     "<<pythia.event[d2_no].e()<<"     "<<
			     sortedJets[0].px()<<"     "<<sortedJets[0].py()<<"     "<<sortedJets[0].pz()<<"     "<<sortedJets[0].E()<<"     "<<
			     sortedJets[1].px()<<"     "<<sortedJets[1].py()<<"     "<<sortedJets[1].pz()<<"     "<<sortedJets[1].E()<<endl;
			
	}
	
 //true values od higgsa i podaci za b i bbar direktno iz pythia outputa

if(signal)
{
  /*if(d1_no != -1 && d2_no != -1)
	  {
		//dat<<iEvent<<"\t"<<pythia.event[d1_no].id()<<"\t"<<pythia.event[d2_no].id()<<endl;
		dat<<pythia.event[h_no].px()<<"     "<<pythia.event[h_no].py()<<"     "<<pythia.event[h_no].pz()<<"     "<<pythia.event[h_no].e()<<"     "<<pythia.event[h_no].m()<<"     "<<pythia.event[d1_no].px()<<"     "<<pythia.event[d1_no].py()<<"     "<<pythia.event[d1_no].pz()<<"     "<<pythia.event[d1_no].e()<<"     "<<pythia.event[d1_no].m()<<"     "<<pythia.event[d2_no].px()<<"     "<<pythia.event[d2_no].py()<<"     "<<pythia.event[d2_no].pz()<<"     "<<pythia.event[d2_no].e()<<"     "<<pythia.event[d2_no].m()<<endl;
		
		
          }*/
}

else
{
   if(d1_no != -1 && d2_no != -1)
	  {
		dat<<pythia.event[d1_no].px()<<"     "<<pythia.event[d1_no].py()<<"     "<<pythia.event[d1_no].pz()<<"     "<<pythia.event[d1_no].e()<<"     "<<pythia.event[d1_no].m()
	  <<"     "<<pythia.event[d2_no].px()<<"     "<<pythia.event[d2_no].py()<<"     "<<pythia.event[d2_no].pz()<<"     "<<pythia.event[d2_no].e()<<"     "<<pythia.event[d2_no].m()<<endl;
		
		
          }
}
	
  // End of event loop. Statistics. Histogram. Done.
  test.close();

  //validacija<<pythia.event[d1_no].px()<<"\t"<<pythia.event[d1_no].py()<<"\t"<<pythia.event[d2_no].px()<<"\t"<<pythia.event[d2_no].py()<<endl;
 
  }

  cout<<h_no<<endl;
  
  //cout<<pythia.event[h_no].m()<<endl;
  validacija.close();
  dat.close();
  jet_anti_kt.close();
  dipole_kt_dat.close();
  jet_kt.close();
  kt_hist.close();
  
  dipole_hist.close();
  antikt_hist.close();
  pythia.stat();
  
  //cout << mult;
  return 0;
}
