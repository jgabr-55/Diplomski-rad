#include "Analyzer.h"

void Analyzer::citanje_iz_datoteke(string filename)
{
  TCanvas *c = new TCanvas("c","c",1000,500);
  c->Divide(2,1);
  higgs_m_rekonstr = new TH1F("masa_rekonstr_higgsa","",50,110,130);		//od 2 b kvarka
  higgs_pt_rekonstr = new TH1F("trans_mom_rekonstr_higgsa","",50,-50,150);

  higgs_m_pyt_outpt= new TH1F("masa_original_higgsa","",50,110,130);		//iz pythia outpt direktno
  higgs_pt_pyt_outpt = new TH1F("trans_mom_original_higgsa","",50,-50,150);

  higgs_pt_true = new TH1F("true_higgs_pt","",50,-50,150);
  higgs_m_true = new TH1F("true_higgs_m","",50,100,150);			//iz 2 b kvarka sa >20 GeV

  ifstream myReadFile;
  myReadFile.open(filename.c_str());
  string line;

  _skipFirstLine = true;

  if (myReadFile.is_open())
  {
    // Read the file line by line
    while(getline(myReadFile, line))
    {
        stringstream   linestream(line);

        // Since we know that first line of the file only describes data skip reading it into variables
        if (_skipFirstLine)
        {
          _skipFirstLine = false;
          continue;
        }

        // Read output and send it to dedicated variables
        linestream >> Higgs_px >> Higgs_py >> Higgs_pz >> Higgs_en >> Higgs_m >> particle1_px >> particle1_py >> particle1_pz >> particle1_en >> particle1_m >> particle2_px >> particle2_py >> particle2_pz >> particle2_en >> particle2_m;

	
	b.SetPxPyPzE(particle1_px,particle1_py,particle1_pz,particle1_en);
	bbar.SetPxPyPzE(particle2_px,particle2_py,particle2_pz,particle2_en);
	
	
	/////////////////
        higgs_pyt_outpt.SetPxPyPzE(Higgs_px,Higgs_py,Higgs_pz,Higgs_en);

	higgs_rekonstr = b + bbar;

        higgs_m_rekonstr->Fill(higgs_rekonstr.M());
	higgs_pt_rekonstr->Fill(higgs_rekonstr.Pt());

		
	higgs_m_pyt_outpt->Fill(higgs_pyt_outpt.M());
        higgs_pt_pyt_outpt->Fill(higgs_pyt_outpt.Pt());
	/////////////////

        if(b.Pt()>20 && bbar.Pt()>20)
          {
		higgs_true = b + bbar; //za true vrijednost uzmemo ono što dobijemo od 2 b kvarka

  		higgs_pt_true->Fill(higgs_true.Pt());
		higgs_m_true->Fill(higgs_true.M());
	  }

    }
  }

  myReadFile.close();
  
  

  //cout<<particle2_en<<endl;
  c->cd(1);
  higgs_m_rekonstr->GetXaxis()->SetTitle("Masa");
  higgs_m_rekonstr->SetTitle("Masa higgsa");
  higgs_m_rekonstr->Draw();
  higgs_m_pyt_outpt->SetLineColor(2);
  higgs_m_pyt_outpt->Draw("same");
  c->cd(2);
  higgs_pt_rekonstr->GetXaxis()->SetTitle("Transverzalni moment");
  higgs_pt_rekonstr->SetTitle("Trans. moment higgsa");
  higgs_pt_rekonstr->Draw();
  higgs_pt_pyt_outpt->SetLineColor(2);
  higgs_pt_pyt_outpt->Draw("same");
  c->SaveAs("Masa_rekonstruiranog_higgsa.png");

  
}


void Analyzer::anti_kt_histo(string filename)
{
  TCanvas *ca = new TCanvas("ca","ca",1000,500);
  ca->Divide(2,1);

   histo_akt_pt = new TH1F("anti_kt_alg_pt","",50,-50,150);
   histo_akt_m = new TH1F("anti_kt_alg_,m","",50,100,150);
	
  ifstream myReadFile;
  myReadFile.open(filename.c_str());
  string line;

  _skipFirstLine = true;

  if (myReadFile.is_open())
  {
    // Read the file line by line
    while(getline(myReadFile, line))
    {
        stringstream   linestream(line);
      	
	if (_skipFirstLine)
        {
          _skipFirstLine = false;
          continue;
        }


        // Read output and send it to dedicated variables
        linestream >> particle1_px >> particle1_py >> particle1_pz >> particle1_en >> particle2_px >> particle2_py >> particle2_pz >> particle2_en;

	b.SetPxPyPzE(particle1_px,particle1_py,particle1_pz,particle1_en);
	bbar.SetPxPyPzE(particle2_px,particle2_py,particle2_pz,particle2_en);

	higgs_rekonstr = b + bbar;
	

	histo_akt_pt->Fill(higgs_rekonstr.Pt());
	histo_akt_m->Fill(higgs_rekonstr.M());

    }
  }

 ca->cd(1);
 histo_akt_pt->SetTitle("anti_kt Pt vs True Pt");
 higgs_pt_true->SetLineColor(2);
 TLegend *legend1 = new TLegend(0.1,0.8,0.3,0.9);

  legend1->AddEntry(histo_akt_pt,"anti_kt","l");
  legend1->AddEntry(higgs_pt_true,"True values","l");

 histo_akt_pt->GetXaxis()->SetTitle("Pt");
 histo_akt_pt->GetYaxis()->SetRangeUser(0.,30.);
 histo_akt_pt->Draw();
 higgs_pt_true->Draw("same");
 legend1->Draw("same");

 //higgs_pt_original->GetYaxis()->SetRange(0,30);
 ca->cd(2);
 histo_akt_m->SetTitle("anti_kt Mass vs True Mass");
 higgs_m_true->SetLineColor(2);
 TLegend *legend2 = new TLegend(0.1,0.8,0.3,0.9);

  legend2->AddEntry(histo_akt_m,"anti_kt","l");
  legend2->AddEntry(higgs_m_true,"True values","l");

 histo_akt_m->GetXaxis()->SetTitle("Mass");
 histo_akt_m->Draw();
 histo_akt_m->GetYaxis()->SetRangeUser(0.,60.);
 higgs_m_true->Draw("same");
 legend2->Draw("same");
 ca->SaveAs("anti_kt.png");
}
