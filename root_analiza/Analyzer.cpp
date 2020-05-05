#include "Analyzer.h"

void Analyzer::citanje_iz_datoteke(string filename)
{
  TCanvas *c = new TCanvas("c","c",1000,500);
  c->Divide(2,1);
  higgs_m_rekonstr = new TH1F("masa_rekonstr_higgsa","",50,100,150);
  higgs_pt_rekonstr = new TH1F("trans_mom_rekonstr_higgsa","",50,90,140);
  higgs_m_original = new TH1F("masa_original_higgsa","",50,100,150);
  higgs_pt_original = new TH1F("trans_mom_original_higgsa","",50,90,140);

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

        higgs_original.SetPxPyPzE(Higgs_px,Higgs_py,Higgs_pz,Higgs_en);

	higgs_rekonstr = b + bbar;

        higgs_m_rekonstr->Fill(higgs_rekonstr.M());
	higgs_pt_rekonstr->Fill(higgs_rekonstr.Pt());
	
	higgs_m_original->Fill(Higgs_m);
        higgs_pt_original->Fill(sqrt(Higgs_px*Higgs_px+Higgs_py*Higgs_py));

    }
  }

  myReadFile.close();

  //cout<<particle2_en<<endl;
  c->cd(1);
  higgs_m_rekonstr->GetXaxis()->SetTitle("Masa");
  higgs_m_rekonstr->SetTitle("Masa higgsa");
  higgs_m_rekonstr->Draw();
  higgs_m_original->SetLineColor(2);
  higgs_m_original->Draw("same");
  c->cd(2);
  higgs_pt_rekonstr->GetXaxis()->SetTitle("Transverzalni moment");
  higgs_pt_rekonstr->SetTitle("Trans. moment higgsa");
  higgs_pt_rekonstr->Draw();
  higgs_pt_original->SetLineColor(2);
  higgs_pt_original->Draw("same");
  c->SaveAs("Masa_rekonstruiranog_higgsa.png");

  
}

void Analyzer::anti_kt_histo(string filename)
{
  TCanvas *c = new TCanvas("c","c",1000,500);

   histo_akt = new TH1F("anti_kt_alg","",20,0,10);
	
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
     

        // Read output and send it to dedicated variables
        linestream >> temp;

	histo_akt->Fill(temp);

    }
  }

 
 histo_akt->Draw();
 c->SaveAs("anti_kt.png");
}
