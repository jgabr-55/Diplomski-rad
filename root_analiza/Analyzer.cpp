#include "Analyzer.h"

void Analyzer::citanje_iz_datoteke(string filename) //pripremit 6 histograma za crtanje --> gotovo!!!!
{
  /*TCanvas *c = new TCanvas("c","c",1000,500);
  c->Divide(2,1);*/
  higgs_m_rekonstr = new TH1F("masa_rekonstr_higgsa","",32,122,128);		//od 2 b kvarka
  higgs_pt_rekonstr = new TH1F("trans_mom_rekonstr_higgsa","",60,990,1030);         //////////////////////////// ovo variramo
  higgs_rap_rekonstr = new TH1F("rapiditet_rekonstr_higgsa","",50,-10,10); // .h

  higgs_m_pyt_outpt= new TH1F("masa_original_higgsa","",32,122,128);		//iz pythia outpt direktno
  higgs_pt_pyt_outpt = new TH1F("trans_mom_original_higgsa","",60,990,1030);	////////////////////////////// ovo variramo
  higgs_rap_pyt_outpt = new TH1F("rapiditet_original_higgsa","",50,-10,10);

 			
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
	higgs_rap_rekonstr->Fill(higgs_rekonstr.Rapidity());

		
	higgs_m_pyt_outpt->Fill(higgs_pyt_outpt.M());
        higgs_pt_pyt_outpt->Fill(higgs_pyt_outpt.Pt());
	higgs_rap_pyt_outpt->Fill(higgs_pyt_outpt.Rapidity());
	/////////////////

        /*if(b.Pt()>10 && bbar.Pt()>10)
          {
		higgs_true = b + bbar; //za true vrijednost uzmemo ono što dobijemo od 2 b kvarka

  		higgs_pt_true->Fill(higgs_true.Pt());
		higgs_m_true->Fill(higgs_true.M());
	  }*/

    }
  }

  myReadFile.close();
  
  

  //cout<<particle2_en<<endl;
 /*c->cd(1);
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
  c->SaveAs("Masa_rekonstruiranog_higgsa.png");*/

  
}


void Analyzer::anti_kt_histo(string filename) // ovu ne diramo
{
  /*TCanvas *ca = new TCanvas("ca","ca",1000,500);
  ca->Divide(2,1);*/

   histo_akt_pt = new TH1F("anti_kt_alg_pt","",80,0.9,1.04);
   histo_akt_rap = new TH1F("anti_kt_alg_rap","",80,0.94,1.06);
	
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
        linestream >> particle1_px_true >> particle1_py_true >> particle1_pz_true >> particle1_en_true >> particle2_px_true >> particle2_py_true >> particle2_pz_true >> particle2_en_true >> 
		      particle1_px >> particle1_py >> particle1_pz >> particle1_en >> particle2_px >> particle2_py >> particle2_pz >> particle2_en;

	b_true.SetPxPyPzE(particle1_px_true,particle1_py_true,particle1_pz_true,particle1_en_true);
	bbar_true.SetPxPyPzE(particle2_px_true,particle2_py_true,particle2_pz_true,particle2_en_true);

	b.SetPxPyPzE(particle1_px,particle1_py,particle1_pz,particle1_en);
	bbar.SetPxPyPzE(particle2_px,particle2_py,particle2_pz,particle2_en);

	higgs_rekonstr = b + bbar;
	higgs_true = b_true + bbar_true;
	

	//histo_akt_pt->Fill(higgs_rekonstr.Pt()/higgs_true.Pt());
	histo_akt_pt -> Fill(higgs_rekonstr.Pt()/higgs_true.Pt());

	histo_akt_rap->Fill(higgs_rekonstr.Rapidity()/higgs_true.Rapidity());

    }
  }
  
   myReadFile.close();

  histo_akt_pt -> Scale(1/histo_akt_pt->Integral());
  histo_akt_rap -> Scale(1/histo_akt_rap -> Integral());

 /*ca->cd(1);
 histo_akt_pt->SetTitle("anti_kt Pt vs True Pt");
 histo_akt_pt->GetXaxis()->SetTitle("anti_kt Pt / true Pt"); 
 histo_akt_pt->Draw("histo");

 ca->cd(2);
 histo_akt_rap->SetTitle("anti_kt Rapidity vs True Rapidity");
 histo_akt_rap->GetXaxis()->SetTitle("anti_kt Rapidity / true Rapidity");
 histo_akt_rap->Draw("histo");
 ca->SaveAs("anti_kt.png");*/
}

void Analyzer::dipole_kt_histo(string filename) // ovu ne diramo
{
   /*TCanvas *cb = new TCanvas("cb","cb",1000,500);
   cb->Divide(2,1);*/

   histo_dkt_pt = new TH1F("dipole_kt_alg_pt","",80,0.9,1.04);  //0.9 - 1.04
   histo_dkt_rap = new TH1F("dipole_kt_alg_rap","",80,0.94,1.06);
	

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
        linestream >> particle1_px_true >> particle1_py_true >> particle1_pz_true >> particle1_en_true >> particle2_px_true >> particle2_py_true >> particle2_pz_true >> particle2_en_true >> 
		      particle1_px >> particle1_py >> particle1_pz >> particle1_en >> particle2_px >> particle2_py >> particle2_pz >> particle2_en;


	b_true.SetPxPyPzE(particle1_px_true,particle1_py_true,particle1_pz_true,particle1_en_true);
	bbar_true.SetPxPyPzE(particle2_px_true,particle2_py_true,particle2_pz_true,particle2_en_true);

	b.SetPxPyPzE(particle1_px,particle1_py,particle1_pz,particle1_en);
	bbar.SetPxPyPzE(particle2_px,particle2_py,particle2_pz,particle2_en);

	higgs_rekonstr = b + bbar;
	higgs_true = b_true + bbar_true;


	histo_dkt_pt->Fill(higgs_rekonstr.Pt()/higgs_true.Pt());
	//histo_dkt_pt -> Fill(bbar.Pt()/bbar_true.Pt());

	histo_dkt_rap->Fill(higgs_rekonstr.Rapidity()/higgs_true.Rapidity());
	}
   }

   myReadFile.close();

  histo_dkt_pt -> Scale(1/histo_dkt_pt->Integral());
  histo_dkt_rap -> Scale(1/histo_dkt_rap -> Integral());


/*cb->cd(1);
histo_dkt_pt->SetTitle("dipole_kt Pt vs True Pt");
histo_dkt_pt->GetXaxis()->SetTitle("dipole_kt Pt / true Pt"); 
histo_dkt_pt->Draw("histo");

cb->cd(2);
histo_dkt_rap->SetTitle("dipole_kt Rapidity vs True Rapidity");
histo_dkt_rap->GetXaxis()->SetTitle("dipole_kt Rapidity / true Rapidity");
histo_dkt_rap->Draw("histo");

cb->SaveAs("dipole_kt_rekonstr.png");*/
}


void Analyzer::kt_histo(string filename) // ovu ne diramo
{

   histo_kt_pt = new TH1F("kt_alg_pt","",80,0.9,1.04);
   histo_kt_rap = new TH1F("kt_alg_rap","",80,0.94,1.06);
	

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
        linestream >> particle1_px_true >> particle1_py_true >> particle1_pz_true >> particle1_en_true >> particle2_px_true >> particle2_py_true >> particle2_pz_true >> particle2_en_true >> 
		      particle1_px >> particle1_py >> particle1_pz >> particle1_en >> particle2_px >> particle2_py >> particle2_pz >> particle2_en;


	b_true.SetPxPyPzE(particle1_px_true,particle1_py_true,particle1_pz_true,particle1_en_true);
	bbar_true.SetPxPyPzE(particle2_px_true,particle2_py_true,particle2_pz_true,particle2_en_true);

	b.SetPxPyPzE(particle1_px,particle1_py,particle1_pz,particle1_en);
	bbar.SetPxPyPzE(particle2_px,particle2_py,particle2_pz,particle2_en);

	higgs_rekonstr = b + bbar;
	higgs_true = b_true + bbar_true;


	histo_kt_pt->Fill(higgs_rekonstr.Pt()/higgs_true.Pt());

	histo_kt_rap->Fill(higgs_rekonstr.Rapidity()/higgs_true.Rapidity());
	}
   }

   myReadFile.close();

  histo_kt_pt -> Scale(1/histo_kt_pt->Integral());
  histo_kt_rap -> Scale(1/histo_kt_rap -> Integral());


/*cb->cd(1);
histo_dkt_pt->SetTitle("dipole_kt Pt vs True Pt");
histo_dkt_pt->GetXaxis()->SetTitle("dipole_kt Pt / true Pt"); 
histo_dkt_pt->Draw("histo");

cb->cd(2);
histo_dkt_rap->SetTitle("dipole_kt Rapidity vs True Rapidity");
histo_dkt_rap->GetXaxis()->SetTitle("dipole_kt Rapidity / true Rapidity");
histo_dkt_rap->Draw("histo");

cb->SaveAs("dipole_kt_rekonstr.png");*/
}


void Analyzer::histogram_rek_jetova(string filename1, string filename2, string filename3)
{
	TCanvas *cc = new TCanvas("cc","cc",600,600);
	//cc->Divide(2,1);
	histo_dipole = new TH1F("dipole_jetovi","",4,0,4);
	histo_anti_kt = new TH1F("antikt_jetovi","",4,0,4);
	histo_kt = new TH1F("kt_jetovi","",4,0,4);

	int br1;

	ifstream myReadFile1;

 ///citamo jetove za dipole_kt alg
  myReadFile1.open(filename1.c_str());
  string line1;

  _skipFirstLine = true;

  if (myReadFile1.is_open())
  {
    // Read the file line by line
    while(getline(myReadFile1, line1))
    {
        stringstream   linestream(line1);
      	
        // Read output and send it to dedicated variables
        linestream >> br1;
	histo_dipole->Fill(br1);
}}

myReadFile1.close();


 ///citamo jetove za anti_kt alg
  myReadFile1.open(filename2.c_str());
  
  _skipFirstLine = true;

  if (myReadFile1.is_open())
  {
   
    while(getline(myReadFile1, line1))
    {
        stringstream   linestream(line1);
      	
        linestream >> br1;
	histo_anti_kt->Fill(br1);
}}

myReadFile1.close();

 ///citamo jetove za kt alg
 myReadFile1.open(filename3.c_str());
  
  _skipFirstLine = true;

  if (myReadFile1.is_open())
  {
   
    while(getline(myReadFile1, line1))
    {
        stringstream   linestream(line1);
      	
        linestream >> br1;
	histo_kt->Fill(br1);
}}

myReadFile1.close();

//cc->cd(1);
//histo_dipole->SetTitle("Broj rekonstruiranih jet-ova");
histo_dipole->GetXaxis()->SetTitle("N");
histo_dipole->SetLineWidth(2);
histo_dipole->SetLineColor(2);
histo_dipole->Draw();
//cc->cd(2);
//histo_test1->SetTitle("anti_kt alg (10 GeV)");
histo_kt->SetLineColor(210);
histo_kt->SetLineWidth(2);
histo_kt->SetLineStyle(9);
histo_kt->Draw("same");
//histo_anti_kt->SetLineColor(2);
histo_anti_kt->SetLineWidth(2);
histo_anti_kt->SetLineStyle(5);
histo_anti_kt->Draw("same");

TLegend *legend2 = new TLegend(0.1,0.8,0.35,0.9);

legend2->AddEntry(histo_dipole,"dipol-k_{t} algoritam","l");
legend2->AddEntry(histo_anti_kt,"anti-k_{t} algoitam","l");
legend2->AddEntry(histo_kt,"k_{t} algoritam","l");

legend2->Draw("same");
gStyle->SetOptStat(0);


cc->SaveAs("rekonstr_jetovi_tisucu_GeV.png"); /////////// variramo
}

void Analyzer::Crtanje()
{
	TCanvas *cf = new TCanvas("cf","cf",1600,1000);
   	cf->Divide(2,1);


cf->cd(1);

histo_dkt_pt->GetYaxis()->SetTitleOffset(1.7);
histo_dkt_pt->GetXaxis()->SetTitle("p_{#perp}  /  p_{#perp 0}"); 
 histo_dkt_pt->GetYaxis()->SetTitle("Relativna frekvencija");
 //histo_akt_pt->GetYaxis()->SetLabelSize(0.1);
 histo_dkt_pt->GetYaxis()->SetRangeUser(0,0.6); ////////////////
 if (gPad) gPad->SetLeftMargin(0.15);

 histo_dkt_pt->SetLineColor(2);
histo_dkt_pt->SetLineWidth(2);
histo_dkt_pt->Draw("histo");

histo_kt_pt->SetLineColor(210);
histo_kt_pt->SetLineWidth(2);
histo_kt_pt->SetLineStyle(9);
histo_kt_pt->Draw("histo, same");


 histo_akt_pt->SetLineWidth(2);
 histo_akt_pt->SetLineStyle(5);
 histo_akt_pt->Draw("histo, same");
 


TLegend *legend1 = new TLegend(0.15,0.8,0.45,0.9);

legend1->AddEntry(histo_akt_pt,"anti-k_{t} inkluzivni, R = 0.3","l");
legend1->AddEntry(histo_dkt_pt,"dipole-k_{t} ekskluzivni","l");
legend1->AddEntry(histo_kt_pt,"k_{t} ekskluzivni, R = 0.3","l");

legend1->Draw("same");
gStyle->SetOptStat(0);

cf->cd(2);

histo_dkt_rap->GetYaxis()->SetTitleOffset(1.7);
 histo_dkt_rap->GetXaxis()->SetTitle("y / y_{0}");
 histo_dkt_rap->SetLineWidth(2);
 histo_dkt_rap->GetYaxis()->SetTitle("Relativna frekvencija");
 histo_dkt_rap->GetYaxis()->SetRangeUser(0,0.4); //////////////////
  if (gPad) gPad->SetLeftMargin(0.15);

histo_dkt_rap->SetLineColor(2);
histo_dkt_rap->SetLineWidth(2);
histo_dkt_rap->Draw("histo");

histo_kt_rap->SetLineColor(210);
histo_kt_rap->SetLineWidth(2);
histo_kt_rap->SetLineStyle(9);
histo_kt_rap->Draw("histo, same");

//histo_akt_rap->SetTitle("anti_kt Rapidity vs dipole_kt Rapidity vs kt Rapidity");
 
 histo_akt_rap->SetLineWidth(2);
 histo_akt_rap->SetLineStyle(5);
 histo_akt_rap->Draw("histo, same");




TLegend *legend2 = new TLegend(0.15,0.8,0.45,0.9);

legend2->AddEntry(histo_akt_rap,"anti-k_{t} inkluzivni, R = 0.3","l");
legend2->AddEntry(histo_dkt_rap,"dipole-k_{t} ekskluzivni","l");
legend2->AddEntry(histo_kt_rap,"k_{t} ekskluzivni, R = 0.3","l");
legend2->Draw("same");
gStyle->SetOptStat(0);

cf->SaveAs("tisucu_Gev_Higgs_usporedba.png");   ///////////// variramo


 TCanvas *ctrue = new TCanvas ("ctrue","ctrue",2000,1000);
 
 ctrue->Divide(3,1);
 
 TLegend *mtrue = new TLegend(0.1,0.8,0.3,0.9);
 TLegend *pttrue = new TLegend(0.1,0.8,0.3,0.9);
 TLegend *raptrue = new TLegend(0.1,0.8,0.3,0.9);

 ctrue->cd(1);
 
 higgs_m_rekonstr->SetTitle("Masa");
 higgs_m_rekonstr->GetXaxis()->SetTitle("m");
 higgs_m_rekonstr->GetYaxis()->SetRangeUser(0,600); //////////
 higgs_m_rekonstr->SetLineWidth(2);
 higgs_m_rekonstr->Draw();
 higgs_m_pyt_outpt->SetLineColor(2);
 higgs_m_pyt_outpt->SetLineWidth(2);
 higgs_m_pyt_outpt->SetLineStyle(5);
 higgs_m_pyt_outpt->Draw("same");
 
 mtrue->AddEntry(higgs_m_rekonstr,"b + #bar{b}", "l");
 mtrue->AddEntry(higgs_m_pyt_outpt,"Higgs", "l");
 mtrue->Draw("same");

 ctrue->cd(2);
 higgs_pt_rekonstr->SetTitle("Transverzalni moment");
 higgs_pt_rekonstr->GetXaxis()->SetTitle("p_{#perp}");
 higgs_pt_rekonstr->SetLineWidth(2);
 higgs_pt_rekonstr->Draw();
 higgs_pt_pyt_outpt->SetLineColor(2);
 higgs_pt_pyt_outpt->SetLineWidth(2);
 higgs_pt_pyt_outpt->SetLineStyle(5);
 higgs_pt_pyt_outpt->Draw("same");

 pttrue->AddEntry(higgs_pt_rekonstr,"b + #bar{b}","l");
 pttrue->AddEntry(higgs_pt_pyt_outpt,"Higgs", "l");
 pttrue->Draw("same");

 ctrue->cd(3);
 higgs_rap_rekonstr->SetTitle("Rapiditet");
 higgs_rap_rekonstr->GetXaxis()->SetTitle("y");
 higgs_rap_rekonstr->SetLineWidth(2);
 higgs_rap_rekonstr->Draw();
 higgs_rap_pyt_outpt->SetLineColor(2);
 higgs_rap_pyt_outpt->SetLineWidth(2);
 higgs_rap_pyt_outpt->SetLineStyle(5);
 higgs_rap_pyt_outpt->Draw("same");

 raptrue->AddEntry(higgs_rap_rekonstr,"b + #bar{b}","l");
 raptrue->AddEntry(higgs_rap_pyt_outpt,"Higgs", "l");
 raptrue->Draw("same");

 ctrue->SaveAs("bbar_vs_Higgs_tisucu.png"); //////////////// variramo


 
}

