#include "Analyzer.h"

void Analyzer::citanje_iz_datoteke(string filename, bool signal) 
{

  higgs_m_rekonstr = new TH1F("masa_rekonstr_higgsa","",32,122,128);		//od 2 b kvarka
  higgs_pt_rekonstr = new TH1F("trans_mom_rekonstr_higgsa","",60,890,930);         //////////////////////////// ovo variramo
  higgs_rap_rekonstr = new TH1F("rapiditet_rekonstr_higgsa","",50,-10,10); // .h

  higgs_m_pyt_outpt= new TH1F("masa_original_higgsa","",32,122,128);		//iz pythia outpt direktno
  higgs_pt_pyt_outpt = new TH1F("trans_mom_original_higgsa","",60,890,930);	////////////////////////////// ovo variramo
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
if(signal)
{
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
}

else
{
	linestream >> particle1_px >> particle1_py >> particle1_pz >> particle1_en >> particle1_m >> particle2_px >> particle2_py >> particle2_pz >> particle2_en >> particle2_m;
	
	b.SetPxPyPzE(particle1_px,particle1_py,particle1_pz,particle1_en);
	bbar.SetPxPyPzE(particle2_px,particle2_py,particle2_pz,particle2_en);

	higgs_rekonstr = b + bbar;

	higgs_m_rekonstr->Fill(higgs_rekonstr.M());

}

    }
  }

  myReadFile.close();
  
}


void Analyzer::anti_kt_histo(string filename, bool signal) // ovu ne diramo
{
if(signal)
{
   histo_akt_pt = new TH1F("anti_kt_alg_pt","",80,0.9,1.04); //0.9 1.04
   histo_akt_rap = new TH1F("anti_kt_alg_rap","",80,-0.1,0.1);
   higgs_m_nakon_akt = new TH1F("masa_nakon_akt","",80,100,150);
}

else
{
   histo_akt_pt = new TH1F("anti_kt_alg_pt","",80,0.8,1.2);
   histo_akt_rap = new TH1F("anti_kt_alg_rap","",80,-0.1,0.2);
   histo_akt_pt1 = new TH1F("anti_kt_alg_pt1","",80,0.8,1.2);
   histo_akt_rap1 = new TH1F("anti_kt_alg_rap1","",80,-0.1,0.2);

   higgs_m_nakon_akt = new TH1F("masa_nakon_akt","",200,1000,3000);
}	
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



if(signal)
{
	higgs_rekonstr = b + bbar;
	higgs_true = b_true + bbar_true;
	
	higgs_m_nakon_akt->Fill(higgs_rekonstr.M());

	//histo_akt_pt->Fill(higgs_rekonstr.Pt()/higgs_true.Pt());
	histo_akt_pt -> Fill(higgs_rekonstr.Pt()/higgs_true.Pt());

	histo_akt_rap->Fill(fabs(higgs_rekonstr.Rapidity() - higgs_true.Rapidity()));
}
else
{

//	cout<<b_true.Pt()<<"\t"<<bbar_true.Pt()<<"\t"<<b.Pt()<<"\t"<<bbar.Pt()<<endl;
	

	if(b.Pt()>=bbar.Pt()){
		leading_jet = b;
		second_jet = bbar;
		}
	else{
		leading_jet = bbar;
		second_jet = b;
		}
		
	if(b_true.Pt()>=bbar_true.Pt()){
		leading_jet_true = b_true;
		second_jet_true = bbar_true;
		}
	else{
		leading_jet_true = bbar_true;
		second_jet_true = b_true;
		}

	histo_akt_pt->Fill(leading_jet.Pt()/leading_jet_true.Pt());
	histo_akt_pt1->Fill(second_jet.Pt()/second_jet_true.Pt());

	//cout<<leading_jet_true.Pt()<<"\t"<<second_jet_true.Pt()<<"\t"<<leading_jet.Pt()<<"\t"<<second_jet.Pt()<<endl;	

	histo_akt_rap->Fill(fabs(leading_jet.Rapidity() - leading_jet_true.Rapidity()));
	histo_akt_rap1->Fill(fabs(second_jet.Rapidity() - second_jet_true.Rapidity()));

	zbroj_pozadina = leading_jet + second_jet;

	higgs_m_nakon_akt->Fill(zbroj_pozadina.M());

}

    }
  }
  
   myReadFile.close();

if(signal)
{
  histo_akt_pt -> Scale(1/histo_akt_pt->Integral());
  histo_akt_rap -> Scale(1/histo_akt_rap -> Integral());
  higgs_m_nakon_akt -> Scale(1/higgs_m_nakon_akt -> Integral());
}

else
{
	histo_akt_pt->Scale(1/histo_akt_pt->Integral());
	histo_akt_pt1->Scale(1/histo_akt_pt1->Integral());
	histo_akt_rap->Scale(1/histo_akt_rap->Integral());
	histo_akt_rap1->Scale(1/histo_akt_rap1->Integral());
  	higgs_m_nakon_akt -> Scale(1/higgs_m_nakon_akt -> Integral());
}

}

void Analyzer::dipole_kt_histo(string filename, bool signal) // ovu ne diramo
{
int pt_lead_out = 0, pt_sec_out = 0, rap_lead_out = 0, rap_sec_out = 0;
if(signal)
{
   histo_dkt_pt = new TH1F("dipole_kt_alg_pt","",80,0.9,1.04);  //0.9 - 1.04
   histo_dkt_rap = new TH1F("dipole_kt_alg_rap","",80,-0.1,0.1);
   higgs_m_nakon_dkt = new TH1F("masa_nakon_dkt","",80,100,150);
}
else
{
   histo_dkt_pt = new TH1F("dipole_kt_alg_pt","",100,0,2);  //0.9 - 1.04
   histo_dkt_rap = new TH1F("dipole_kt_alg_rap","",105,-3,3);	
   histo_dkt_pt1 = new TH1F("dipole_kt_alg_pt1","",100,0,2);
   histo_dkt_rap1 = new TH1F("dipole_kt_alg_rap1","",105,-3,3);

   higgs_m_nakon_dkt = new TH1F("masa_nakon_dkt","",200,1000,3000);
}

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


	b_true.SetPxPyPzE(particle1_px_true,particle1_py_true,particle1_pz_true,particle1_en_true); // true 1 gluon
	bbar_true.SetPxPyPzE(particle2_px_true,particle2_py_true,particle2_pz_true,particle2_en_true); // true 2. gluon

	b.SetPxPyPzE(particle1_px,particle1_py,particle1_pz,particle1_en); // 1. rek gluon
	bbar.SetPxPyPzE(particle2_px,particle2_py,particle2_pz,particle2_en); // 2. rek gluon

	//cout<<b_true.Rapidity()<<"\t"<<bbar_true.Rapidity()<<endl;   //rapiditeti razliciti

if(signal)
{
	higgs_rekonstr = b + bbar;
	higgs_true = b_true + bbar_true;

	higgs_m_nakon_dkt->Fill(higgs_rekonstr.M());

	histo_dkt_pt->Fill(higgs_rekonstr.Pt()/higgs_true.Pt());
	//histo_dkt_pt -> Fill(bbar.Pt()/bbar_true.Pt());

	histo_dkt_rap->Fill(fabs(higgs_rekonstr.Rapidity() - higgs_true.Rapidity()));
}
else
{
	if(b.Pt()>=bbar.Pt()){
		leading_jet = b;
		second_jet = bbar;
		}
	else{
		leading_jet = bbar;
		second_jet = b;
		}
		
	if(b_true.Pt()>=bbar_true.Pt()){
		leading_jet_true = b_true;
		second_jet_true = bbar_true;
		}
	else{
		leading_jet_true = bbar_true;
		second_jet_true = b_true;
		}

	histo_dkt_pt->Fill(leading_jet.Pt()/leading_jet_true.Pt());
	histo_dkt_pt1->Fill(second_jet.Pt()/second_jet_true.Pt());
	
	if((leading_jet.Pt()/leading_jet_true.Pt())>2.0 || (leading_jet.Pt()/leading_jet_true.Pt())<0.0)
		pt_lead_out++;
	if((second_jet.Pt()/second_jet_true.Pt())>2.0 || (second_jet.Pt()/second_jet_true.Pt())<0.0)
		pt_sec_out++;
	


	if(fabs(b.Rapidity())<=fabs(bbar.Rapidity()))
		{
		leading_jet = b;
		second_jet = bbar;
		}
	else{
		leading_jet = bbar;
		second_jet = b;
		}
		
	if(fabs(b_true.Rapidity())<=fabs(bbar_true.Rapidity())){
		leading_jet_true = b_true;
		second_jet_true = bbar_true;
		}
	else{
		leading_jet_true = bbar_true;
		second_jet_true = b_true;
		}

	histo_dkt_rap->Fill(fabs(leading_jet.Rapidity()) - fabs(leading_jet_true.Rapidity()));
	histo_dkt_rap1->Fill(fabs(second_jet.Rapidity()) - fabs(second_jet_true.Rapidity()));

	if((fabs(leading_jet.Rapidity()) - fabs(leading_jet_true.Rapidity()))>3.0 || (fabs(leading_jet.Rapidity()) - fabs(leading_jet_true.Rapidity()))<-3.0)
		rap_lead_out++;
		//cout<<fabs(leading_jet.Rapidity()) - fabs(leading_jet_true.Rapidity())<<endl;}
	if((fabs(second_jet.Rapidity()) - fabs(second_jet_true.Rapidity()))>3.0 ||(fabs(second_jet.Rapidity()) - fabs(second_jet_true.Rapidity()))<-3.0)
		rap_sec_out++;


	zbroj_pozadina = leading_jet + second_jet;

	higgs_m_nakon_dkt->Fill(zbroj_pozadina.M());
}
	}
   }

   myReadFile.close();

if(signal)
{
  histo_dkt_pt -> Scale(1/histo_dkt_pt->Integral());
  histo_dkt_rap -> Scale(1/histo_dkt_rap -> Integral());
  higgs_m_nakon_dkt -> Scale(1/higgs_m_nakon_dkt -> Integral());
}
else
{
	/*histo_dkt_pt->Scale(1/histo_dkt_pt->Integral());
	histo_dkt_pt1->Scale(1/histo_dkt_pt1->Integral());
	histo_dkt_rap->Scale(1/histo_dkt_rap->Integral());
	histo_dkt_rap1->Scale(1/histo_dkt_rap1->Integral());

	higgs_m_nakon_dkt -> Scale(1/higgs_m_nakon_dkt -> Integral());*/
	cout<<pt_lead_out<<"\t"<<pt_sec_out<<"\t"<<rap_lead_out<<"\t"<<rap_sec_out<<endl;
}

}


void Analyzer::kt_histo(string filename, bool signal) // ovu ne diramo
{
int pt_lead_out = 0, pt_sec_out = 0, rap_lead_out = 0, rap_sec_out = 0;
if(signal)
{
   histo_kt_pt = new TH1F("kt_alg_pt","",80,0.9,1.04);
   histo_kt_rap = new TH1F("kt_alg_rap","",80,-0.1,0.1);
   higgs_m_nakon_kt = new TH1F("masa_nakon_kt","",80,100,150);
}

else
{
   histo_kt_pt = new TH1F("kt_alg_pt","",100,0,2);
   histo_kt_rap = new TH1F("kt_alg_rap","",105,-3,3);	
   histo_kt_pt1 = new TH1F("kt_alg_pt1","",100,0,2);
   histo_kt_rap1 = new TH1F("kt_alg_rap1","",105,-3,3);

   higgs_m_nakon_kt = new TH1F("masa_nakon_kt","",200,1000,3000);
}

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


	b_true.SetPxPyPzE(particle1_px_true,particle1_py_true,particle1_pz_true,particle1_en_true); // 1. gluon true
	bbar_true.SetPxPyPzE(particle2_px_true,particle2_py_true,particle2_pz_true,particle2_en_true); //  2. gluon true

	b.SetPxPyPzE(particle1_px,particle1_py,particle1_pz,particle1_en); // 1. rek gluon
	bbar.SetPxPyPzE(particle2_px,particle2_py,particle2_pz,particle2_en); // 2. rek gluon

if(signal)
{
	higgs_rekonstr = b + bbar;
	higgs_true = b_true + bbar_true;


	higgs_m_nakon_kt->Fill(higgs_rekonstr.M());

	histo_kt_pt->Fill(higgs_rekonstr.Pt()/higgs_true.Pt());

	histo_kt_rap->Fill(fabs(higgs_rekonstr.Rapidity() - higgs_true.Rapidity()));
}
else
{
	if(b.Pt()>=bbar.Pt()){
		leading_jet = b;
		second_jet = bbar;
		}
	else{
		leading_jet = bbar;
		second_jet = b;
		}
		
	if(b_true.Pt()>=bbar_true.Pt()){
		leading_jet_true = b_true;
		second_jet_true = bbar_true;
		}
	else{
		leading_jet_true = bbar_true;
		second_jet_true = b_true;
		}

	histo_kt_pt->Fill(leading_jet.Pt()/leading_jet_true.Pt());
	histo_kt_pt1->Fill(second_jet.Pt()/second_jet_true.Pt());
	
	if((leading_jet.Pt()/leading_jet_true.Pt())>2.0 || (leading_jet.Pt()/leading_jet_true.Pt())<0.0)
		pt_lead_out++;
	if((second_jet.Pt()/second_jet_true.Pt())>2.0 || (second_jet.Pt()/second_jet_true.Pt())<0.0)
		pt_sec_out++;
	
	
	if(fabs(b.Rapidity())<=fabs(bbar.Rapidity()))
		{
		leading_jet = b;
		second_jet = bbar;
		}
	else{
		leading_jet = bbar;
		second_jet = b;
		}
		
	if(fabs(b_true.Rapidity())<=fabs(bbar_true.Rapidity())){
		leading_jet_true = b_true;
		second_jet_true = bbar_true;
		}
	else{
		leading_jet_true = bbar_true;
		second_jet_true = b_true;
		}

	histo_kt_rap->Fill(fabs(leading_jet.Rapidity()) - fabs(leading_jet_true.Rapidity()));
	histo_kt_rap1->Fill(fabs(second_jet.Rapidity()) - fabs(second_jet_true.Rapidity()));

	if((fabs(leading_jet.Rapidity()) - fabs(leading_jet_true.Rapidity()))>3.0 || (fabs(leading_jet.Rapidity()) - fabs(leading_jet_true.Rapidity()))<-3.0)
		rap_lead_out++;
		//cout<<fabs(leading_jet.Rapidity()) - fabs(leading_jet_true.Rapidity())<<endl;}
	if((fabs(second_jet.Rapidity()) - fabs(second_jet_true.Rapidity()))>3.0 ||(fabs(second_jet.Rapidity()) - fabs(second_jet_true.Rapidity()))<-3.0)
		rap_sec_out++;


	zbroj_pozadina = leading_jet + second_jet;

        higgs_m_nakon_kt->Fill(zbroj_pozadina.M());
}
	}
   }

   myReadFile.close();

if(signal)
{
  histo_kt_pt -> Scale(1/histo_kt_pt->Integral());
  histo_kt_rap -> Scale(1/histo_kt_rap -> Integral());
  higgs_m_nakon_kt -> Scale(1/higgs_m_nakon_kt -> Integral());
}

else
{
	/*histo_kt_pt->Scale(1/histo_kt_pt->Integral());
	histo_kt_pt1->Scale(1/histo_kt_pt1->Integral());
	histo_kt_rap->Scale(1/histo_kt_rap->Integral());
	histo_kt_rap1->Scale(1/histo_kt_rap1->Integral());

	higgs_m_nakon_kt -> Scale(1/higgs_m_nakon_kt -> Integral());*/

	cout<<pt_lead_out<<"\t"<<pt_sec_out<<"\t"<<rap_lead_out<<"\t"<<rap_sec_out<<endl;
}


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


cc->SaveAs("rekonstr_jetovi_900_GeV_R_03_hadr_off.png"); /////////// variramo
}

void Analyzer::Crtanje_signal()
{

///////// CRTAMO PT/PT_TRUE I RAP - RAP_TRUE ZA SVA 3 ALGORITMA /////////////

TCanvas *cf = new TCanvas("cf","cf",1600,1000);
cf->Divide(2,1);


cf->cd(1);
histo_dkt_pt->SetTitle("R_08_1000_GeV_hadr___sa_leptonima");
histo_dkt_pt->GetYaxis()->SetTitleOffset(1.7);
histo_dkt_pt->GetXaxis()->SetTitle("p_{#perp}  /  p_{#perp 0}"); 
histo_dkt_pt->GetYaxis()->SetTitle("Relativna frekvencija");
//histo_akt_pt->GetYaxis()->SetLabelSize(0.1);
//histo_dkt_pt->GetYaxis()->SetRangeUser(0,0.2); ////////////////
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
histo_dkt_rap->GetXaxis()->SetTitle("|y - y_{0}|");
histo_dkt_rap->SetLineWidth(2);
histo_dkt_rap->GetYaxis()->SetTitle("Relativna frekvencija");
histo_dkt_rap->GetYaxis()->SetRangeUser(0,0.85); //////////////////
if (gPad) gPad->SetLeftMargin(0.15);

histo_dkt_rap->SetLineColor(2);
histo_dkt_rap->SetLineWidth(2);
histo_dkt_rap->Draw("histo");

histo_kt_rap->SetLineColor(210);
histo_kt_rap->SetLineWidth(2);
histo_kt_rap->SetLineStyle(9);
histo_kt_rap->Draw("histo, same");

histo_akt_rap->SetLineWidth(2);
histo_akt_rap->SetLineStyle(5);
histo_akt_rap->Draw("histo, same");


TLegend *legend2 = new TLegend(0.15,0.8,0.45,0.9);

legend2->AddEntry(histo_akt_rap,"anti-k_{t} inkluzivni, R = 0.3","l");
legend2->AddEntry(histo_dkt_rap,"dipole-k_{t} ekskluzivni","l");
legend2->AddEntry(histo_kt_rap,"k_{t} ekskluzivni, R = 0.3","l");
legend2->Draw("same");
gStyle->SetOptStat(0);

cf->SaveAs("900_Gev_Higgs_usporedba_R_03_hadr_off.png");   ///////////// variramo


////////////// CRTAMO MASU, PT, RAPIDITET ZA PODATKE IZ PYTHIA OUTPUTA /////////////////

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
 raptrue->Draw("histo, same");

 ctrue->SaveAs("bbar_vs_Higgs_900_R_03_hadr_off.png"); //////////////// variramo



 ///////// CRTAMO MASU B + BBAR NAKON REKONSTRUKCIJE ZA SVA 3 ALGORITMA ///////////////

 TCanvas *cm = new TCanvas("cm","cm",600,600);
 higgs_m_nakon_dkt->SetTitle("R_08_1000_GeV_hadr___sa_leptonima");
 higgs_m_nakon_dkt->GetXaxis()->SetTitle("m");
 higgs_m_nakon_dkt->GetYaxis()->SetTitle("Relativna frekvencija");
 if (gPad) gPad->SetLeftMargin(0.15);
 higgs_m_nakon_dkt->GetYaxis()->SetRangeUser(0,0.35);
 higgs_m_nakon_dkt->SetLineColor(2);
  higgs_m_nakon_dkt->SetLineWidth(2);
 higgs_m_nakon_dkt->Draw("histo");

 higgs_m_nakon_kt->SetLineColor(209);
 higgs_m_nakon_kt->SetLineStyle(9);
 higgs_m_nakon_kt->SetLineWidth(2);
  higgs_m_nakon_kt->Draw("histo, same");

  higgs_m_nakon_akt->SetLineWidth(2);
  higgs_m_nakon_akt->SetLineStyle(5);
  higgs_m_nakon_akt->Draw("histo, same");

 TLegend *m_nakon_alg = new TLegend(0.15,0.8,0.45,0.9);

 m_nakon_alg->AddEntry(higgs_m_nakon_dkt,"dipole k_{t}","l");
 m_nakon_alg->AddEntry(higgs_m_nakon_akt,"anti k_{t}","l");
 m_nakon_alg->AddEntry(higgs_m_nakon_kt,"k_{t}","l");

 m_nakon_alg->Draw("same");
 
 

 cm->SaveAs("masa_nakon_rekonstr_900_R_03_hadr_off.png");


 
}

void Analyzer::Crtanje_pozadina()
{

TCanvas *c_poz_pt_rap = new TCanvas("c_poz_pt_rap","c_poz_pt_rap",1000,900);
c_poz_pt_rap->Divide(2,2);

//crtamo za leading jet pt
c_poz_pt_rap->cd(1);

histo_kt_pt->SetTitle("1st jet pt");
histo_kt_pt->GetYaxis()->SetRangeUser(0,100);
histo_kt_pt->Draw("histo");
//histo_dkt_pt->SetLineWidth(2);
histo_dkt_pt->SetLineColor(2);
histo_dkt_pt->SetLineStyle(5);
histo_dkt_pt->Draw("histo, same");

//histo_kt_pt->SetLineWidth(2);
//histo_kt_pt->SetLineColor(209);
//histo_kt_pt->SetLineStyle(9);


/*histo_akt_pt->SetLineWidth(2);
histo_akt_pt->SetLineStyle(5);
histo_akt_pt->Draw("histo, same");*/

//crtamo za second jet pt
c_poz_pt_rap->cd(2);

histo_kt_pt1->SetTitle("2nd jet pt");
histo_kt_pt1->GetYaxis()->SetRangeUser(0,100);
histo_kt_pt1->Draw("histo,");
//histo_dkt_pt1->SetLineWidth(2);
histo_dkt_pt1->SetLineColor(2);
histo_dkt_pt1->SetLineStyle(5);
histo_dkt_pt1->Draw("histo, same");

//histo_kt_pt1->SetLineWidth(2);
//histo_kt_pt1->SetLineColor(209);
//histo_kt_pt1->SetLineStyle(9);


/*histo_akt_pt1->SetLineWidth(2);
histo_akt_pt1->SetLineStyle(5);
histo_akt_pt1->Draw("histo, same");*/

//crtamo za leading rapiditet
c_poz_pt_rap->cd(3);

histo_kt_rap->SetTitle("1st jet rap");
histo_kt_rap->GetYaxis()->SetRangeUser(0,150);
histo_kt_rap->Draw("histo");
//histo_dkt_rap->SetLineWidth(2);
histo_dkt_rap->SetLineColor(2);
histo_dkt_rap->SetLineStyle(5);
histo_dkt_rap->Draw("histo, same");

//histo_kt_rap->SetLineWidth(2);
//histo_kt_rap->SetLineColor(209);
//histo_kt_rap->SetLineStyle(9);


/*histo_akt_rap->SetLineWidth(2);
histo_akt_rap->SetLineStyle(5);
histo_akt_rap->Draw("histo, same");*/

//crtamo za second rapiditet
c_poz_pt_rap->cd(4);

histo_kt_rap1->SetTitle("2nd jet rap");
histo_kt_rap1->GetYaxis()->SetRangeUser(0,150);
histo_kt_rap1->Draw("histo");
//histo_dkt_rap1->SetLineWidth(2);
histo_dkt_rap1->SetLineColor(2);
histo_dkt_rap1->SetLineStyle(5);
histo_dkt_rap1->Draw("histo, same");

//histo_kt_rap1->SetLineWidth(2);
//histo_kt_rap1->SetLineColor(209);
//histo_kt_rap1->SetLineStyle(9);


/*histo_akt_rap1->SetLineWidth(2);
histo_akt_rap1->SetLineStyle(5);
histo_akt_rap1->Draw("histo, same");*/


/*TLegend *pt_rap_poz = new TLegend(0.1,0.8,0.3,0.9);

pt_rap_poz->AddEntry(histo_dkt_pt,"dipole_kt b","l");
pt_rap_poz->AddEntry(histo_dkt_pt1,"dipole_kt bbar","l");
pt_rap_poz->AddEntry(histo_kt_pt,"kt b","l");
pt_rap_poz->AddEntry(histo_kt_pt1,"kt bbar","l");
pt_rap_poz->AddEntry(histo_akt_pt,"anti_kt b","l");
pt_rap_poz->AddEntry(histo_akt_pt1,"anti_kt bbar","l");

pt_rap_poz->Draw("same");*/


c_poz_pt_rap->SaveAs("100_pozadina_pt_rap_gluoni.png");

////////// CRTAMO ZBROJ MASA JETOVA ZA POZADINU ////////////	

TCanvas *c_masa_poz = new TCanvas("c_masa_poz","c_masa_poz", 700, 700);
c_masa_poz->cd(1);

  higgs_m_nakon_dkt->SetTitle("Masa za pozadinu");
  higgs_m_nakon_dkt->GetXaxis()->SetTitle("m");
  if (gPad) gPad->SetLeftMargin(0.15);
  higgs_m_nakon_dkt->GetYaxis()->SetTitle("Relativna frekvencija");
  higgs_m_nakon_dkt->SetLineWidth(2);
  higgs_m_nakon_dkt->SetLineColor(2);
  higgs_m_nakon_dkt->Draw("histo");  

  higgs_m_nakon_kt->SetLineWidth(2);
  higgs_m_nakon_kt->SetLineColor(209);
  higgs_m_nakon_kt->SetLineStyle(9);
  higgs_m_nakon_kt->Draw("histo, same");

  /*higgs_m_nakon_akt->SetLineWidth(2);
  higgs_m_nakon_akt->SetLineStyle(5);
  higgs_m_nakon_akt->Draw("histo, same");*/

TLegend *masa_poz = new TLegend(0.15,0.8,0.35,0.9);

masa_poz->AddEntry(higgs_m_nakon_dkt,"dipole_kt","l");
masa_poz->AddEntry(higgs_m_nakon_kt,"kt","l");
//masa_poz->AddEntry(higgs_m_nakon_akt,"anti_kt","l");

masa_poz->Draw("same");

c_masa_poz->SaveAs("ukupna_masa_pozadina.png");



}

void Analyzer::ProvjeraKoda(string filename_dkt, string filename_kt)
{
	histo_dkt_pt = new TH1F("dipole_kt_alg_pt","",120,0,2);  //0.9 - 1.04
   	histo_dkt_rap = new TH1F("dipole_kt_alg_rap","",120,-3,3);

	histo_kt_pt = new TH1F("kt_alg_pt","",120,0,2);
        histo_kt_rap = new TH1F("kt_alg_rap","",120,-3,3);

	int br1 = 0, br2 = 0;

	

  ifstream myReadFile;
  myReadFile.open(filename_dkt.c_str());
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
        linestream >> particle1_px_true >> particle1_py_true >> particle1_pz_true >> particle1_en_true >> 
		      particle1_px >> particle1_py >> particle1_pz >> particle1_en;


	b_true.SetPxPyPzE(particle1_px_true,particle1_py_true,particle1_pz_true,particle1_en_true);

	b.SetPxPyPzE(particle1_px,particle1_py,particle1_pz,particle1_en);
	
	histo_dkt_pt->Fill(b.Pt()/b_true.Pt());

	if(b.Pt()/b_true.Pt() > 2.0 || b.Pt()/b_true.Pt()<0.0)
		br1++;

	histo_dkt_rap->Fill(fabs(b.Rapidity())-fabs(b_true.Rapidity()));
	
	/*if((b.Pt()/b_true.Pt())>10)
	{
		cout<<b.Pt()/b_true.Pt()<<"\t"<<b.Pt()<<"\t"<<b_true.Pt()<<endl;
	}*/
   }}

 myReadFile.close();



  //ifstream myReadFile;
  myReadFile.open(filename_kt.c_str());
  //string line;

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
        linestream >> particle2_px_true >> particle2_py_true >> particle2_pz_true >> particle2_en_true >> 
		      particle2_px >> particle2_py >> particle2_pz >> particle2_en;


	b_true.SetPxPyPzE(particle2_px_true,particle2_py_true,particle2_pz_true,particle2_en_true);

	b.SetPxPyPzE(particle2_px,particle2_py,particle2_pz,particle2_en);
	
	histo_kt_pt->Fill(b.Pt()/b_true.Pt());

	if(b.Pt()/b_true.Pt() > 2.0 || b.Pt()/b_true.Pt()<0.0)
		br2++;


	histo_kt_rap->Fill(fabs(b.Rapidity())-fabs(b_true.Rapidity()));

	//if((b.Pt()/b_true.Pt())>10)
	//{
		//cout<<b.Pt()/b_true.Pt()<<"\t"<<b.Pt()<<"\t"<<b_true.Pt()<<endl;
	//}
   }}

 myReadFile.close();

 cout<<br1<<"\t"<<br2<<endl;
 
 TCanvas *c_test = new TCanvas("c_test","c_test", 1600, 800);
 c_test->Divide(2,1);

 c_test->cd(1);

 histo_dkt_pt->GetXaxis()->SetTitle("pt/pt_true");
 histo_dkt_pt->SetLineWidth(2);
 histo_dkt_pt->SetLineColor(2);
 histo_dkt_pt->SetLineStyle(9);
 histo_dkt_pt->Draw();
 histo_kt_pt->Draw("same");

 c_test->cd(2);

 histo_dkt_rap->GetXaxis()->SetTitle("y - y_true");
 histo_dkt_rap->SetLineWidth(2);
 histo_dkt_rap->SetLineColor(2);
 histo_dkt_rap->SetLineStyle(9);
 histo_dkt_rap->Draw();
 histo_kt_rap->Draw("same");

c_test->SaveAs("1_jet_test_01.png");
 


}
