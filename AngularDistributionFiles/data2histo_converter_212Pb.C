{

#include "TROOT.h"
#include "TStyle.h"
#include "../Analysis/LibPerson.C"
#include "TVector3.h"

gROOT->SetStyle("Plain");
gStyle->SetOptStat(0);

const double R2D = 57.29577951;
const double D2R = 0.017453292;
const double PI = 3.141592653589793238;
const double Mp = 0.938279; // GeV
const double Mn = 0.939565; // GeV

TH1F *hdsigma1 = new TH1F("hdsigma1","hdsigma1",200,-0.5,200);
TH2F *hdsigma1VsTheta = new TH2F("hdsigma1VsTheta","d#sigma / d#omega vs #theta for 212Pb(d,p)213Pb ",180,0,180,1000,0.,2);
TH2F *hdsigma2VsTheta = new TH2F("hdsigma2VsTheta","d#sigma_{2} / d#omega vs #theta ",180,0,180,1000,0.,2);
TH2F *hdsigma3VsTheta = new TH2F("hdsigma3VsTheta","d#sigma_{3} / d#omega vs #theta ",180,0,180,1000,0.,2);
TH2F *hdsigma4VsTheta = new TH2F("hdsigma4VsTheta","d#sigma_{4} / d#omega vs #theta ",180,0,180,1000,0.,2);

//*************************************** Read Momenta from the file ******************************************************************//

Int_t counter = 0;
Float_t theta[200], sigma1[200], sigma2[200], sigma3[200], sigma4[200],dummy=0.0;

ifstream infile_1;

infile_1.open("212Pb_dp.dat"); 
//infile_1.open("Mg30PereyFrescoSpline_p.dat"); 

while(!infile_1.eof())
{
  
     infile_1 >> theta[counter] >> sigma1[counter] >> dummy >> sigma2[counter] >> dummy >> sigma3[counter] >> dummy >> sigma4[counter];
      //infile_1 >> theta[counter] >> sigma1[counter] >> dummy >> dummy >> dummy >> dummy;
 counter++;
            //if(counter>0) cout<<"counter: "<<counter<< "\t" <<theta[counter-1]<< "\t" <<sigma1[counter-1]<< "\t" <<sigma2[counter-1]<< "\t" <<sigma3[counter-1]<< "\t" <<sigma4[counter-1]<<endl;

              //if(counter==5) break;

}

	infile_1.close();

	for(int i=1;i<counter-1;i++){

		cout<<" i: "<<i<<"\t" <<theta[i]<< "\t" <<sigma1[i]<< "\t" <<sigma2[i]<< "\t" <<sigma3[i]<< "\t" <<sigma4[i]<<endl;                   
		hdsigma1VsTheta->Fill(theta[i],sigma1[i]);
		hdsigma2VsTheta->Fill(theta[i],sigma2[i]);
		hdsigma3VsTheta->Fill(theta[i],sigma3[i]);
		hdsigma4VsTheta->Fill(theta[i],sigma4[i]);
	}
		
  TCanvas *c1 = new TCanvas("c1","z",800,700);
  c1->Divide(1,1);
  c1->cd(1);
  //c1->cd(1)->SetGridx();
  //c1->cd(1)->SetGridy();
  c1->cd(1)->SetTickx();
  c1->cd(1)->SetTicky();  
  //c1->cd(1)->SetLogy(); 
  hdsigma1VsTheta->Draw(""); 
  hdsigma1VsTheta->SetMarkerStyle(4); 
  hdsigma1VsTheta->SetMarkerColor(2); 
  hdsigma2VsTheta->Draw("same"); 
  hdsigma2VsTheta->SetMarkerStyle(4); 
  hdsigma2VsTheta->SetMarkerColor(3);
  hdsigma3VsTheta->Draw("same"); 
  hdsigma3VsTheta->SetMarkerStyle(4); 
  hdsigma3VsTheta->SetMarkerColor(1);
  hdsigma4VsTheta->Draw("same"); 
  hdsigma4VsTheta->SetMarkerStyle(4); 
  hdsigma4VsTheta->SetMarkerColor(4);
  
  TLegend * L_all = new TLegend(0.40, 0.60, 0.85, 0.85);
  L_all->AddEntry(hdsigma1VsTheta,"g.s. l=0","lp");
  L_all->AddEntry(hdsigma2VsTheta,"805 keV 3/2+ state l=2","lp");
  L_all->AddEntry(hdsigma3VsTheta,"1117 keV 3/2- state l=1","lp");
  L_all->AddEntry(hdsigma4VsTheta,"1277 keV 3/2- state l=1","lp");
  L_all->SetFillColor(kWhite);
  L_all->SetTextSize(0.05);
  L_all->SetTextFont(132);
  L_all->SetBorderSize(0);
  L_all->Draw();
  
//###################################################################

  // store
  /*c1->Print("simtohist_converter_beampos.pdf(");
  c2->Print("simtohist_converter_beampos.pdf");
  c6->Print("simtohist_converter_beampos.pdf)");*/

  c1->Print("30Mg_dp_angularDistribution.pdf");

  //gSystem->Exec("mv 30Mg_dp_angularDistribution.pdf /home/latar/devinTrex/TRexGeant4/dropbox");

}
