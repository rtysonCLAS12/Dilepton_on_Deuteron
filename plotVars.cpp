#include <cstdlib>
#include <iostream>
#include <TROOT.h>
#include <TH1.h>
#include <TChain.h>
#include <TCanvas.h>
#include "TH2.h"
#include <cstring>
#include <TLatex.h>
#include <TDatabasePDG.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <fstream>
#include <cmath>
#include <chrono>
#include <TBenchmark.h>

using namespace std;

void plotVars(){

  gStyle->SetOptStat("");
  gStyle->SetOptFit(0012);

  string fileLoc = "/w/work/clas12/tyson/plots/c12_scripts/eed/IM/";

  string treeLocRoot = "/w/work/clas12/tyson/data_repo/c12scripts_out/eed/eedFS_allRGB.root";

  string endName="_allRGB"; //_noCuts

  auto file = new TFile((treeLocRoot).c_str());
  //auto output =(TTree*) file->Get("FINALOUTTREE");
  auto output =(TTree*) file->Get("eed");

  gStyle->SetOptStat("");

  string cut="";

  cut+="abs(elStatus)>=2000 && abs(elStatus)<4000 && abs(poStatus)>=2000 && abs(poStatus)<4000";
  cut+="&& Q2<5 && abs(MM2)<1";
  cut+="&& deutChi2PID<5";

  string cutFD="&& abs(deutStatus)>=2000 && abs(deutStatus)<4000";
  string cutCD="&& abs(deutStatus)>=4000";

  string cut_tighter=cut;
  cut_tighter+="&& elTriangCut==1 && poTriangCut==1"; //e+ e- ID
  cut_tighter+="&& Combis==0";
  
 
  string cutIM=cut_tighter+"&& abs(MM2)<1. && Q2<1.0";

  TF1* f2 = new TF1("Polynomial Background and Gaussian Signal2","[0]*0.398942*0.0333*TMath::Exp(-0.5*((x-[1])/([2]))*((x-[1])/([2])))/TMath::Abs([2]) -[3]*(x-[1]) - [4]*(x-[1])*(x-[1]) - [5]*(x-[1])*(x-[1])*(x-[1]) + [6]");
  f2->SetParameters(1, 3.097, 1, 1, 1,1,1);
  f2->SetParNames("J/#psi Yield", "Mean", "#sigma", "1st order coef", "2nd order coef", "3rd order coef", "offset");
  f2->SetLineColor(kBlack);
  f2->SetLineWidth(2);
  f2->ReleaseParameter(1);
  f2->SetParLimits(1, 3.06, 3.11);
  f2->ReleaseParameter(2);
  f2->SetParLimits(2, 0.02, 0.1);

  TF1* gauss2 = new TF1("Gaussian Signal2","[0]*0.398942*0.0333*TMath::Exp(-0.5*((x-[1])/([2]))*((x-[1])/([2])))/TMath::Abs([2])");
  gauss2->SetLineColor(kRed); 
  gauss2->SetRange(2.5, 3.5);
  gauss2->SetLineStyle(2);
  gauss2->SetLineWidth(2);
  //gauss->SetParameters(f->GetParameter(0),f->GetParameter(1),f->GetParameter(2));


  TF1* bg2 = new TF1("Polynomial Background2","-[1]*(x-[0]) - [2]*(x-[0])*(x-[0]) - [3]*(x-[0])*(x-[0])*(x-[0]) + [4]");
  bg2->SetRange(2.5, 3.5);
  bg2->SetLineStyle(2);
  bg2->SetLineWidth(2);
  bg2->SetLineColor(kBlack);

  TCanvas cIMFit;
  TH1F *hIMFit=new TH1F("hIMFit","e+ e- Invariant Mass",30,2.5,3.5);
  hIMFit->SetTitle("e+ e- Invariant Mass ");
  hIMFit->GetXaxis()->SetTitle("Invariant Mass [GeV]");
  hIMFit->SetLineWidth(2);
  output->Draw("IM>>hIMFit", cutIM.c_str(),"");
  hIMFit->Fit(f2, "", "", 2.5, 3.5);
  hIMFit->SetFillColor(kAzure-9);
  gauss2->SetParameters(f2->GetParameter(0),f2->GetParameter(1),f2->GetParameter(2));
  bg2->SetParameters(f2->GetParameter(1),f2->GetParameter(3),f2->GetParameter(4),f2->GetParameter(5),f2->GetParameter(6));
  hIMFit->Draw();
  gauss2->Draw("same");
  bg2->Draw("same");
  cIMFit.Draw();
  cIMFit.SaveAs((fileLoc+"Vars"+endName+".pdf(").c_str());

  TCanvas cEgamma;
  TH1F *hEgamma=new TH1F("hEgamma","E_{#gamma}",60,5,11);
  hEgamma->GetXaxis()->SetTitle("E_{#gamma} [GeV]");
  output->Draw("Egamma>>hEgamma", (cut).c_str());
  cEgamma.Draw();
  cEgamma.SaveAs((fileLoc+"Vars"+endName+".pdf").c_str());

  TCanvas cMM2;
  TH1F *hMM2=new TH1F("hMM2","Missing Mass Squared",100,-1,1);
  hMM2->GetXaxis()->SetTitle("Missing Mass Squared [GeV^{2}]");
  output->Draw("MM2>>hMM2", (cut).c_str());
  cMM2.Draw();
  cMM2.SaveAs((fileLoc+"Vars"+endName+".pdf").c_str());

  TCanvas cQ2;
  TH1F *hQ2=new TH1F("hQ2","Q^{2}",100,0,0.5);
  hQ2->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
  output->Draw("Q2>>hQ2",(cut).c_str());
  cQ2.Draw();
  cQ2.SaveAs((fileLoc+"Vars"+endName+".pdf").c_str());

  TCanvas cIM;
  TH1F *hIM=new TH1F("hIM","e^{+} e^{-} Invariant Mass",100,0,3.5);
  hIM->GetXaxis()->SetTitle("Invariant Mass [GeV]");
  output->Draw("IM>>hIM", (cutIM+"").c_str());
  cIM.Draw();
  cIM.SaveAs((fileLoc+"Vars"+endName+".pdf").c_str());

  TCanvas cCombis;
  TH1F *hCombis=new TH1F("hCombis","Number of Permutations per Event",11,-0.5,10.5);
  hCombis->GetXaxis()->SetTitle("Number of Permutations per Event");
  output->Draw("Combis>>hCombis", (cut+"").c_str());
  cCombis.SetLogy();
  cCombis.Draw();
  cCombis.SaveAs((fileLoc+"Vars"+endName+".pdf").c_str());

  string titlePart[3];
  titlePart[0]="e-";
  titlePart[1]="e+";
  titlePart[2]="d";
  
  string part[3];
  part[0]="el";
  part[1]="po";
  part[2]="deut";

  for(int i=0; i<3;i++){

    string title = titlePart[i];
    string pName= part[i];

    double upLim=10;
    if(pName=="deut"){
      upLim=2;
    }

    TCanvas cP;
    TH1F *hP=new TH1F(("hP"+title).c_str(),(title+" Momentum").c_str(),100,0,upLim);
    hP->GetXaxis()->SetTitle("Momentum [GeV]");
    output->Draw((pName+"P>>hP"+title).c_str(),cut.c_str(),"colz");
    cP.Draw();
    cP.SaveAs((fileLoc+"Vars"+endName+".pdf").c_str());
  
    upLim=45;
    if(pName=="deut"){
      upLim=125;
    } else {
      TCanvas cRadCor;
      TH1F *hRadCor=new TH1F(("hRadCor"+title).c_str(),(title+" Momentum Correction (Radiated Photons)").c_str(),100,0,1.5);
      hRadCor->GetXaxis()->SetTitle("Momentum Correction [GeV]");
      output->Draw((pName+"RadCor>>hRadCor"+title).c_str(),cut.c_str(),"colz");
      cRadCor.SetLogy();
      cRadCor.Draw();
      cRadCor.SaveAs((fileLoc+"Vars"+endName+".pdf").c_str());

    }

    TCanvas ceTh;
    TH1F *heTh=new TH1F(("heTh"+title).c_str(),(title+" Theta").c_str(),100,0,upLim);
    heTh->GetXaxis()->SetTitle("Theta [Degrees]");
    output->Draw((pName+"Theta*(180./3.14)>>heTh"+title).c_str(),cut.c_str(),"colz");
    ceTh.Draw();
    ceTh.SaveAs((fileLoc+"Vars"+endName+".pdf").c_str());
    
    TCanvas cePh;
    TH1F *hePh=new TH1F(("hePh"+title).c_str(),(title+" Phi").c_str(),70,-180,180);
    hePh->GetXaxis()->SetTitle("Phi [Degrees]");
    output->Draw((pName+"Phi*(180./3.14)>>hePh"+title).c_str(),cut.c_str(),"colz");
    cePh.Draw();
    cePh.SaveAs((fileLoc+"Vars"+endName+".pdf").c_str());


    TCanvas ceSFETot;
    TH2F *heSFETot=new TH2F(("heSFETot"+title).c_str(),(title+" Sampling Fraction vs E_{Total}").c_str(),100,0.01,2,100,0.05,0.35);
    heSFETot->GetXaxis()->SetTitle("E_{Total}");
    heSFETot->GetYaxis()->SetTitle("Sampling Fraction");
    output->Draw((pName+"SF:"+pName+"EDep>>heSFETot"+title).c_str(),cut.c_str(),"colz");
    ceSFETot.Draw();
    ceSFETot.SaveAs((fileLoc+"Vars"+endName+".pdf").c_str());

  }

  TCanvas celPCALSFECinSF;
  TH2F *helPCALSFECinSF=new TH2F("helPCALSFECinSF","e^{-} E_{PCAL}/P vs E_{ECin}/P",100,0.01,0.2,100,0.01,0.3);
  helPCALSFECinSF->GetXaxis()->SetTitle("E_{ECin}/P");
  helPCALSFECinSF->GetYaxis()->SetTitle("E_{PCAL}/P");
  output->Draw("elSFPCAL:elSFECIN>>helPCALSFECinSF",cut.c_str(),"colz");
  celPCALSFECinSF.Draw();
  celPCALSFECinSF.SaveAs((fileLoc+"Vars"+endName+".pdf").c_str());

  TCanvas cpoPCALSFECinSF;
  TH2F *hpoPCALSFECinSF=new TH2F("hpoPCALSFECinSF","e^{+} E_{PCAL}/P vs E_{ECin}/P",100,0.01,0.2,100,0.01,0.3);
  hpoPCALSFECinSF->GetXaxis()->SetTitle("E_{ECin}/P");
  hpoPCALSFECinSF->GetYaxis()->SetTitle("E_{PCAL}/P");
  output->Draw("poSFPCAL:poSFECIN>>hpoPCALSFECinSF",cut.c_str(),"colz");
  cpoPCALSFECinSF.Draw();
  cpoPCALSFECinSF.SaveAs((fileLoc+"Vars"+endName+".pdf").c_str());

  TF1* deutMass = new TF1("deutMass","x/sqrt(x*x+1.875612*1.875612)");
  deutMass->SetLineColor(kBlack);
  deutMass->SetLineWidth(2);
  deutMass->SetLineStyle(9);
  deutMass->SetRange(0,5.0);

  TCanvas cdeutBetaP;
  TH2F *hdeutBetaP=new TH2F("hdeutBetaP","d #beta vs P",100,0.,5.0,100,0.,1.0);
  hdeutBetaP->GetXaxis()->SetTitle("P [GeV]");
  hdeutBetaP->GetYaxis()->SetTitle("#beta");
  output->Draw("deutBeta:deutP>>hdeutBetaP",cut.c_str(),"colz");
  deutMass->Draw("same");
  cdeutBetaP.Draw();
  cdeutBetaP.SaveAs((fileLoc+"Vars"+endName+".pdf)").c_str());




}
