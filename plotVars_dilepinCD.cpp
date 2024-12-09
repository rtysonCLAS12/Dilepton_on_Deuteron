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

void plotVars_dilepinCD(){

  gStyle->SetOptStat("");
  gStyle->SetOptFit(0012);

  string fileLoc = "/w/work/clas12/tyson/plots/c12_scripts/eed/IM/";

  string treeLocRoot = "/w/work/clas12/tyson/data_repo/c12scripts_out/eed/eedFS_allRGB_dilepinCD.root";

  string endNameBase="_allRGB"; //_noCuts

  auto file = new TFile((treeLocRoot).c_str());
  //auto output =(TTree*) file->Get("FINALOUTTREE");
  auto output =(TTree*) file->Get("eed");

  gStyle->SetOptStat("");

  string cutBase="";

  cutBase+="Q2<0.5 && abs(MM2)<0.2 && abs(elStatus)>=2000 && abs(poStatus)>=2000";
  cutBase+="&& deutChi2PID<5";

  string cutFD="&& abs(deutStatus)>=2000 && abs(deutStatus)<4000";
  string cutCD="&& abs(deutStatus)>=4000";

  string cut_tighter=cutBase;
  cut_tighter+="&& elTriangCut==1 && poTriangCut==1"; //e+ e- ID
  cut_tighter+="&& Combis==0";
  
  string cutIMBase=cut_tighter;

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

  string endNameC[5];
  endNameC[0]="_dilepInBothCDAndFD";
  endNameC[1]="_dilepOnlyInCD";
  endNameC[2]="_dilepOnlyInFD";
  endNameC[3]="_emFD_epCD";
  endNameC[4]="_emCD_epFD";

  string statusReq[5];
  statusReq[0]="";
  statusReq[1]="&& abs(elStatus)>=4000 && abs(poStatus)>=4000 "; // in CD
  statusReq[2]="&& abs(elStatus)<4000 && abs(poStatus)<4000"; //in FD
  statusReq[3]="&& abs(elStatus)<4000 && abs(poStatus)>=4000 "; //e- in FD, e+ in CD
  statusReq[4]="&& abs(elStatus)>=4000 && abs(poStatus)<4000"; // e- in CD, e+ in FD

  TF1* deutMass = new TF1("deutMass","x/sqrt(x*x+1.875612*1.875612)");
  deutMass->SetLineColor(kBlack);
  deutMass->SetLineWidth(2);
  deutMass->SetLineStyle(9);
  deutMass->SetRange(0,5.0);

  for(int cutnb=0;cutnb<5;cutnb++){ 
    string endName=endNameBase+endNameC[cutnb];
    string cut=cutBase+statusReq[cutnb];
    string cutIM=cutIMBase+statusReq[cutnb];

    TCanvas cIMFit;
    TH1F *hIMFit=new TH1F(("hIMFit"+endName).c_str(),"e+ e- Invariant Mass",30,2.5,3.5);
    hIMFit->SetTitle("e+ e- Invariant Mass ");
    hIMFit->GetXaxis()->SetTitle("Invariant Mass [GeV]");
    hIMFit->SetLineWidth(2);
    output->Draw(("IM>>hIMFit"+endName).c_str(), cutIM.c_str(),"");
    hIMFit->Fit(f2, "", "", 2.5, 3.5);
    hIMFit->SetFillColor(kAzure-9);
    gauss2->SetParameters(f2->GetParameter(0),f2->GetParameter(1),f2->GetParameter(2));
    bg2->SetParameters(f2->GetParameter(1),f2->GetParameter(3),f2->GetParameter(4),f2->GetParameter(5),f2->GetParameter(6));
    hIMFit->Draw();
    gauss2->Draw("same");
    bg2->Draw("same");
    cIMFit.Draw();
    cIMFit.SaveAs((fileLoc+"Vars"+endName+".pdf(").c_str());
    delete hIMFit;

    TCanvas cEgamma;
    TH1F *hEgamma=new TH1F(("hEgamma"+endName).c_str(),"E_{#gamma}",60,5,11);
    hEgamma->GetXaxis()->SetTitle("E_{#gamma} [GeV]");
    output->Draw(("Egamma>>hEgamma"+endName).c_str(), (cut).c_str());
    cEgamma.Draw();
    cEgamma.SaveAs((fileLoc+"Vars"+endName+".pdf").c_str());
    delete hEgamma;

    TCanvas cMM2;
    TH1F *hMM2=new TH1F(("hMM2"+endName).c_str(),"Missing Mass Squared",100,-1,1);
    hMM2->GetXaxis()->SetTitle("Missing Mass Squared [GeV^{2}]");
    output->Draw(("MM2>>hMM2"+endName).c_str(), (cut).c_str());
    cMM2.Draw();
    cMM2.SaveAs((fileLoc+"Vars"+endName+".pdf").c_str());
    delete hMM2;

    TCanvas cQ2;
    TH1F *hQ2=new TH1F(("hQ2"+endName).c_str(),"Q^{2}",100,0,0.5);
    hQ2->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    output->Draw(("Q2>>hQ2"+endName).c_str(),(cut).c_str());
    cQ2.Draw();
    cQ2.SaveAs((fileLoc+"Vars"+endName+".pdf").c_str());
    delete hQ2;

    TCanvas cIM;
    TH1F *hIM=new TH1F(("hIM"+endName).c_str(),"e^{+} e^{-} Invariant Mass",100,0,3.5);
    hIM->GetXaxis()->SetTitle("Invariant Mass [GeV]");
    output->Draw(("IM>>hIM"+endName).c_str(), (cutIM+"").c_str());
    cIM.Draw();
    cIM.SaveAs((fileLoc+"Vars"+endName+".pdf").c_str());
    delete hIM;

    TCanvas cCombis;
    TH1F *hCombis=new TH1F(("hCombis"+endName).c_str(),"Number of Permutations per Event",11,-0.5,10.5);
    hCombis->GetXaxis()->SetTitle("Number of Permutations per Event");
    output->Draw(("Combis>>hCombis"+endName).c_str(), (cut+"").c_str());
    cCombis.SetLogy();
    cCombis.Draw();
    cCombis.SaveAs((fileLoc+"Vars"+endName+".pdf").c_str());
    delete hCombis;

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
      TH1F *hP=new TH1F(("hP"+title+endName).c_str(),(title+" Momentum").c_str(),100,0,upLim);
      hP->GetXaxis()->SetTitle("Momentum [GeV]");
      output->Draw((pName+"P>>hP"+title+endName).c_str(),cut.c_str(),"colz");
      cP.Draw();
      cP.SaveAs((fileLoc+"Vars"+endName+".pdf").c_str());
      delete hP;
    
      upLim=45;
      if(pName=="deut"){
        upLim=125;
      } else {
        TCanvas cRadCor;
        TH1F *hRadCor=new TH1F(("hRadCor"+title+endName).c_str(),(title+" Momentum Correction (Radiated Photons)").c_str(),100,0,1.5);
        hRadCor->GetXaxis()->SetTitle("Momentum Correction [GeV]");
        output->Draw((pName+"RadCor>>hRadCor"+title+endName).c_str(),cut.c_str(),"colz");
        cRadCor.SetLogy();
        cRadCor.Draw();
        cRadCor.SaveAs((fileLoc+"Vars"+endName+".pdf").c_str());
        delete hRadCor;

      }

      TCanvas ceTh;
      TH1F *heTh=new TH1F(("heTh"+title+endName).c_str(),(title+" Theta").c_str(),100,0,upLim);
      heTh->GetXaxis()->SetTitle("Theta [Degrees]");
      output->Draw((pName+"Theta*(180./3.14)>>heTh"+title+endName).c_str(),cut.c_str(),"colz");
      ceTh.Draw();
      ceTh.SaveAs((fileLoc+"Vars"+endName+".pdf").c_str());
      delete heTh;
      
      TCanvas cePh;
      TH1F *hePh=new TH1F(("hePh"+title+endName).c_str(),(title+" Phi").c_str(),70,-180,180);
      hePh->GetXaxis()->SetTitle("Phi [Degrees]");
      output->Draw((pName+"Phi*(180./3.14)>>hePh"+title+endName).c_str(),cut.c_str(),"colz");
      cePh.Draw();
      cePh.SaveAs((fileLoc+"Vars"+endName+".pdf").c_str());
      delete hePh;


      TCanvas ceSFETot;
      TH2F *heSFETot=new TH2F(("heSFETot"+title+endName).c_str(),(title+" Sampling Fraction vs E_{Total}").c_str(),100,0.01,2,100,0.05,0.35);
      heSFETot->GetXaxis()->SetTitle("E_{Total}");
      heSFETot->GetYaxis()->SetTitle("Sampling Fraction");
      output->Draw((pName+"SF:"+pName+"EDep>>heSFETot"+title+endName).c_str(),cut.c_str(),"colz");
      ceSFETot.Draw();
      ceSFETot.SaveAs((fileLoc+"Vars"+endName+".pdf").c_str());
      delete heSFETot;

    }

    TCanvas celPCALSFECinSF;
    TH2F *helPCALSFECinSF=new TH2F(("helPCALSFECinSF"+endName).c_str(),"e^{-} E_{PCAL}/P vs E_{ECin}/P",100,0.01,0.2,100,0.01,0.3);
    helPCALSFECinSF->GetXaxis()->SetTitle("E_{ECin}/P");
    helPCALSFECinSF->GetYaxis()->SetTitle("E_{PCAL}/P");
    output->Draw(("elSFPCAL:elSFECIN>>helPCALSFECinSF"+endName).c_str(),cut.c_str(),"colz");
    celPCALSFECinSF.Draw();
    celPCALSFECinSF.SaveAs((fileLoc+"Vars"+endName+".pdf").c_str());
    delete helPCALSFECinSF;

    TCanvas cpoPCALSFECinSF;
    TH2F *hpoPCALSFECinSF=new TH2F(("hpoPCALSFECinSF"+endName).c_str(),"e^{+} E_{PCAL}/P vs E_{ECin}/P",100,0.01,0.2,100,0.01,0.3);
    hpoPCALSFECinSF->GetXaxis()->SetTitle("E_{ECin}/P");
    hpoPCALSFECinSF->GetYaxis()->SetTitle("E_{PCAL}/P");
    output->Draw(("poSFPCAL:poSFECIN>>hpoPCALSFECinSF"+endName).c_str(),cut.c_str(),"colz");
    cpoPCALSFECinSF.Draw();
    cpoPCALSFECinSF.SaveAs((fileLoc+"Vars"+endName+".pdf").c_str());
    delete hpoPCALSFECinSF;

    TCanvas cdeutBetaP;
    TH2F *hdeutBetaP=new TH2F(("hdeutBetaP"+endName).c_str(),"d #beta vs P",100,0.,5.0,100,0.,1.0);
    hdeutBetaP->GetXaxis()->SetTitle("P [GeV]");
    hdeutBetaP->GetYaxis()->SetTitle("#beta");
    output->Draw(("deutBeta:deutP>>hdeutBetaP"+endName).c_str(),cut.c_str(),"colz");
    deutMass->Draw("same");
    cdeutBetaP.Draw();
    cdeutBetaP.SaveAs((fileLoc+"Vars"+endName+".pdf)").c_str());
    delete hdeutBetaP;
  }




}
