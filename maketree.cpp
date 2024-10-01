#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"

using namespace clas12;
using namespace std;


void initNames(string* varNames);

int getRow(string key,string * varNames,int nVars){
  for(int i=0;i<nVars;i++){
    if(varNames[i]==key){
      return i;
    }
  }
  cout<<"Variable name < "<<key<<" > not found..."<<endl;
  return -1;
}

int TriangCut(region_part_ptr part);

void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp,double M){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
	      rp->par()->getPz(),M);

}

Float_t GetMeanSF(Float_t  Edep, const TableOfDoubles_t& ccdbPhSF){
    return ccdbPhSF[0][3]*(ccdbPhSF[0][4]+(ccdbPhSF[0][5]/Edep)+(ccdbPhSF[0][6]/(Edep*Edep)));
}

void elRadCor(const TableOfDoubles_t& ccdbPhSF, clas12::region_part_ptr elrp,vector<clas12::region_part_ptr> nrps,TLorentzVector& el);

void getEs(region_part_ptr part, float * Es);

bool pass_Dead_Paddle_PCAL(region_part_ptr part);

void readConfig(string configFileName, string &RCDBPath, string &CCDBPath, string &treename,int &useGoldenRunsOnly);

//run with eg:
// clas12root -b 'maketree.cpp("/cache/clas12/rg-b/production/recon/spring2019/torus-1/pass2/v0/dst/train/jpsi/jpsi_006334.hipo","config.dat","eedFS_test.root")'
//to cat output to log file:
//clas12root -b 'maketree.cpp("/cache/clas12/rg-b/production/recon/spring2019/torus-1/pass2/v0/dst/train/jpsi/jpsi_006334.hipo","config.dat","eedFS_test.root")' > eedFS_test.log
//also works with wildcards in data path

void maketree(string dataPath, string configFileName, string outLoc){

    //Timer
    auto start = std::chrono::high_resolution_clock::now(); 

    //if no config file is provided, will use default
    if(dataPath==""){
      cout<<"\nNeed data file"<<endl;
      cout<<"Bye bye...\n"<<endl;
      exit(0);
    }

    //string configFileName="config.dat";
    string RCDBPath="",CCDBPath="",treeName="";
    int useGoldenRunsOnly=0;
    readConfig(configFileName, RCDBPath, CCDBPath, treeName,useGoldenRunsOnly);

    clas12databases::SetCCDBLocalConnection(CCDBPath.c_str());
    clas12databases::SetRCDBRootConnection(RCDBPath.c_str());

    cout<<"\nOutput tree file path: "<<outLoc<<"\n"<<endl;

    // Create a ROOT file to store the tree
    TFile outFile(outLoc.c_str(), "RECREATE");
    // Create a TTree and set its structure
    TTree* outTree = new TTree(treeName.c_str(), treeName.c_str());

    int nVars=67;
    string varNames[nVars];
    initNames(varNames);
    vector<Double_t> branchValues;
    TBranch* Branches[nVars];

    //branch points to element of vectors
    //memory location of vector changes when add to vector
    //so we add as many elements as needed then don't add more to vector
    for(int irow=0; irow<nVars;irow++){
        branchValues.push_back(-999);
    }
    //do not do this in single loop!
    //otherwise the branch won't have correct address!
    for(int irow=0; irow<nVars;irow++){
        Branches[irow]=outTree->Branch((varNames[irow]).c_str(), &branchValues[irow], (varNames[irow]+"/D").c_str());
    }

    //just used to get list of files
    TChain fake("hipo");
    //fake.Add("/cache/clas12/rg-b/production/recon/spring2019/torus-1/pass2/v0/dst/train/jpsi/*.hipo");
    fake.Add(dataPath.c_str());
    auto listfiles=fake.GetListOfFiles();

    //masses and lorentz vectors
    double massDeut=1.875612, massE=0.000510;
    TLorentzVector target(0, 0, 0, massDeut);
    TLorentzVector beam(0, 0, 10.6, 10.6);

    TLorentzVector el(0, 0, 0, massE);
    TLorentzVector po(0, 0, 0, massE);
    TLorentzVector deut(0, 0, 0, massDeut);

    gBenchmark->Start("timer");
    int counter = 0, entries=0;
    double accCharge=0;

    //iterate over list of files
    for (Int_t i = 0; i < listfiles->GetEntries(); i++)
    {
        cout<<"File "<<i+1<<"/"<<listfiles->GetEntries()<<endl;
        // create the event reader
        clas12reader c12(listfiles->At(i)->GetTitle());
        clas12databases c12db;
        c12.connectDataBases(&c12db);
        auto& rcdbData= c12.rcdb()->current();
        //only golden files
        if(useGoldenRunsOnly==1){
          c12.db()->qadb_requireGolden(true);
        }
        c12.applyQA();

        //beam_energy in MeV
        float beamE=(rcdbData.beam_energy)/1000.;
        //beam won't change on run per run basis
        beam.SetXYZM(0,0,beamE,massE);

        while (c12.next() == true)
        {

            double st=c12.event()->getStartTime();
            auto electrons = c12.getByID(11);
            auto deuterons = c12.getByID(45);
            auto positrons = c12.getByID(-11);
            auto neutrals = c12.getByCharge(0);
            
            if(electrons.size()>0 && positrons.size() > 0 && deuterons.size()>0){
              int combis=0;
              for(int el_combis=0;el_combis<electrons.size();el_combis++){
                  for(int po_combis=0;po_combis<positrons.size();po_combis++){
                    for(int deut_combis=0;deut_combis<deuterons.size();deut_combis++){
                      
                      for(int irow=0;irow<nVars;irow++){
                        branchValues[irow]=-999;
                      }

                      // set the particle momentum
                      SetLorentzVector(el, electrons[el_combis],massE); 
                      SetLorentzVector(po, positrons[po_combis],massE);  
                      SetLorentzVector(deut, deuterons[deut_combis],massDeut);        

                      double elPBefore=el.P();
                      double poPBefore=po.P();
                      elRadCor(c12.ccdb()->requestTableDoubles("/calibration/eb/photon_sf"), electrons[el_combis],neutrals,el);
                      double elRC=el.P()-elPBefore;
                      elRadCor(c12.ccdb()->requestTableDoubles("/calibration/eb/photon_sf"), positrons[po_combis],neutrals,po);
                      double poRC=po.P()-poPBefore;

                      TLorentzVector miss = beam + target - el - po - deut;
                      TLorentzVector inv = el + po;
                      double missP = miss.P();
                      double missPx = miss.Px();
                      double missPy = miss.Py();
                      double missPxP = missPx/missP;
                      double missPyP = missPy/missP;
                      double missPtP = sqrt(missPxP*missPxP + missPyP*missPyP);	 
                      double q2 = 2*beam.E()*missP*(1-cos(miss.Theta()));
                      double Ega =el.E() + po.E() + deut.E() - massDeut;
                      double t=-1*(deut-target).M2();

                      /*cout<<"\nBeamE "<<beamE<<" from vect "<<beam.E()<<" "<<beam.Pz()<<" "<<beam.M()<<endl;
                      cout<<"target "<<target.E()<<" "<<target.Pz()<<" "<<target.M()<<endl;
                      cout<<"el "<<el.E()<<" "<<el.Pz()<<" "<<el.M()<<endl;
                      cout<<"po "<<po.E()<<" "<<po.Pz()<<" "<<po.M()<<endl;
                      cout<<"deut "<<deut.E()<<" "<<deut.Pz()<<" "<<deut.M()<<endl;
                      cout<<"MM2 "<<miss.M2()<<" IM "<<inv.M()<<" Q2 "<<q2<<endl;*/

                      int elPd=0;
                      if(pass_Dead_Paddle_PCAL(electrons[el_combis])){elPd=1;}
                      int poPd=0;
                      if(pass_Dead_Paddle_PCAL(positrons[po_combis])){poPd=1;}

                      float elEs[8];
                      float poEs[8];
                      float deutEs[8];
                      getEs(electrons[el_combis],elEs);
                      getEs(positrons[po_combis],poEs);
                      getEs(deuterons[deut_combis],deutEs);

                      //use getRow function to avoid mixing up rows
                      branchValues[getRow("elSector", varNames,nVars)] = electrons[el_combis]->getSector();
                      branchValues[getRow("elP", varNames,nVars)] = el.P();
                      branchValues[getRow("elTheta", varNames,nVars)] = el.Theta();
                      branchValues[getRow("elPhi", varNames,nVars)] = el.Phi();
                      branchValues[getRow("elPCALLU", varNames,nVars)] = electrons[el_combis]->cal(PCAL)->getLu();
                      branchValues[getRow("elPCALLV", varNames,nVars)] = electrons[el_combis]->cal(PCAL)->getLv();
                      branchValues[getRow("elPCALLW", varNames,nVars)] = electrons[el_combis]->cal(PCAL)->getLw();
                      branchValues[getRow("elSF", varNames,nVars)] = elEs[4];
                      branchValues[getRow("elEDep", varNames,nVars)] = elEs[0];
                      branchValues[getRow("elSFPCAL", varNames,nVars)] = elEs[5];
                      branchValues[getRow("elSFECIN", varNames,nVars)] = elEs[6];
                      branchValues[getRow("elSFECOUT", varNames,nVars)] = elEs[7];
                      branchValues[getRow("elStatus", varNames,nVars)] = electrons[el_combis]->par()->getStatus();
                      branchValues[getRow("elPassDeadPaddlePCAL", varNames,nVars)] = elPd;
                      branchValues[getRow("elTriangCut", varNames,nVars)] = TriangCut(electrons[el_combis]);
                      branchValues[getRow("elVz", varNames,nVars)] = electrons[el_combis]->par()->getVz();
                      branchValues[getRow("elRadCor", varNames,nVars)] = elRC;

                      branchValues[getRow("poSector", varNames,nVars)] = positrons[po_combis]->getSector();
                      branchValues[getRow("poP", varNames,nVars)] = po.P();
                      branchValues[getRow("poTheta", varNames,nVars)] = po.Theta();
                      branchValues[getRow("poPhi", varNames,nVars)] = po.Phi();
                      branchValues[getRow("poPCALLU", varNames,nVars)] = positrons[po_combis]->cal(PCAL)->getLu();
                      branchValues[getRow("poPCALLV", varNames,nVars)] = positrons[po_combis]->cal(PCAL)->getLv();
                      branchValues[getRow("poPCALLW", varNames,nVars)] = positrons[po_combis]->cal(PCAL)->getLw();
                      branchValues[getRow("poSF", varNames,nVars)] = poEs[4];
                      branchValues[getRow("poEDep", varNames,nVars)] = poEs[0];
                      branchValues[getRow("poSFPCAL", varNames,nVars)] = poEs[5];
                      branchValues[getRow("poSFECIN", varNames,nVars)] = poEs[6];
                      branchValues[getRow("poSFECOUT", varNames,nVars)] = poEs[7];
                      branchValues[getRow("poStatus", varNames,nVars)] = positrons[po_combis]->par()->getStatus();
                      branchValues[getRow("poPassDeadPaddlePCAL", varNames,nVars)] = poPd;
                      branchValues[getRow("poTriangCut", varNames,nVars)] = TriangCut(positrons[po_combis]);
                      branchValues[getRow("poVz", varNames,nVars)] = positrons[po_combis]->par()->getVz();
                      branchValues[getRow("poRadCor", varNames,nVars)] = poRC;

                      branchValues[getRow("deutSector", varNames,nVars)] = deuterons[deut_combis]->getSector();
                      branchValues[getRow("deutP", varNames,nVars)] = deut.P();
                      branchValues[getRow("deutTheta", varNames,nVars)] = deut.Theta();
                      branchValues[getRow("deutPhi", varNames,nVars)] = deut.Phi();
                      branchValues[getRow("deutSF", varNames,nVars)] = deutEs[4];
                      branchValues[getRow("deutEDep", varNames,nVars)] = deutEs[0];
                      branchValues[getRow("deutDeltaEnergy", varNames,nVars)] = deuterons[deut_combis]->getDeltaEnergy();

                      if(abs(deuterons[deut_combis]->par()->getStatus()) >= 2000 && deuterons[deut_combis]->par()->getStatus() < 4000) {
                          branchValues[getRow("deutTOFEnergy", varNames,nVars)] = deuterons[deut_combis]->sci(FTOF1B)->getEnergy();
                          branchValues[getRow("deutTOFPath", varNames,nVars)] = deuterons[deut_combis]->sci(FTOF1B)->getPath();
                          branchValues[getRow("deutTOFTime", varNames,nVars)] = deuterons[deut_combis]->sci(FTOF1B)->getTime();
                      } else {
                          branchValues[getRow("deutTOFEnergy", varNames,nVars)] = deuterons[deut_combis]->sci(CTOF)->getEnergy();
                          branchValues[getRow("deutTOFPath", varNames,nVars)] = deuterons[deut_combis]->sci(CTOF)->getPath();
                          branchValues[getRow("deutTOFTime", varNames,nVars)] = deuterons[deut_combis]->sci(CTOF)->getTime();
                      }

                      branchValues[getRow("deutBeta", varNames,nVars)] = deuterons[deut_combis]->par()->getBeta();
                      branchValues[getRow("deutStatus", varNames,nVars)] = deuterons[deut_combis]->par()->getStatus();
                      branchValues[getRow("deutVz", varNames,nVars)] = deuterons[deut_combis]->par()->getVz();
                      branchValues[getRow("deutChi2PID", varNames,nVars)] = deuterons[deut_combis]->par()->getChi2Pid();

                      branchValues[getRow("elpoVzDif", varNames,nVars)] = positrons[po_combis]->par()->getVz() - electrons[el_combis]->par()->getVz();
                      branchValues[getRow("eldeutVzDif", varNames,nVars)] = deuterons[deut_combis]->par()->getVz() - electrons[el_combis]->par()->getVz();
                      branchValues[getRow("podeutVzDif", varNames,nVars)] = deuterons[deut_combis]->par()->getVz() - positrons[po_combis]->par()->getVz();
                      branchValues[getRow("elpoSectorDif", varNames,nVars)] = positrons[po_combis]->getSector() - electrons[el_combis]->getSector();
                      branchValues[getRow("eldeutSectorDif", varNames,nVars)] = deuterons[deut_combis]->getSector() - electrons[el_combis]->getSector();
                      branchValues[getRow("podeutSectorDif", varNames,nVars)] = deuterons[deut_combis]->getSector() - positrons[po_combis]->getSector();
                      branchValues[getRow("elpoHTCCTimeDif", varNames,nVars)] = positrons[po_combis]->che(HTCC)->getTime() - electrons[el_combis]->che(HTCC)->getTime();

                      branchValues[getRow("IM", varNames,nVars)] = inv.M();
                      branchValues[getRow("MM2", varNames,nVars)] = miss.M2();
                      branchValues[getRow("MM", varNames,nVars)] = miss.M();
                      branchValues[getRow("Q2", varNames,nVars)] = q2;
                      branchValues[getRow("Egamma", varNames,nVars)] = Ega;
                      branchValues[getRow("PtP", varNames,nVars)] = missPtP;
                      branchValues[getRow("t", varNames,nVars)] = t;
                      branchValues[getRow("BeamE", varNames,nVars)] = beamE;
                      branchValues[getRow("Combis", varNames,nVars)] = combis;
                      branchValues[getRow("RunNb", varNames,nVars)] = c12.getRunNumber();
                      branchValues[getRow("EventNb", varNames,nVars)] = c12.runconfig()->getEvent();
                      branchValues[getRow("UID", varNames,nVars)] = entries;

                      /*for(int irow=0;irow<nVars;irow++){
                        cout<<"branch "<<varNames[irow]<<" "<<branchValues[irow]<<endl;
                      }*/

                      //outFile.cd();
                      //write out data
                      outTree->Fill();
                      combis++;
                      entries++;
                  }//iterate over deuterons
                }//iterate over positrons
              }//iterate over electrons   
              counter++;
            }//at least 1 electron 1 positron 1 deuteron
            accCharge+=c12.db()->qa()->getAccCharge();
        }//while read file

   }//loop over files

    outFile.cd();
    outTree->Write();

    // Close the file
    outFile.Close();

   gBenchmark->Stop("timer");
   gBenchmark->Print("timer");
  
   auto finish = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> elapsed = finish - start;
   std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<<" with "<<entries<<" entries in total\n";
   cout<<"Accumulated charge past QA: "<<accCharge<<" nC"<<endl;

   //leave root command line interpreter
   exit(1);

}


void elRadCor(const TableOfDoubles_t& ccdbPhSF, clas12::region_part_ptr elrp,vector<clas12::region_part_ptr> nrps,TLorentzVector& el){
    
    for(auto nrp:nrps){
        Float_t nTheta=nrp->getTheta();
        Float_t nPhi=nrp->getPhi();
        Float_t diffTheta= (nTheta - elrp->getTheta())*TMath::RadToDeg(); //difference in Theta

        Float_t EdepPCAL = nrp->cal(clas12::PCAL)->getEnergy();
        Float_t EdepECIN = nrp->cal(clas12::ECIN)->getEnergy();
        Float_t EdepECOUT = nrp->cal(clas12::ECOUT)->getEnergy();

        Float_t dist=0;
        //check which calorimeter neutral has hit in
        if(EdepPCAL>0.01){
          Float_t elLU=elrp->cal(clas12::PCAL)->getLu();
          Float_t elLV=elrp->cal(clas12::PCAL)->getLv();
          Float_t elLW=elrp->cal(clas12::PCAL)->getLw();
          Float_t nrpLU=nrp->cal(clas12::PCAL)->getLu();
          Float_t nrpLV=nrp->cal(clas12::PCAL)->getLv();
          Float_t nrpLW=nrp->cal(clas12::PCAL)->getLw();
          dist=sqrt( (elLU-nrpLU)*(elLU-nrpLU) + (elLV-nrpLV)*(elLV-nrpLV) +(elLW-nrpLW)*(elLW-nrpLW) );
        } else if(EdepECIN>0.01){
          Float_t elLU=elrp->cal(clas12::ECIN)->getLu();
          Float_t elLV=elrp->cal(clas12::ECIN)->getLv();
          Float_t elLW=elrp->cal(clas12::ECIN)->getLw();
          Float_t nrpLU=nrp->cal(clas12::ECIN)->getLu();
          Float_t nrpLV=nrp->cal(clas12::ECIN)->getLv();
          Float_t nrpLW=nrp->cal(clas12::ECIN)->getLw();
          dist=sqrt( (elLU-nrpLU)*(elLU-nrpLU) + (elLV-nrpLV)*(elLV-nrpLV) +(elLW-nrpLW)*(elLW-nrpLW) );
        } else if(EdepECOUT>0.01){
          Float_t elLU=elrp->cal(clas12::ECOUT)->getLu();
          Float_t elLV=elrp->cal(clas12::ECOUT)->getLv();
          Float_t elLW=elrp->cal(clas12::ECOUT)->getLw();
          Float_t nrpLU=nrp->cal(clas12::ECOUT)->getLu();
          Float_t nrpLV=nrp->cal(clas12::ECOUT)->getLv();
          Float_t nrpLW=nrp->cal(clas12::ECOUT)->getLw();
          dist=sqrt( (elLU-nrpLU)*(elLU-nrpLU) + (elLV-nrpLV)*(elLV-nrpLV) +(elLW-nrpLW)*(elLW-nrpLW) );
        }

        //theta diff will be small for photons radiated by electrons
        //don't want split offs one cluster split into two (so req dist>30cm)
        if(abs(diffTheta)<0.7 && dist>30){
            
            if(nrp->getPid()!=22){
                //here we have a photon that was IDed as a neutron
                //need to recalculate its momentum based on sampling fraction parametrisation in CCDB
                //as is done for photons, then add this momentum back to the electron
                Float_t EdepTot = EdepPCAL+EdepECIN+EdepECOUT;
                Float_t reP = EdepTot/GetMeanSF(EdepTot,ccdbPhSF);
                TLorentzVector newNP4(reP*sin(nTheta)*cos(nPhi), reP*sin(nTheta)*sin(nPhi), reP*cos(nTheta), 0);
                el=el+newNP4;
            } else{
                //this is a photon IDed as a photon
                //momentum is already calculated correctly
                Float_t reP=nrp->getP();
                TLorentzVector newNP4(reP*sin(nTheta)*cos(nPhi), reP*sin(nTheta)*sin(nPhi), reP*cos(nTheta), 0);
                el=el+newNP4;
            }

        }

    }
}


bool pass_Dead_Paddle_PCAL(region_part_ptr part) {
    int sector=part->getSector();
    double U=part->cal(PCAL)->getLu();
    double V=part->cal(PCAL)->getLv();
    double W=part->cal(PCAL)->getLw();
    if (sector == 1) {
        return !((W > 72. && W < 93.) || (W > 210. && W < 231.));
    } else if (sector == 2) {
        return !((V > 100. && V < 115.));
    } else if (sector == 4) {
        return !((V > 228. && V < 242.) || (V < 15.));
    } else if (sector == 6) {
        return !((W > 170. && W < 194.));
    } else {
        return true;
    }
}

int TriangCut(region_part_ptr part){
    double ePCAL=part->cal(PCAL)->getEnergy();
    double eECIN=part->cal(ECIN)->getEnergy();
    double eECOUT=part->cal(ECOUT)->getEnergy();
    double SFPCAL=ePCAL/part->getP();
    double SFECIN=eECIN/part->getP();
    double SFECOUT=eECOUT/part->getP();
    int pass_triangle=0;
    if (part->getP() < 4.5)
    {
        pass_triangle = 1;
    }
    else
    {
        if(SFECIN > (0.2 - SFPCAL)){
            pass_triangle=1;
        }
    }
    return pass_triangle;
}

void getEs(region_part_ptr part, float * Es){
    double ePCAL=part->cal(PCAL)->getEnergy();
    double eECIN=part->cal(ECIN)->getEnergy();
    double eECOUT=part->cal(ECOUT)->getEnergy();
    double SFPCAL=ePCAL/part->getP();
    double SFECIN=eECIN/part->getP();
    double SFECOUT=eECOUT/part->getP();
    double edep=ePCAL+eECIN+eECOUT;
    double SF=edep/part->getP();

    Es[0]=edep;
    Es[1]=ePCAL;
    Es[2]=eECIN;
    Es[3]=eECOUT;
    Es[4]=SF;
    Es[5]=SFPCAL;
    Es[6]=SFECIN;
    Es[7]=SFECOUT;

}

void readConfig(string configFileName, string &RCDBPath, string &CCDBPath, string &treename,int &useGoldenRunsOnly){

    ifstream inpconfig(configFileName.c_str());
    
    string gap="\n";
    //load settings
    map<std::string, std::string> m_Settings;
    if( inpconfig.is_open() ){
        while( !inpconfig.eof() ){
            std::string Key;
            std::string Val;
            inpconfig>>Key;
            inpconfig>>Val;
            m_Settings[Key] = Val;
            //cout<<setw(10)<<Key<<setw(20)<<m_Settings[Key]<<endl;
        }
    } else{
        cout<<"\nCannot open config file: "<<configFileName<<endl;
        cout<<"Using default settings"<<endl;
        gap="";
    }

    //initialise default
    RCDBPath="/w/work/clas12/tyson/analysis_code_jlab/databases/rcdb.root";
    CCDBPath="/w/work/clas12/tyson/analysis_code_jlab/databases/ccdb.sqlite";
    treename="eed";
    useGoldenRunsOnly=0;
    
    for( map<std::string, std::string>::iterator it =  m_Settings.begin(); it!= m_Settings.end(); it++ ){
    
      std::string key = (*it).first;
      std::string val = (*it).second;

      if( key.compare("RCDBPath") == 0 ){
          RCDBPath = val;
      } else if( key.compare("CCDBPath") == 0 ){
          CCDBPath = val;
      } else if( key.compare("treeName") == 0 ){
          treename = val;
      } else if( key.compare("GoldenRuns") == 0 ){
          useGoldenRunsOnly = atoi(val.c_str());
      }
        
    }

    cout<<gap<<"RCDBPath: "<<RCDBPath<<endl;
    cout<<"CCDBPath: "<<CCDBPath<<endl;
    cout<<"ROOT tree name: "<<treename<<endl;
    cout<<"Golden Runs: "<<useGoldenRunsOnly<<"\n"<<endl;
}

void initNames(string* varNames){

    varNames[0]="elSector";
    varNames[1]="elP";
    varNames[2]="elTheta";
    varNames[3]="elPhi";
    varNames[4]="elPCALLU";
    varNames[5]="elPCALLV";
    varNames[6]="elPCALLW";
    varNames[7]="elSF";
    varNames[8]="elEDep";
    varNames[9]="elSFPCAL";
    varNames[10]="elSFECIN";
    varNames[11]="elSFECOUT";
    varNames[12]="elStatus";
    varNames[13]="elPassDeadPaddlePCAL";


    varNames[14]="poSector";
    varNames[15]="poP";
    varNames[16]="poTheta";
    varNames[17]="poPhi";
    varNames[18]="poPCALLU";
    varNames[19]="poPCALLV";
    varNames[20]="poPCALLW";
    varNames[21]="poSF";
    varNames[22]="poEDep";
    varNames[23]="poSFPCAL";
    varNames[24]="poSFECIN";
    varNames[25]="poSFECOUT";
    varNames[26]="poStatus";


    varNames[27]="poPassDeadPaddlePCAL";
    varNames[28]="elTriangCut";
    varNames[29]="poTriangCut";
    varNames[30]="elRadCor";
    varNames[31]="poRadCor";
    varNames[32]="elVz";
    varNames[33]="poVz";


    varNames[34]="deutSector";
    varNames[35]="deutP";
    varNames[36]="deutTheta";
    varNames[37]="deutPhi";
    varNames[38]="deutSF";
    varNames[39]="deutEDep";
    varNames[40]="deutDeltaEnergy";
    varNames[41]="deutTOFEnergy";
    varNames[42]="deutTOFPath";
    varNames[43]="deutTOFTime";
    varNames[44]="deutBeta";
    varNames[45]="deutStatus";
    varNames[46]="deutVz";
    varNames[47]="deutChi2PID";


    varNames[48]="MM2";
    varNames[49]="MM";
    varNames[50]="Q2";
    varNames[51]="PtP";
    varNames[52]="IM";
    varNames[53]="Egamma";
    varNames[54]="t";
    varNames[55]="BeamE";
    varNames[56]="Combis";
    varNames[57]="RunNb";
    varNames[58]="EventNb";
    varNames[59]="UID";

    varNames[60]="elpoVzDif";
    varNames[61]="eldeutVzDif";
    varNames[62]="podeutVzDif";
    varNames[63]="elpoSectorDif";
    varNames[64]="eldeutSectorDif";
    varNames[65]="podeutSectorDif";
    varNames[66]="elpoHTCCTimeDif";

}

