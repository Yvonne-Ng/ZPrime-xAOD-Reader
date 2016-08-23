#include <iostream>
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODJet/JetContainer.h"
#include "TFile.h"
#include "TTree.h"
#include "xAODEgamma/Photon.h"

using namespace std;

int main(int argc, char** argv){
  if (argc<2){
    std::cerr<< "Please specify the file name!"<< std ::endl;
    return 1;
  }

  string input_filename= argv[1];
  string output_filename = "jout.root";

  if (argc>2) {
    output_filename = argv[2];
  }

  std::cout << "Will write output to " << output_filename << "!" << std::endl;


  // TODO: add a 3rd argument to check if output file is specified and use that to name output file
  
  std::cout <<"Reading from file: " << input_filename<<std::endl;
  xAOD::Init();

  TFile *input_file = TFile::Open(input_filename.c_str());

  xAOD::TEvent* evt= new xAOD::TEvent(input_file);

  int nevt = evt->getEntries();
  std::cout<<"There are "<<nevt <<" event in the input file." << std::endl;

  std::cout<<"========================================================="<<std::endl;

  // important: create output file after loading the input file
  TFile* output_file = new TFile(output_filename.c_str(), "recreate");

  // create output tree after creating output file
  TTree* output_tree = new TTree("truth", "truth");

  // make vectors to hold variables
  vector<float> kt4_pt;
  vector<float> kt4_e;
  vector<float> kt4_phi;
  vector<float> kt4_eta;
  vector<float> kt4_rapidity;
  vector<float> kt4_m;

  vector<float> jet_pt;
  vector<float> jet_e;
  vector<float> jet_phi;
  vector<float> jet_eta;
  vector<float> jet_rapidity;
  vector<float> jet_m;
  vector<float> jet_tau1;
  vector<float> jet_tau2;
  vector<float> jet_tau21;


  // add branches to the tree
  output_tree->Branch("kt4_pt", &kt4_pt);
  output_tree->Branch("kt4_e", &kt4_e);
  output_tree->Branch("kt4_phi", &kt4_phi);
  output_tree->Branch("kt4_eta", &kt4_eta); 
  output_tree->Branch("kt4_rapidity", &kt4_rapidity);
  output_tree->Branch("kt4_m", &kt4_m);

  output_tree->Branch("jet_pt", &jet_pt);
  output_tree->Branch("jet_e", &jet_e);
  output_tree->Branch("jet_phi", &jet_phi);
  output_tree->Branch("jet_eta", &jet_eta);  
  output_tree->Branch("jet_rapidity", &jet_rapidity);
  output_tree->Branch("jet_m", &jet_m);
  output_tree->Branch("jet_tau1", &jet_tau1);
  output_tree->Branch("jet_tau2", &jet_tau2);
  output_tree->Branch("jet_tau21", &jet_tau21);
    
    
    
  for (int ievt=0; ievt<nevt; ++ievt){
    evt->getEntry(ievt);

    // important!! clear the branches on every loop
    kt4_pt.clear();
    kt4_e.clear();
    kt4_phi.clear();
    kt4_eta.clear();
    kt4_rapidity.clear();
    kt4_m.clear();

    jet_pt.clear();
    jet_e.clear();
    jet_phi.clear();
    jet_eta.clear();
    jet_rapidity.clear();
    jet_m.clear();
    jet_tau1.clear();
    jet_tau2.clear();
    jet_tau21.clear();

    
    //Declaration of variable (each event)
/*
    if (ievt % 500 ==0) {
      std::cout << "Processing event #" << ievt << "/" << nevt <<std::endl;
    }
*/ 
    const xAOD::JetContainer* kt4jets(nullptr);
    evt->retrieve(kt4jets, "AntiKt4TruthJets");
   
    std::cout<<"==================Event ID: "<<ievt<<"=================="<<std::endl;
    std::cout<<std::endl;
        
    //TODO: 1. skip this event if less than one photon
    //continue
    //TODO:  note the kinematics of the _leading photon
    //hint: ph->p4() ref (photons->size()<1) {
    //
    //         continue;
    //                }utnrs a TLorentzVector
    //TODO: read the "TrimmedAntiKt10TruthJets" contrainer from the XAOD
    //TODO: Loop through the jets and print out the pT and mass of the _leaing_et which is more than dealtaR >1.0 away formt he leading photon
    //hint: TLorenetzVector has a method DeltaR()
    cout<<"00000   Kt4  INFO   00000"<<endl;
    std::cout<<kt4jets->size()<<"kt4 in this event: " <<std::endl;
  
    // Discard events that has no photons. Otherwise print out leading photon kinematics
    if (kt4jets->size()<1) {
      std::cout<<"Event skipped: kt4 jet"<<std::endl;  
       continue;
                            }
    else {

    //Finding the leading photon ID of the event  
      //xAOD::Photon* lead_photon = photons->at(0);
      int kt4_idx = -1;
      TLorentzVector kt4_0;
      for (auto kt4: *kt4jets){ 
          kt4_idx++;
        cout<<"kt4 jet index: "<<kt4_idx<<"  kt4 pt: "<<kt4->pt()<<endl;
      //Find the leading photon index 
        if (kt4->pt() > kt4_0.Pt()){
          kt4_0 = kt4->p4();
                                    }    
                              } 

      // add the leading photon to the output tree
      kt4_pt.push_back(kt4_0.Pt());
      kt4_e.push_back(kt4_0.E());
      kt4_phi.push_back(kt4_0.Phi());
      kt4_eta.push_back(kt4_0.Eta());
      kt4_rapidity.push_back(kt4_0.Rapidity());
      kt4_m.push_back(kt4_0.M());

     const xAOD::JetContainer* jets(nullptr);
     evt->retrieve(jets, "TrimmedAntiKt10TruthJets");
      
      

      cout<<endl;
      cout<<"00000   JET INFO   00000"<<endl;
      cout<<jets->size()<<" jets in the event:"<<endl;
      cout<<"------------------------------"<<endl;
      

      //Reading the TrimmedAntiKt10TruthJets
     int jet_idx = -1;
     TLorentzVector j0;
     const xAOD::Jet* j0_obj(nullptr);
     //j0.Pt()=-1;//initializing the value
     bool goodjet=false; //Check to see if we actually have jets outside of the exluded region
     for (auto j : *jets){
       jet_idx++;
       cout<<"Jet Index: "<<jet_idx<<"  Jet pt: "<<j->pt()<<endl;

       if ((kt4_0.DeltaR(j->p4()) > 1.0) && (j->pt()>j0.Pt())){
         j0= j->p4();
         j0_obj = j;
         goodjet=true;
                                                          }
                          }
     if (goodjet==false){
       cout<<"Event skipped, there are no jet outside of the exlusion region"<<endl;
       continue;
                       }
/*
      TLorentzVector jet0;
      for (j : jets) {
        if (ph0.DeltaR(j->p4()) > 1.0
            && j->pt() > jet0.Pt()) {
          jet0 = j->p4();
        }
        }
      }
*/
   //writing the jet information into the vector
     jet_pt.push_back(j0.Pt());
     jet_e.push_back(j0.E());
     jet_phi.push_back(j0.Phi());
     jet_eta.push_back(j0.Eta());
     jet_rapidity.push_back(j0.Rapidity());
     jet_m.push_back(j0.M());

     float tau1 = j0_obj->auxdata<float>("Tau1_wta");
     float tau2 = j0_obj->auxdata<float>("Tau2_wta");
     float tau21 = (tau1>0) ? (tau2/tau1) : -1;

     jet_tau1.push_back(tau1);
     jet_tau2.push_back(tau2);
     jet_tau21.push_back(tau21);


   //printing  leading photon info of the event
          std::cout<<"======Leading kt4 jet  Kinematics====="<<std::endl;
          std::cout<<"    Momentum and energy Kinematics"<<std::endl;
          cout<<"-----------------------------------------"<<endl;
          std::cout<<"The Total Energy of the leading photon is: " <<kt4_0.E()<<std::endl; 
          std::cout<<"The PT of the leading photon is: "<<kt4_0.Pt()<<std::endl;
          std::cout<<std::endl;
          std::cout<<"    Direction and Position Kinematics"<<std::endl;
          cout<<"-----------------------------------------"<<endl;
          std::cout<<"The phi of the leading photon is: " <<kt4_0.Phi()<<std::endl;
          std::cout<<"The pseudo-rapidity of the photon is: " <<kt4_0.Eta()<<std::endl;
          std::cout<<"The rapidity of the photon is: " <<kt4_0.Rapidity()<<std::endl;

    // Printing out leading jet info of the event
          cout<<endl;
          cout<<"=====Leading JET Kinematics====="<<endl;
          cout<<"    Momentum and energy kinenmatics"<<endl;
          cout<<"-----------------------------------------"<<endl;
          cout<<"The Total Energy of the leading Jet: "<<j0.E()<<endl;
          cout<<"The PT of the leading Jet: "<<j0.Pt()<<endl;
          cout<<endl;
          cout<<"    Direction and position kinematics"<<endl;
          cout<<"-----------------------------------------"<<endl;
          cout<<"The phi of the leading jet"<<j0.Phi()<<endl;
          cout<<"The pseudo-rapidity of the leading jet"<< j0.Eta()<<std::endl;
          cout<<"The rapidity of the leading jet"<<j0.Rapidity()<<endl;
          cout<<endl;
    

    // IMPORTANT! save the latest entries in the tree
    output_tree->Fill();
          }
  }   

  output_file->Write();
  output_file->Close();
  
  return 0;
}
