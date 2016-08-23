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

  // TODO: add a 3rd argument to check if output file is specified and use that to name output file
  string output_filename = "out.root";
  
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
  vector<float> ph_pt;
  vector<float> ph_e;
  vector<float> ph_phi;
  vector<float> ph_eta;
  vector<float> ph_rapidity;
  vector<float> ph_m;

  vector<float> jet_pt;
  vector<float> jet_e;
  vector<float> jet_phi;
  vector<float> jet_eta;
  vector<float> jet_rapidity;
  vector<float> jet_m;


  // add branches to the tree
  output_tree->Branch("ph_pt", &ph_pt);
  output_tree->Branch("ph_e", &ph_e);
  output_tree->Branch("ph_phi", &ph_phi);
  output_tree->Branch("ph_eta", &ph_eta); 
  output_tree->Branch("ph_rapidity", &ph_rapidity);
  output_tree->Branch("ph_m", &ph_m);

  output_tree->Branch("jet_pt", &jet_pt);
  output_tree->Branch("jet_e", &jet_e);
  output_tree->Branch("jet_phi", &jet_phi);
  output_tree->Branch("jet_eta", &jet_eta);  
  output_tree->Branch("jet_rapidity", &jet_rapidity);
  output_tree->Branch("jet_m", &jet_m);
    
    
    
  for (int ievt=0; ievt<nevt; ++ievt){
    evt->getEntry(ievt);

    // important!! clear the branches on every loop
    ph_pt.clear();
    ph_e.clear();
    ph_phi.clear();
    ph_eta.clear();
    ph_rapidity.clear();
    ph_m.clear();

    jet_pt.clear();
    jet_e.clear();
    jet_phi.clear();
    jet_eta.clear();
    jet_rapidity.clear();
    jet_m.clear();

    
    //Declaration of variable (each event)
/*
    if (ievt % 500 ==0) {
      std::cout << "Processing event #" << ievt << "/" << nevt <<std::endl;
    }
*/ 
    const xAOD::TruthParticleContainer* photons(nullptr);
    evt->retrieve(photons, "Truth3Photons");
   
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
    cout<<"00000   PHOTON INFO   00000"<<endl;
    std::cout<<photons->size()<<"photons in this event: " <<std::endl;
  
    // Discard events that has no photons. Otherwise print out leading photon kinematics
    if (photons->size()<1) {
      std::cout<<"Event skipped: No photon"<<std::endl;  
       continue;
                            }
    else {

    //Finding the leading photon ID of the event  
      //xAOD::Photon* lead_photon = photons->at(0);
      int pho_idx = -1;
      TLorentzVector ph0;
      for (auto ph: *photons){ 
          pho_idx++;
        cout<<"Photon index: "<<pho_idx<<"  Photon Barcode: "<<ph->barcode()<<"  photon pt: "<<ph->pt()<<endl;
      //Find the leading photon index 
        if (ph->pt() > ph0.Pt()){
          ph0 = ph->p4();
                                    }    
                              } 

      // add the leading photon to the output tree
      ph_pt.push_back(ph0.Pt());
      ph_e.push_back(ph0.E());
      ph_phi.push_back(ph0.Phi());
      ph_eta.push_back(ph0.Eta());
      ph_rapidity.push_back(ph0.Rapidity());
      ph_m.push_back(ph0.M());

     const xAOD::JetContainer* jets(nullptr);
     evt->retrieve(jets, "TrimmedAntiKt10TruthJets");
      
      

      cout<<endl;
      cout<<"00000   JET INFO   00000"<<endl;
      cout<<jets->size()<<" jets in the event:"<<endl;
      cout<<"------------------------------"<<endl;
      

      //Reading the TrimmedAntiKt10TruthJets
     int jet_idx = -1;
     TLorentzVector j0;
     //j0.Pt()=-1;//initializing the value
     bool goodjet=false; //Check to see if we actually have jets outside of the exluded region
     for (auto j : *jets){
       jet_idx++;
       cout<<"Jet Index: "<<jet_idx<<"  Jet pt: "<<j->pt()<<endl;

       if ((ph0.DeltaR(j->p4()) > 1.0) && (j->pt()>j0.Pt())){
         j0= j->p4();
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

   //printing  leading photon info of the event
          std::cout<<"======Leading PHOTON Kinematics====="<<std::endl;
          std::cout<<"    Momentum and energy Kinematics"<<std::endl;
          cout<<"-----------------------------------------"<<endl;
          std::cout<<"The Total Energy of the leading photon is: " <<ph0.E()<<std::endl; 
          std::cout<<"The PT of the leading photon is: "<<ph0.Pt()<<std::endl;
          std::cout<<std::endl;
          std::cout<<"    Direction and Position Kinematics"<<std::endl;
          cout<<"-----------------------------------------"<<endl;
          std::cout<<"The phi of the leading photon is: " <<ph0.Phi()<<std::endl;
          std::cout<<"The pseudo-rapidity of the photon is: " <<ph0.Eta()<<std::endl;
          std::cout<<"The rapidity of the photon is: " <<ph0.Rapidity()<<std::endl;

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
