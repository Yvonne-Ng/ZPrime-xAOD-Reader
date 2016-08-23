#include <iostream>
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODJet/JetContainer.h"
#include "TFile.h"
#include "TTree.h"
#include "xAODEgamma/Photon.h"


//using jtruth_reader.cxx and tau2/tau1 to find the QCD and the resonance jets 

using namespace std;

int main(int argc, char** argv){
  if (argc<2){
    std::cerr<< "Please specify the file name!"<< std ::endl;
    return 1;
  }

  string input_filename= argv[1];
  string output_filename = "j2out.root";

//Edits for j2
  int swapCount=0;

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

  vector<float> Ejet0_pt;
  vector<float> Ejet0_e;
  vector<float> Ejet0_phi;
  vector<float> Ejet0_eta;
  vector<float> Ejet0_rapidity;
  vector<float> Ejet0_m;
  vector<float> Ejet0_tau1;
  vector<float> Ejet0_tau2;
  vector<float> Ejet0_tau21;

  vector<float> Ejet1_pt;
  vector<float> Ejet1_e;
  vector<float> Ejet1_phi;
  vector<float> Ejet1_eta;
  vector<float> Ejet1_rapidity;
  vector<float> Ejet1_m;
  vector<float> Ejet1_tau1;
  vector<float> Ejet1_tau2;
  vector<float> Ejet1_tau21;

// Ekt4 kt4 found after edit: tau21 comparison between 2 kt10 jets
  vector<float> Ekt40_pt;
  vector<float> Ekt40_e;
  vector<float> Ekt40_phi;
  vector<float> Ekt40_eta;
  vector<float> Ekt40_rapidity;
  vector<float> Ekt40_m;
  /*
  vector<float> Ekt40_tau1;
  vector<float> Ekt40_tau2;
  vector<float> Ekt40_tau21;
*/


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

  output_tree->Branch("Ejet0_pt", &Ejet0_pt);
  output_tree->Branch("Ejet0_e", &Ejet0_e);
  output_tree->Branch("Ejet0_phi", &Ejet0_phi);
  output_tree->Branch("Ejet0_eta", &Ejet0_eta);  
  output_tree->Branch("Ejet0_rapidity", &Ejet0_rapidity);
  output_tree->Branch("Ejet0_m", &Ejet0_m);
  output_tree->Branch("Ejet0_tau1", &Ejet0_tau1);
  output_tree->Branch("Ejet0_tau2", &Ejet0_tau2);
  output_tree->Branch("Ejet0_tau21", &Ejet0_tau21);

  output_tree->Branch("Ejet1_pt", &Ejet1_pt);
  output_tree->Branch("Ejet1_e", &Ejet1_e);
  output_tree->Branch("Ejet1_phi", &Ejet1_phi);
  output_tree->Branch("Ejet1_eta", &Ejet1_eta);  
  output_tree->Branch("Ejet1_rapidity", &Ejet1_rapidity);
  output_tree->Branch("Ejet1_m", &Ejet1_m);
  output_tree->Branch("Ejet1_tau1", &Ejet1_tau1);
  output_tree->Branch("Ejet1_tau2", &Ejet1_tau2);
  output_tree->Branch("Ejet1_tau21", &Ejet1_tau21);

  output_tree->Branch("Ekt40_pt", &Ekt40_pt);
  output_tree->Branch("Ekt40_e", &Ekt40_e);
  output_tree->Branch("Ekt40_phi", &Ekt40_phi);
  output_tree->Branch("Ekt40_eta", &Ekt40_eta); 
  output_tree->Branch("Ekt40_rapidity", &Ekt40_rapidity);
  output_tree->Branch("Ekt40_m", &Ekt40_m);
  /*
  output_tree->Branch("Ekt40_tau1", &Ekt40_tau1);
  output_tree->Branch("Ekt40_tau2", &Ekt40_tau2);
  output_tree->Branch("Ekt40_tau21", &Ekt40_tau21);
   */ 
    
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

    Ejet0_pt.clear();
    Ejet0_e.clear();
    Ejet0_phi.clear();
    Ejet0_eta.clear();
    Ejet0_rapidity.clear();
    Ejet0_m.clear();
    Ejet0_tau1.clear();
    Ejet0_tau2.clear();
    Ejet0_tau21.clear();


    Ejet1_pt.clear();
    Ejet1_e.clear();
    Ejet1_phi.clear();
    Ejet1_eta.clear();
    Ejet1_rapidity.clear();
    Ejet1_m.clear();
    Ejet1_tau1.clear();
    Ejet1_tau2.clear();
    Ejet1_tau21.clear();

    Ekt40_pt.clear();
    Ekt40_e.clear();
    Ekt40_phi.clear();
    Ekt40_eta.clear();
    Ekt40_rapidity.clear();
    Ekt40_m.clear();
    /*
    Ekt40_tau1.clear();
    Ekt40_tau2.clear();
    Ekt40_tau21.clear();
*/
    
    
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
        
    //TODO: 1. skip this event if less than one leading Kt4 jet
    //continue
    //TODO:  note the kinematics of the _leading leading Kt4 jet
    //hint: ph->p4() ref (leading Kt4 jets->size()<1) {
    //
    //         continue;
    //                }utnrs a TLorentzVector
    //TODO: read the "TrimmedAntiKt10TruthJets" contrainer from the XAOD
    //TODO: Loop through the jets and print out the pT and mass of the _leaing_et which is more than dealtaR >1.0 away formt he leading leading Kt4 jet
    //hint: TLorenetzVector has a method DeltaR()
    cout<<"00000   Kt4  INFO   00000"<<endl;
    std::cout<<kt4jets->size()<<"kt4 in this event: " <<std::endl;
  
    // Discard events that has no leading Kt4 jets. Otherwise print out leading leading Kt4 jet kinematics
    if (kt4jets->size()<1) {
      std::cout<<"Event skipped: kt4 jet"<<std::endl;  
       continue;
                            }
    else {

    //Finding the leading leading Kt4 jet ID of the event  
      //xAOD::Photon* lead_leading Kt4 jet = leading Kt4 jets->at(0);
      int kt4_idx = -1;
      TLorentzVector kt4_0;
      for (auto kt4: *kt4jets){ 
          kt4_idx++;
        cout<<"kt4 jet index: "<<kt4_idx<<"  kt4 pt: "<<kt4->pt()<<endl;
      //Find the leading leading Kt4 jet index 
        if (kt4->pt() > kt4_0.Pt()){
          kt4_0 = kt4->p4();
                                    }    
                              } 

      // add the leading leading Kt4 jet to the output tree
      kt4_pt.push_back(kt4_0.Pt());
      kt4_e.push_back(kt4_0.E());
      kt4_phi.push_back(kt4_0.Phi());
      kt4_eta.push_back(kt4_0.Eta());
      kt4_rapidity.push_back(kt4_0.Rapidity());
      kt4_m.push_back(kt4_0.M());

     const xAOD::JetContainer* jets(nullptr);
     evt->retrieve(jets, "TrimmedAntiKt10TruthJets");
      
      

//j2 tau2/tau1 edits 
    if (jets->size()<2){
      cout<<"less than two KT10 jets, event skipped"<<endl;
      continue;
                       }



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


//j2
//Continue if there are less than 2 kt10 jets in the event
  if (jets->size()<2){
    cout<<"Event skipped, there are less than 2 kt10 jets in the event"<<endl;
    continue;
  }
//j2

//edits for j2
    int Ejet_idx = -1;
    TLorentzVector Ej0, Ej1;
    
     const xAOD::Jet* Ej0_obj(nullptr);
     const xAOD::Jet* Ej1_obj(nullptr);
//j2: finding leading jet and the subleading jet kt10
     for (auto j : *jets){
       Ejet_idx++;
       if (j->pt()>Ej0.Pt()){
         Ej1=Ej0;
         Ej1_obj=Ej0_obj;
       
         Ej0=j->p4();
         Ej0_obj=j;
        
       }
       else if (j->pt()>Ej1.Pt()){
          Ej1=j->p4();
           Ej1_obj=j;
       }
       }

        //finding tau2 and tau 1
     float j0_tau1 = j0_obj->auxdata<float>("Tau1_wta");
     float j0_tau2 = j0_obj->auxdata<float>("Tau2_wta");
     float j0_tau21 = (j0_tau1>0) ? (j0_tau2/j0_tau1) : -1;

//edits for j2

     float Ej0_tau1 = Ej0_obj->auxdata<float>("Tau1_wta");
     float Ej0_tau2 = Ej0_obj->auxdata<float>("Tau2_wta");
     float Ej0_tau21 = (Ej0_tau1>0) ? (Ej0_tau2/Ej0_tau1): -1;


     float Ej1_tau1 = Ej1_obj->auxdata<float>("Tau1_wta");
     float Ej1_tau2 = Ej1_obj->auxdata<float>("Tau2_wta");
     float Ej1_tau21 = (Ej1_tau1>0) ? (Ej1_tau2/Ej1_tau1): -1;


//j2: Making j0 the variable the one with a smaller tau2/tau1
//j2: print out the leading and subleading jets after selection
    cout<<"j2: The leading kt10 jet PT: "<<Ej0.Pt()<<" tau1: "<<Ej0_tau1<<" tau2: "<<Ej0_tau2<<" tau21: "<<Ej0_tau21<<endl;
    cout<<"j2: The subleading kt10 jet PT: "<<Ej1.Pt()<<" tau1: "<<Ej1_tau1<<" tau2: "<<Ej1_tau2<<" tau21: "<<Ej1_tau21<<endl;
  
      


      //swap
      if (&Ej0_obj != &Ej1_obj){
        if (Ej0_tau21>Ej1_tau21){
          auto t_obj=Ej0_obj;
          auto t=Ej0;
          Ej0_obj=Ej1_obj;
          Ej0=Ej1;
          Ej1_obj=t_obj;
          Ej1=t;
          cout<<"swap happened!"<<endl;
          swapCount++;
        }
      }
      else {
        cout<<"Error! Ej0 and Ej1 shares the same address"<<endl;
      }


      //
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

     Ejet_idx = -1;
   //  TLorentzVector j0;
   //  const xAOD::Jet* j0_obj(nullptr);
     //j0.Pt()=-1;//initializing the value
     bool Egoodjet=false; //Check to see if we actually have jets outside of the exluded region
     //define a leading kt4 jet.
     TLorentzVector Ekt40;
     //   Ekt40.Pt()=0;

//     const xAOD::Jet* Ekt40_obj(nullptr);

     for (auto j : *kt4jets){
       Ejet_idx++;
       cout<<"kt4 Jet Index: "<<Ejet_idx<<" kt4 Jet pt: "<<j->pt()<<endl;

       if ((Ej0.DeltaR(j->p4()) > 1.0) && (j->pt()>Ekt40.Pt())){
         Ekt40= j->p4();
       //  Ekt40_obj = j;
         Egoodjet=true;
                                                          }
      }
     if (Egoodjet==false){
       cout<<"j2: Event skipped, there are no jet outside of the exlusion region"<<endl;
       continue;
                       }
/*
     float Ekt40tau1 = Ekt40_obj->auxdata<float>("Tau1_wta");
     float Ekt40tau2 = Ekt40_obj->auxdata<float>("Tau2_wta");
     float Ekt40tau21 = (Ekt40tau1>0) ? (Ekt40tau2/Ekt40tau1): -1;
*/


   //writing the jet information into the vector
     jet_pt.push_back(j0.Pt());
     jet_e.push_back(j0.E());
     jet_phi.push_back(j0.Phi());
     jet_eta.push_back(j0.Eta());
     jet_rapidity.push_back(j0.Rapidity());
     jet_m.push_back(j0.M());

     jet_tau1.push_back(j0_tau1);
     jet_tau2.push_back(j0_tau2);
     jet_tau21.push_back(j0_tau21);
//edits for j2
//

     Ejet0_pt.push_back(Ej0.Pt());
     Ejet0_e.push_back(Ej0.E());
     Ejet0_phi.push_back(Ej0.Phi());
     Ejet0_eta.push_back(Ej0.Eta());
     Ejet0_rapidity.push_back(Ej0.Rapidity());
     Ejet0_m.push_back(Ej0.M());

     Ejet0_tau1.push_back(Ej0_tau1);
     Ejet0_tau2.push_back(Ej0_tau2);
     Ejet0_tau21.push_back(Ej0_tau21);


     Ejet1_pt.push_back(Ej1.Pt());
     Ejet1_e.push_back(Ej1.E());
     Ejet1_phi.push_back(Ej1.Phi());
     Ejet1_eta.push_back(Ej1.Eta());
     Ejet1_rapidity.push_back(Ej1.Rapidity());
     Ejet1_m.push_back(Ej1.M());

     Ejet1_tau1.push_back(Ej1_tau1);
     Ejet1_tau2.push_back(Ej1_tau2);
     Ejet1_tau21.push_back(Ej1_tau21);
     

    Ekt40_pt.push_back(Ekt40.Pt());
    Ekt40_e.push_back(Ekt40.E());
    Ekt40_phi.push_back(Ekt40.Phi());
    Ekt40_eta.push_back(Ekt40.Eta());
    Ekt40_rapidity.push_back(Ekt40.Rapidity());
    Ekt40_m.push_back(Ekt40.M());
/*
    Ekt40_tau1.push_back(Ekt40tau1);
    Ekt40_tau2.push_back(Ekt40tau2);
    Ekt40_tau21.push_back(Ekt40tau21);
*/

//edits for j2

     
   //printing  leading leading Kt4 jet info of the event
          std::cout<<"======Leading kt4 jet  Kinematics====="<<std::endl;
          std::cout<<"    Momentum and energy Kinematics"<<std::endl;
          cout<<"-----------------------------------------"<<endl;
          std::cout<<"The Total Energy of the leading Kt4 jet is: " <<kt4_0.E()<<std::endl; 
          std::cout<<"The PT of the leading Kt4 jet is: "<<kt4_0.Pt()<<std::endl;
          std::cout<<std::endl;
          std::cout<<"    Direction and Position Kinematics"<<std::endl;
          cout<<"-----------------------------------------"<<endl;
          std::cout<<"The phi of the leading Kt4 jet is: " <<kt4_0.Phi()<<std::endl;
          std::cout<<"The pseudo-rapidity of the leading Kt4 jet is: " <<kt4_0.Eta()<<std::endl;
          std::cout<<"The rapidity of the leading Kt4 jet is: " <<kt4_0.Rapidity()<<std::endl;

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
  
 // cout<<"% of j2 swap:"<< swapCount/nevt<<endl; 
  return 0;
}
