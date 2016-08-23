#include <iostream>
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODJet/JetContainer.h"
#include "TFile.h"
#include "TTree.h"
#include "xAODEgamma/Photon.h"
#include "TLorentzVector.h"
#include <math.h>
#include <cmath>

using namespace std;
//Reading an input file and a treefile with a branch where it can be determined whether the event contain at least one particle of >330 GeV

using namespace std;

int main(int argc, char** argv){
  if (argc<2){
    std::cerr<< "Please specify the file name!"<< std ::endl;
    return 1;
  }


  string input_filename= argv[1];

  // TODO: add a 3rd argument to check if output file is specified and use that to name output file
  string output_filename = "out_mR90.root";
  
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
 // vector<float> cut330;
  vector<float> cut330;
  vector<TLorentzVector> jet;
  vector<float> jet_pt;
  vector<float> jet_pt_below330;
  vector<float> jet_pt_above330;
  vector<float> minJetPt;

  vector<float> alljets_pt;
  vector<float> alljets_pt_below330;
  vector<float> alljets_pt_above330;

  vector<TLorentzVector> kt4;
  vector<float> kt4_pt;
  vector<float> kt4_pt_below330;
  vector<float> kt4_pt_above330;
  vector<float> minKt4Pt;

  vector<float> allkt4_pt;
  vector<float> allkt4_pt_below330;
  vector<float> allkt4_pt_above330;


  // add branches to the tree
  output_tree->Branch("cut330", &cut330);    
  output_tree->Branch("jet_pt", &jet_pt);
  output_tree->Branch("jet_pt_below330", &jet_pt_below330);
  output_tree->Branch("jet_pt_above330", &jet_pt_above330);
  output_tree->Branch("minJetPt", &minJetPt);

  output_tree->Branch("alljets_pt", &alljets_pt);
  output_tree->Branch("alljets_pt_below330", &alljets_pt_below330);
  output_tree->Branch("alljets_pt_above330", &alljets_pt_above330);


  output_tree->Branch("kt4_pt", &kt4_pt);
  output_tree->Branch("kt4_pt_below330", &kt4_pt_below330);
  output_tree->Branch("kt4_pt_above330", &kt4_pt_above330);
  output_tree->Branch("minKt4Pt", &minKt4Pt);

  output_tree->Branch("allkt4_pt", &allkt4_pt);
  output_tree->Branch("allkt4_pt_below330", &allkt4_pt_below330);
  output_tree->Branch("allkt4_pt_above330", &allkt4_pt_above330);

    
  for (int ievt=0; ievt<nevt; ++ievt){
    evt->getEntry(ievt);

    // important!! clear the branches on every loop
    cut330.clear();
    jet.clear();
    jet_pt.clear();
    jet_pt_below330.clear();
    jet_pt_above330.clear();
    minJetPt.clear();

    alljets_pt.clear();
    alljets_pt_below330.clear();
    alljets_pt_above330.clear();


    kt4.clear();
    kt4_pt.clear();
    kt4_pt_below330.clear();
    kt4_pt_above330.clear();
    minKt4Pt.clear();

    allkt4_pt.clear();
    allkt4_pt_below330.clear();
    allkt4_pt_above330.clear();


    
    //Declaration of variable (each event)
/*
    if (ievt % 500 ==0) {
      std::cout << "Processing event #" << ievt << "/" << nevt <<std::endl;
    }
*/ 
 //   const xAOD::TruthParticleContainer* photons(nullptr);
   //Grabing the pdgID of the particle  
   //Grabing all particles 
    const xAOD::TruthParticleContainer* particles(nullptr);
    evt->retrieve(particles, "TruthParticles");

    const xAOD::JetContainer *kt4s(nullptr);
    evt->retrieve(kt4s, "AntiKt4TruthJets");
  
        
    float leadPt=0;
    float leadPtAbv330=0;
    float leadPtBel330=0;

    float leadKt4Pt=0;
    float subLeadKt4Pt=0;
    float subSubLeadKt4Pt=0;
    float leadKt4PtAbv330=0;
    float leadKt4PtBel330=0;


    int partoncount=0;
    for (auto par: *particles){
      
   //   cout<<"pdgID: "<<par->absPdgId()<<" status id: "<<par->status()<<endl;
//cout<<"Check1"<<endl;

      if (( (par->absPdgId() <7) || (par->absPdgId()==21)) && ((par->status()>21)&&(par->status()<30))){
        partoncount++;
      }

    }

    //Skiping events that has less than 3 parton jets
    if (partoncount<3){
      cout<<"Event skipped, less than 3 partons"<<endl;
      cout<<"# of partons in this event"<< partoncount<<endl;
      continue;
  //    continue;
    }

    
    
   
    std::cout<<"==================Event ID: "<<ievt<<"=================="<<std::endl;
    std::cout<<std::endl;
        
    //TODO: 1. skip this event if less than one photon
    //contii
    //TODO:  note the kinematics of the _leading photon
    //hint: ph->p4() ref (photons->size()<1) {
    //
    //         continue;
    //                }utnrs a TLorentzVector
    //TODO: read the "TrimmedAntiKt10TruthJets" contrainer from the XAOD
    //TODO: Loop through the jets and print out the pT and mass of the _leaing_et which is more than dealtaR >1.0 away formt he leading photon
    //hint: TLorenetzVector has a method DeltaR()
 //   cout<<"00000   PHOTON INFO   00000"<<endl;n
  //  std::cout<<photons->size()<<"photons in this event: " <<std::endl;
  
    // Discard events that has no photons. Otherwise print out leading photon kinematics

//DOING ALL THE ANTIKT4 CALCULATIONS HERE
  int checkKt=0;
  int kt4count=0;
  for (auto ktFour: *kt4s){
    kt4count++;
  }

  if (kt4count<3){
    cout<<"Less than 3 kt4 jets. Event skipped!"<<endl;
    continue;
  }



    for ( auto ktFour: *kt4s){
      TLorentzVector k=ktFour->p4();
        allkt4_pt.push_back(k.Pt());
        cout<<"kT4 PT: "<<k.Pt()<<endl;
      if (k.Pt()>subSubLeadKt4Pt){
        subSubLeadKt4Pt=k.Pt();
      }
      if (k.Pt()>subLeadKt4Pt){
        subSubLeadKt4Pt=subLeadKt4Pt;
        subLeadKt4Pt=k.Pt();
      }
      if (k.Pt()>leadKt4Pt){
        subLeadKt4Pt=leadKt4Pt;
        leadKt4Pt=k.Pt();
      }
      if (k.Pt()>330000){
        checkKt=1;
      } 
    }

  
    kt4_pt.push_back(leadKt4Pt);
    kt4_pt.push_back(subLeadKt4Pt);
    kt4_pt.push_back(subSubLeadKt4Pt);
  if (kt4s->size()<1){
    cout<<"There are no kt4s in the event. Event Skipped!"<<endl;
  }

  if (checkKt==1){
     for (auto ktFour: *kt4s){
      TLorentzVector k=ktFour->p4();
     //   TLorentzVector j= par;

        allkt4_pt_above330.push_back(k.Pt());
        if (k.Pt()>leadKt4PtAbv330){

          leadKt4PtAbv330=k.Pt();
      //Store jet_pt_below330 for all events with check=0(no jets above 330)
         }
      }
     kt4_pt_above330.push_back(leadKt4PtAbv330);
  }
  else {
    for (auto ktFour: *kt4s){
      TLorentzVector k=ktFour->p4();
      allkt4_pt_below330.push_back(k.Pt());
      if (k.Pt()>leadKt4PtBel330){
        leadKt4PtBel330=k.Pt();
      }
    }  
    kt4_pt_below330.push_back(leadKt4PtBel330);

  }

    int check=0; 
    int jetcount=0;
    
    float minJPt=1000000000;
   // cout<<"The number of particles in the event is: "<<particles->size()<<endl;
/*
    if (particles->size()<1) {
      cout<<"No particles! skipped!"<<endl;
       continue;
    }

*/


//Doing all the parton calculations here. 

    for (auto par: *particles){
      
   //   cout<<"pdgID: "<<par->absPdgId()<<" status id: "<<par->status()<<endl;


      if (( (par->absPdgId() <7) || (par->absPdgId()==21)) && ((par->status()>21)&&(par->status()<30))){
        jetcount++;
        TLorentzVector j=par->p4();
        jet.push_back(j);
        alljets_pt.push_back(j.Pt());
        

        float pt= j.Pt();
        if (pt< minJPt){
          minJPt=pt;
        }
        //store jet_pt for all particles that are jets with status()<=7
       // jet_pt.push_back(pt);
        if (pt>=330000){
          check=1;
        }
        //Print out the jets Pts
        cout<<"jet pT:"<<pt<<" status code: "<<par->status()<<endl;
      }

    }
    cout<<endl;
    
    minJetPt.push_back(minJPt);

  cout<<"There are "<<jetcount<<" of jets in the event."<<endl;
  cout<<"Min jet Pt of this event: "<<minJPt<<endl;

  if (jetcount<1) {
    cout<<"No jets! skipped!"<<endl;
    continue;
  }


  if (check==1){
    cout<<"There are at least one jets in this event that has pT of larger than 330." <<endl;
     for (auto j: jet){


        alljets_pt_above330.push_back(j.Pt());

     //   TLorentzVector j= par;
        if (j.Pt()>leadPtAbv330){
          leadPtAbv330=j.Pt();
      //Store jet_pt_below330 for all events with check=0(no jets above 330)
     
         }
       
     }


    jet_pt_above330.push_back(leadPtAbv330);

  }
  else {
    cout<<"There are no jet in this event that has PT of larger than 330."<<endl;
    for (auto j: jet){
      alljets_pt_below330.push_back(j.Pt()); 
  //      TLorentzVector j= par->p4();
      //Store jet_pt_below330 for all events with check=0(no jets above 330)
        if (j.Pt()>leadPtBel330){
        leadPtBel330=j.Pt();
        }

      
    }


    jet_pt_below330.push_back(leadPtBel330);

  }

 for (auto j: jet){
      
     //   TLorentzVector j= par->p4();
      //Store jet_pt_below330 for all events with check=0(no jets above 330)
        if (j.Pt()>leadPt){
          leadPt=j.Pt();
        }

 }


    jet_pt.push_back(leadPt);
  //Store a tree to keep track of how many events have at least one jet >330
  cut330.push_back(check);
  
    // IMPORTANT! save the latest entries in the tree
       
    
    output_tree->Fill();

        
  }
  output_file->Write();
  output_file->Close();
  
  return 0;
}
