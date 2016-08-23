#include <iostream>
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODJet/JetContainer.h"
#include "xAODEventInfo/EventInfo.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <math.h>
#include <cmath>
#include "OutputTree/OutputTree.h"
#include "truthReader/Retrieve.h"
#include <algorithm>

using namespace std;
//Reading an input file and a treefile with a branch where it can be determined whether the event contain at least one particle of >330 GeV

bool pt_compare(const xAOD::IParticle* p1, const xAOD::IParticle* p2){
  return p1->pt() >p2->pt();
}

void usage(int argc, char** argv) {
  cout << "Usage: " << argv[0] << " input_file output_file" << endl;
}

int main(int argc, char** argv){
  if (argc<3){
    std::cerr<< "Please specify the file name!"<< std ::endl;
    usage(argc, argv);
    return 1;
  }


  string input_filename= argv[1];

  if (input_filename == "-h" || input_filename == "--help") {
    usage(argc, argv);
    return 0;
  }

  //j3 
  //  string input_filename2=argv[2];


  // TODO: add a 3rd argument to check if output file is specified and use that to name output file
  string output_filename = argv[2];

  std::cout <<"Reading from jetm8  file: " << input_filename<<std::endl;

  //j3
  // std::cout<<"Reading from truth3 file: " << input_filename2<<std::endl;

  xAOD::Init();


  TFile *input_file = TFile::Open(input_filename.c_str());
  //j3 
  // TFile *input_file2 = TFile::Open(input_filename2.c_str());


  xAOD::TEvent* evt= new xAOD::TEvent(input_file);


  int nevt = evt->getEntries();
  //j3

  std::cout<<"There are "<<nevt <<" event in the jetm8  input file." << std::endl;
  //j3
  // std::cout<<"There are "<<nevt2<<" events in the truth3 input file."<< std::endl;


  std::cout<<"========================================================="<<std::endl;

  // important: create output file after loading the input file
  TFile* output_file = new TFile(output_filename.c_str(), "recreate");

  // create output tree after creating output file
  //  TTree* output_tree = new TTree("truth", "truth");
  OutputTree *output_tree= new OutputTree("truth");

  // make vectors to hold variables
  // vector<float> cut330;

  int event_number;

  // add branches to the tree
  for (int ievt=0; ievt<nevt; ++ievt){


    evt->getEntry(ievt);
    //j3


    output_tree->clear();

    //Declaration of variable (each event)
    /*
    if (ievt % 500 ==0) {
    std::cout << "Processing event #" << ievt << "/" << nevt <<std::endl;
    }
    */ 
    //   const xAOD::TruthParticleContainer* photons(nullptr);
    //Grabing the pdgID of the particle  
    //Grabing all particles 
    //
    const xAOD::EventInfo* evt_info(nullptr);
    evt->retrieve(evt_info,"EventInfo");


    //const xAOD::JetContainer *kt4s(nullptr);
    //evt->retrieve(kt4s, "AntiKt4TruthJets");
    auto jet4 = CopyRetrieve<xAOD::Jet>(evt,"AntiKt4LCTopoJets");
    sort(jet4.begin(), jet4.end(), pt_compare);
    output_tree->add_jets("jet",jet4);

    //j4
    auto jet10 = CopyRetrieve<xAOD::Jet>(evt,"AntiKt10LCTopoTrimmedPtFrac5SmallR20Jets");
    sort(jet10.begin(),jet10.end(),pt_compare);
    output_tree->add_jets("fjet",jet10);



    //j4
    unsigned long long event_number=evt_info->eventNumber();

    //j3
    if (jet10.size()<2){
    cout<<"Event skipped, less than 2 AKT10 jets."<<endl;
    continue;
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
    //j3 all Antikt10 calculations are done here.

    cout<<"Check1"<<endl;


    // calculate tau21 for the leading/subleading Akt10 jets
    float j100_tau1 = jet10[0]->auxdata<float>("Tau1_wta");
    float j100_tau2 = jet10[0]->auxdata<float>("Tau2_wta");
    float j100_tau21 = (j100_tau1>0) ? (j100_tau2/j100_tau1): -1;
    cout<<"Check2"<<endl;
    float j101_tau1 = jet10[1]->auxdata<float>("Tau1_wta");
    float j101_tau2 = jet10[1]->auxdata<float>("Tau2_wta");
    float j101_tau21 = (j101_tau1>0) ? (j101_tau2/j101_tau1): -1;


    cout<<"Check3"<<endl;

    // sort the leading/subleading Akt10 jets in order of
    // their tau21 value
    vector<const xAOD::Jet* > jet10_tau;
    if (j100_tau21 < j101_tau21) {
      jet10_tau.push_back(jet10[0]);
      jet10_tau.push_back(jet10[1]);
    }
    else {
      jet10_tau.push_back(jet10[1]);
      jet10_tau.push_back(jet10[0]);
    }
    // write the tau21-sorted jets to the tree
    output_tree->add_jets("fjetTau", jet10_tau);

    cout<<"Check4"<<endl;


    // find small-radius jets that are not overlapping with the
    // candidate antikt10 jets
    std::vector<const xAOD::Jet* > jet4_nonoverlap;
    std::vector<const xAOD::Jet* > jet4_nonoverlap_tau;

    TLorentzVector leading_candidate = jet10[0]->p4();
    TLorentzVector resonance_candidate = jet10_tau[0]->p4();
    for (auto j : jet4) {
      if (leading_candidate.DeltaR(j->p4()) > 1.0) {
        jet4_nonoverlap.push_back(j);
      }
      if (resonance_candidate.DeltaR(j->p4()) > 1.0) {
        jet4_nonoverlap_tau.push_back(j);
      }
    }
    // write the nonoverlapping jets to the tree
    output_tree->add_jets("jetNOleading", jet4_nonoverlap);
    output_tree->add_jets("jetNOtau", jet4_nonoverlap_tau);


    cout<<"Check8"<<endl;



    output_tree->Fill();

  }
  output_file->Write();
  output_file->Close();

  return 0;
}
