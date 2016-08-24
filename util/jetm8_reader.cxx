#include <iostream>
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODJet/JetContainer.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include <math.h>
#include <cmath>
#include "OutputTree/OutputTree.h"
#include "ZprimexAOD/Retrieve.h"
#include <algorithm>

using namespace std;

/*
 * Function to compare the pT of xAOD particle objects; useful for sorting.
 */
bool pt_compare(const xAOD::IParticle* p1, const xAOD::IParticle* p2){
  return p1->pt() > p2->pt();
}

float get_tau21(const xAOD::Jet* j) {
  float tau1 = j->auxdata<float>("Tau1_wta");
  if (tau1 == 0) {
    return 0;
  }
  float tau2 = j->auxdata<float>("Tau2_wta");
  return tau2 / tau1;
}


float get_D2(const xAOD::Jet* j) {
  float ECF2 = j->auxdata<float>("ECF2");

  if (ECF2 == 0) {
    return 0;
  }
  float ECF1 = j->auxdata<float>("ECF1");
  float ECF3 = j->auxdata<float>("ECF3");

  return (ECF3 * ECF1*ECF1*ECF1)/(ECF2 *ECF2 * ECF2);
}

pair<float, float> get_event_counts(xAOD::TEvent* evt) {
  const xAOD::CutBookkeeperContainer* cuts(nullptr);
  evt->retrieveMetaInput(cuts, "CutBookkeepers");

  map<int, vector<const xAOD::CutBookkeeper*> > cycles;
  for (auto cb : *cuts) {
      cycles[cb->cycle()].push_back(cb);
  }

  for (auto c : cycles) {
      float nevt, nevt_wt;
      auto cbks = c.second;
      if (cbks.size() < 2) continue;
      bool goes_to_jetm8 = false;
      for (auto cb : cbks) {
          auto outputs = cb->outputStreams();
          goes_to_jetm8 = ( find(begin(outputs), end(outputs), "StreamDAOD_JETM8") != end(outputs) );
          if (cb->name() == "AllExecutedEvents") {
              nevt = cb->nAcceptedEvents();
              nevt_wt = cb->sumOfEventWeights();
          }
      }
      if (goes_to_jetm8) {
          cout << "Found AOD stream to JETM8 in cycle " << c.first << endl;
          return make_pair(nevt, nevt_wt);
      }
  }

  return make_pair(0., 0.);
}

/*
 * Process a single event of the input tree. If the event passes,
 * the result is written to the given OutputTree.
 * Assumes the TEvent is already pointing at the event of interest.
 */
void process_event(xAOD::TEvent* evt, OutputTree* output_tree) {
  const xAOD::EventInfo* evt_info(nullptr);
  evt->retrieve(evt_info,"EventInfo");

  float event_weight = evt_info->mcEventWeight();
  output_tree->add_scalar("event_weight", event_weight);

  auto jet4 = CopyRetrieve<xAOD::Jet>(evt,"AntiKt4LCTopoJets");
  sort(jet4.begin(), jet4.end(), pt_compare);
  output_tree->add_jets("jet",jet4);

  auto jet10 = CopyRetrieve<xAOD::Jet>(evt,"AntiKt10LCTopoTrimmedPtFrac5SmallR20Jets");
  sort(jet10.begin(),jet10.end(),pt_compare);

  output_tree->add_jets("fjet",jet10);

  vector<float> fatjet_tau21;
  vector<float> fatjet_D2;
  for (auto j : jet10) {
    fatjet_tau21.push_back(get_tau21(j));
    fatjet_D2.push_back(get_D2(j));
  }
  output_tree->add_vector("fjet_tau21", fatjet_tau21);
  output_tree->add_vector("fjet_D2", fatjet_D2);


  if (jet10.size()<2){
    cout<<"Event skipped, less than 2 AKT10 jets."<<endl;
    return;
  }


  // sort the leading/subleading Akt10 jets in order of
  // their tau21 value
  vector<const xAOD::Jet* > jet10_tau;
  if (fatjet_tau21[0] < fatjet_tau21[1]) {
    jet10_tau.push_back(jet10[0]);
    jet10_tau.push_back(jet10[1]);
    output_tree->add_vector("fjetTau_tau21", {fatjet_tau21[0], fatjet_tau21[1]});
    output_tree->add_vector("fjetTau_D2", {fatjet_D2[0], fatjet_D2[1]});
  }
  else {
    jet10_tau.push_back(jet10[1]);
    jet10_tau.push_back(jet10[0]);
    output_tree->add_vector("fjetTau_tau21", {fatjet_tau21[1], fatjet_tau21[0]});
    output_tree->add_vector("fjetTau_D2", {fatjet_D2[1], fatjet_D2[0]});
  }
  // write the tau21-sorted jets to the tree
  output_tree->add_jets("fjetTau", jet10_tau);


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

  output_tree->Fill();
}


void usage(int /*argc*/, char** argv) {
  cout << "Usage: " << argv[0] << " input_file [input_file, ...] output_file" << endl;
}

int main(int argc, char** argv){
  if (argc<3){
    std::cerr<< "Please specify the file name!"<< std ::endl;
    usage(argc, argv);
    return 1;
  }

  if (string(argv[1]) == "-h" || string(argv[1]) == "--help") {
    usage(argc, argv);
    return 0;
  }

  vector<string> input_filenames;
  cout << "Input files: " << endl;
  for (int iarg = 1; iarg < (argc-1); ++iarg) {
    input_filenames.push_back(argv[iarg]);
    cout << "  " << argv[iarg] << endl;
  }

  string output_filename = argv[argc-1];
  cout << "Output file: " << endl;
  cout << "  " << output_filename << endl;

  xAOD::Init();

  TFile* output_file = new TFile(output_filename.c_str(), "create");
  if (not output_file->IsOpen()) {
    cerr << "Could not open output file! Abort." << endl;
    return 1;
  }
  output_file->cd();
  OutputTree *output_tree= new OutputTree("zpj");

  TH1F* h_nevt_total = new TH1F("nevt_total", "nevt_total", 1, 0, 1);
    TH1F* h_nevt_total_wt = new TH1F("nevt_total_wt", "nevt_total_wt", 1, 0, 1);

  xAOD::TEvent* evt= new xAOD::TEvent();

  for (string filename : input_filenames) {
    TFile *input_file = TFile::Open(filename.c_str());

    evt->readFrom(input_file);

    float nevt_total, nevt_total_wt;
    tie(nevt_total, nevt_total_wt) = get_event_counts(evt);

    if ( (nevt_total<0) || (nevt_total_wt<0) ) {
        cerr << "Invalid event counts found!!! Skipping file." << endl;
        continue;
    }

    h_nevt_total->Fill(0., nevt_total);
    h_nevt_total_wt->Fill(0., nevt_total_wt);
    cout << "Total events:    " << nevt_total << endl;
    cout << "Weighted events: " << nevt_total_wt << endl;

    int nevt = evt->getEntries();
    std::cout<<"There are "<< nevt <<" event in this jetm8 file." << std::endl;

    for (int ievt=0; ievt<nevt; ++ievt){
      //std::cout<<"==================Event: "<<ievt<<"=================="<<std::endl;
      //std::cout<<std::endl;

      evt->getEntry(ievt);
      output_tree->clear();

      process_event(evt, output_tree);
    }

    input_file->Close();
  }

  output_file->cd();
  output_tree->Write();
  h_nevt_total->Write();
  h_nevt_total_wt->Write();
  output_file->Close();

  return 0;
}
