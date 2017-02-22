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
#include <utility>
#include <LHAPDF/LHAPDF.h>
#include "LHAPDF/Reweighting.h"
#include "xAODTruth/TruthEvent.h"
#include "xAODTruth/versions/TruthEvent_v1.h"


using namespace std;

/*
 * Function to compare the pT of xAOD particle objects; useful for sorting.
 */
bool pt_compare(const xAOD::IParticle* p1, const xAOD::IParticle* p2){
  return p1->pt() > p2->pt();
}
/*
bool deltaR_compare(const xAOD::IParticle* p1, const xAOD::IParticle* p2, const xAOD::IParticle* pStar, int i){
  pStar= ut_jet.at(i);
  return pStar->p4().DeltaR(p1->p4())< pStar->p4().DeltaR(p2->p4());
}
*/
/*
bool deltaR_compare(const xAOD::IParticle* p1, const xAOD::IParticle* p2, const xAOD::IParticle* pStar, int i){
  pStar= ut_jet.at(i);
  return pStar->p4().DeltaR(p1->p4())< pStar->p4().DeltaR(p2->p4());
}
*/
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
void process_event(xAOD::TEvent* evt, OutputTree* output_tree, LHAPDF::PDF* p0, LHAPDF::PDF* p1, LHAPDF::PDF* p2) {

  cout<<"check1"<<endl;

  const xAOD::EventInfo* evt_info(nullptr);
  evt->retrieve(evt_info,"EventInfo");

  //const xAOD::TruthEvent_v1::PdfInfo *pdf_info;
  //evt->retrieve(pdf_info, "PdfInfo");

  const xAOD::TruthEvent_v1::PdfInfo *pdf_info=&xAOD::TruthEvent_v1::pdfInfo();

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

 //Retrieve and sort and add ut_fjet10
  auto ut_jet10 = CopyRetrieve<xAOD::Jet>(evt,"AntiKt10LCTopoJets");
  sort(ut_jet10.begin(), ut_jet10.end(), pt_compare);
  output_tree->add_jets("ut_fjet10",ut_jet10);


  //Kill the program if ut_jet size doesn't match that of jet10 in an event.
  if (ut_jet10.size()!=jet10.size()){
      cout<<"Unmatched event number between jet10 and untrimmed jet10!"<<endl;
      return;
    }

cout<<"Check2"<<endl;

//Adding in the additional CMS trigger requirements   //Need to move this to before the tree items are written in. 
double HT;
double largestPT=0.;
double largestTrimmedMass=0.;
for (int i ; i < jet10.size(); i++){
	HT = HT+jet10.at(i)->p4().Pt();
	if (largestPT < jet10.at(i)->p4().Pt()){
		largestPT = jet10.at(i)->p4().Pt();
	}	
	if (largestTrimmedMass < ut_jet10.at(i)->p4().M()){
		largestTrimmedMass = ut_jet10.at(i)->p4().M();
	}

}
if (((HT>650e3) && (largestTrimmedMass>50e3)) || ((largestPT >350e3) && (largestTrimmedMass > 30e3)) || (HT > 800e3)){
	cout<<"Event passes the CMS trigger"<<endl;

}
else{
	cout<<"Event does not pass the CMS Trigger"<<endl;
return;
}

    vector<int> assoInd;
    for (int i=0; i<jet10.size(); i++){
      vector<pair<double,int> > deltaR_ind;
      for (int j=0; j<ut_jet10.size(); j++){

        
          double deltaR=jet10.at(i)->p4().DeltaR(ut_jet10.at(j)->p4());
          deltaR_ind.emplace_back(deltaR,j);
        }
        sort(deltaR_ind.begin(),deltaR_ind.end());
        //Create an vector of associate index
        assoInd.push_back(deltaR_ind.at(0).second);
      }
    

//cout<<"Check3"<<endl;

//need to still remove the overlapping candidates.


//finding the variable for the asso jets
 	  vector<float> ut_asso_fatjet_Jvt;
	  vector<float> ut_asso_fatjet_tau21;
 	  vector<float> ut_asso_fatjet_D2;
	  vector<float> ut_asso_fatjet_pt;
	  vector<float> ut_asso_fatjet_m;
	  for (auto i : assoInd){
	    ut_asso_fatjet_tau21.push_back(get_tau21(ut_jet10.at(i)));
	    ut_asso_fatjet_D2.push_back(get_D2(ut_jet10.at(i)));
	    ut_asso_fatjet_pt.push_back(ut_jet10.at(i)->p4().Pt());
	    ut_asso_fatjet_m.push_back(ut_jet10.at(i)->p4().M());
	    ut_asso_fatjet_Jvt.push_back(ut_jet10.at(i)->auxdata<float>("Jvt"));
	  }

	  output_tree->add_vector("ut_asso_fjet_tau21", ut_asso_fatjet_tau21);
	  output_tree->add_vector("ut_asso_fjet_D2",ut_asso_fatjet_D2);
	  output_tree->add_vector("ut_asso_fjet_pt", ut_asso_fatjet_tau21);
	  output_tree->add_vector("ut_asso_fjet_m",ut_asso_fatjet_m);
	  output_tree->add_vector("ut_asso_fjet_Jvt", ut_asso_fatjet_Jvt);




	  


//***********************Edit later**************************//
/*
  vector<float> rhoPrime;
  vector<float> ut_tau21;

  for (auto i: assoInd){
      fatjet_tau21.push_back(get_tau21(j));
      rp=fjet
      rhoPrime.push_back()
  }
 // vector<float
*/

cout<<"check4"<<endl;

 // outputtree->add_vector("rho",log(
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

cout<<"check5"<<endl;
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

 //pdf weight calculation 
  
  double pdfwVsCt4 = LHAPDF::weightxxQ2(pdf_info->pdgId1, pdf_info->pdgId2, pdf_info->x1, pdf_info->x2, pdf_info->Q, p0, p1);
  double pdfwVsMMHT = LHAPDF::weightxxQ2(pdf_info->pdgId1, pdf_info->pdgId2, pdf_info->x1, pdf_info->x2, pdf_info->Q, p0, p2);
 

  output_tree->add_scalar("pdfweightVS_Ct4", pdfwVsCt4);
  output_tree->add_scalar("pdfweightVS_MMHT", pdfwVsMMHT);
  output_tree->Fill();
}

pair<vector<string>, string> parse_file_list(int argc, char** argv) {
  vector<string> input_filenames;

  // load in the list of input files
  string firstfile = argv[1];
  if (firstfile.find(",") != string::npos) {
    // if the first argument contains commas, split it up and treat
    // it as a list of files
    size_t pos = 0;
    while ( (pos = firstfile.find(",")) != string::npos) {
      input_filenames.push_back(firstfile.substr(0, pos));
      cout<<"check6"<<endl;

      firstfile.erase(0, pos + 1);
    }
    if (firstfile.size() > 0) {
      input_filenames.push_back(firstfile);
    }
  }
  else {
    // otherwise, read the input files off the command line
    for (int i = 1; i < (argc-1); ++i) {
      input_filenames.push_back(argv[i]);
    }
  }
  // in any case, assume the last argument is the output filename.
  string output_filename = argv[argc-1];

  return make_pair(input_filenames, output_filename);
  cout<<"check7"<<endl;


}
void usage(int /*argc*/, char** argv) {
  cout << "Usage: " << argv[0] << " input_file [input_file ...] output_file" << endl;
  cout << "       " << argv[0] << " input_file[,input_file,...] output_file" << endl;
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

  // load up the input and output file names
  vector<string> input_filenames;
  string output_filename;
  tie(input_filenames, output_filename) = parse_file_list(argc, argv);

  cout << "Input files: " << endl;
  for (auto fn : input_filenames) {
      cout << "  " << fn << endl;
  }
  cout << "Output file: " << endl;
  cout << "  " << output_filename << endl;

  // initialize xAOD framework
  xAOD::Init();
  cout<<"check main 1"<<endl;

  // open up our output file (or throw a fit if it already exists)
  TFile* output_file = new TFile(output_filename.c_str(), "RECREATE");
  if (not output_file->IsOpen()) {
    cerr << "Could not open output file! Abort." << endl;
    return 1;
  }
  cout<<"check main 2"<<endl;

  output_file->cd();
  OutputTree *output_tree= new OutputTree("zpj");

  cout<<"check main 3"<<endl;

  TH1F* h_nevt_total = new TH1F("nevt_total", "nevt_total", 1, 0, 1);
    TH1F* h_nevt_total_wt = new TH1F("nevt_total_wt", "nevt_total_wt", 1, 0, 1);

  xAOD::TEvent* evt= new xAOD::TEvent();

  cout<<"check main 4"<<endl;

  for (string filename : input_filenames) {
    cout<<"Check main 4.1.1"<<endl;
    TFile *input_file = TFile::Open(filename.c_str());

    cout<<"check main 4.1"<<endl;

    evt->readFrom(input_file);

    cout<<"check main 4.2"<<endl;

    float nevt_total, nevt_total_wt;
    tie(nevt_total, nevt_total_wt) = get_event_counts(evt);

    cout<<"check main 4.3"<<endl;

    if ( (nevt_total<0) || (nevt_total_wt<0) ) {
        cerr << "Invalid event counts found!!! Skipping file." << endl;
        continue;
    }

    cout<<"check main 5"<<endl;

    h_nevt_total->Fill(0., nevt_total);
    h_nevt_total_wt->Fill(0., nevt_total_wt);
    cout << "Total events:    " << nevt_total << endl;
    cout << "Weighted events: " << nevt_total_wt << endl;

    cout<<"check main 6"<<endl;

    int nevt = evt->getEntries();
    std::cout<<"There are "<< nevt <<" event in this jetm8 file." << std::endl;

//adding the lha library
    LHAPDF::PDF* p0 = LHAPDF::mkPDF("NNPDF30_lo_as_0130/0");
    LHAPDF::PDF* p1 = LHAPDF::mkPDF("CT14lo/0");
    LHAPDF::PDF* p2 = LHAPDF::mkPDF("MMHT2014lo68cl/0"); 
//Event loop
    for (int ievt=0; ievt<nevt; ++ievt){
      //std::cout<<"==================Event: "<<ievt<<"=================="<<std::endl;
      //std::cout<<std::endl;


      evt->getEntry(ievt);
      output_tree->clear();

      process_event(evt, output_tree,p0,p1,p2);
    }
    cout<<"check main 7"<<endl;

    input_file->Close();
  }

  output_file->cd();
  output_tree->Write();
  h_nevt_total->Write();
  h_nevt_total_wt->Write();
  output_file->Close();

  return 0;
}
