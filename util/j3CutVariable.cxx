#include <iostream>
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODEventInfo/EventInfo.h"
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
  if (argc<4){
    std::cerr<< "Please specify the file name!"<< std ::endl;
    return 1;
  }


  string input_filename= argv[1];

  //j3 
  string input_filename2=argv[2];


  // TODO: add a 3rd argument to check if output file is specified and use that to name output file
  string output_filename = argv[3];
  
  std::cout <<"Reading from truth1 file: " << input_filename<<std::endl;
  
  //j3
  std::cout<<"Reading from truth3 file: " << input_filename2<<std::endl;

  xAOD::Init();


  TFile *input_file = TFile::Open(input_filename.c_str());
  //j3 
  TFile *input_file2 = TFile::Open(input_filename2.c_str());


  xAOD::TEvent* evt= new xAOD::TEvent(input_file);
  //j3
  xAOD::TEvent* evt2 = new xAOD::TEvent(input_file2);


  int nevt = evt->getEntries();
  //j3
  int nevt2 = evt2->getEntries();

  std::cout<<"There are "<<nevt <<" event in the truth1 input file." << std::endl;
  //j3
  std::cout<<"There are "<<nevt2<<" events in the truth3 input file."<< std::endl;


  std::cout<<"========================================================="<<std::endl;

  // important: create output file after loading the input file
  TFile* output_file = new TFile(output_filename.c_str(), "recreate");

  // create output tree after creating output file
  TTree* output_tree = new TTree("truth", "truth");

  // make vectors to hold variables
 // vector<float> cut330;
 
  int event_number;

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

//j3
  vector<float> jet10_pt;
  vector<float> jet10_m;

  vector<float> jet100_pt;
  vector<float> jet100_m;
  vector<float> jet100_tau1;
  vector<float> jet100_tau2;
  vector<float> jet100_tau21;

  vector<float> jet101_pt;
  vector<float> jet101_m;
  vector<float> jet101_tau1;
  vector<float> jet101_tau2;
  vector<float> jet101_tau21;

  vector<float> Ejet100_pt;
  vector<float> Ejet100_m;
  vector<float> Ejet100_tau1;
  vector<float> Ejet100_tau2;
  vector<float> Ejet100_tau21;

  vector<float> Ejet101_pt;
  vector<float> Ejet101_m;
  vector<float> Ejet101_tau1;
  vector<float> Ejet101_tau2;
  vector<float> Ejet101_tau21;

  vector<float> Ekt40_pt;
  vector<float> Ekt40_m;

  // add branches to the tree
  output_tree->Branch("event_number", &event_number, "event_number/I");

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

  //j3
  output_tree->Branch("jet10_pt", &jet10_pt);
  output_tree->Branch("jet10_m", &jet10_m);
    
  output_tree->Branch("jet100_pt", &jet100_pt);
  output_tree->Branch("jet100_m", &jet100_m);
  output_tree->Branch("jet100_tau1", &jet100_tau1);
  output_tree->Branch("jet100_tau2",&jet100_tau2);
  output_tree->Branch("jet100_tau21",&jet100_tau21);

  output_tree->Branch("jet101_pt", &jet101_pt);
  output_tree->Branch("jet101_m", &jet101_m);
  output_tree->Branch("jet101_tau1", &jet101_tau1);
  output_tree->Branch("jet101_tau2",&jet101_tau2);
  output_tree->Branch("jet101_tau21",&jet101_tau21);

  output_tree->Branch("Ejet100_pt", &Ejet100_pt);
  output_tree->Branch("Ejet100_m", &Ejet100_m);
  output_tree->Branch("Ejet100_tau1", &Ejet100_tau1);
  output_tree->Branch("Ejet100_tau2",&Ejet100_tau2);
  output_tree->Branch("Ejet100_tau21",&Ejet100_tau21);

  output_tree->Branch("Ejet101_pt", &Ejet101_pt);
  output_tree->Branch("Ejet101_m", &Ejet101_m);
  output_tree->Branch("Ejet101_tau1", &Ejet101_tau1);
  output_tree->Branch("Ejet101_tau2",&Ejet101_tau2);
  output_tree->Branch("Ejet101_tau21",&Ejet101_tau21);

  for (int ievt=0; ievt<nevt; ++ievt){
    

    evt->getEntry(ievt);
    //j3
    evt2->getEntry(ievt);

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

    jet10_pt.clear();
    jet10_m.clear();

    jet100_pt.clear();
    jet100_m.clear();

    jet101_pt.clear();
    jet101_m.clear();
    jet101_tau1.clear();
    jet101_tau2.clear();
    jet101_tau21.clear();

    Ejet100_pt.clear();
    Ejet100_m.clear();
    Ejet100_tau1.clear();
    Ejet100_tau2.clear();
    Ejet100_tau21.clear();

    Ejet101_pt.clear();
    Ejet101_m.clear();
    Ejet101_tau1.clear();
    Ejet101_tau2.clear();
    Ejet101_tau21.clear();

    Ekt40_pt.clear();
    Ekt40_m.clear();

    
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
    const xAOD::EventInfo* truth1EvntInfo(nullptr);
    evt->retrieve(truth1EvntInfo,"EventInfo");


    const xAOD::EventInfo* truth3EvntInfo(nullptr);
    evt2->retrieve(truth3EvntInfo,"EventInfo");

    const xAOD::TruthParticleContainer* particles(nullptr);
    evt->retrieve(particles, "TruthParticles");

    const xAOD::JetContainer *kt4s(nullptr);
    evt->retrieve(kt4s, "AntiKt4TruthJets");
  
    //j3
    const xAOD::JetContainer *jet10(nullptr);
    evt2->retrieve(jet10, "TrimmedAntiKt10TruthJets");

    float leadPt=0;
    float leadPtAbv330=0;
    float leadPtBel330=0;

    float leadKt4Pt=0;
    float subLeadKt4Pt=0;
    float subSubLeadKt4Pt=0;
    float leadKt4PtAbv330=0;
    float leadKt4PtBel330=0;

    //j3
    float leadKt10Pt=0;
    float subLeadKt10Pt=0;


//j4
    unsigned long long truth1EvtNumber=truth1EvntInfo->eventNumber();
    unsigned long long truth3EvtNumber=truth3EvntInfo->eventNumber();

    //Check
    if (truth1EvtNumber!=truth3EvtNumber){
      cerr << "ERROR! Event numbers do not match between TRUTH1 and TRUTH3. Aborting." << endl;
      return 1;
    }
    event_number = truth1EvtNumber;

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
    }

    //j3
    if (jet10->size()<2){
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
 int jet_idx=-1;
 TLorentzVector j100;
 TLorentzVector j101;
 TLorentzVector Ej100;
 TLorentzVector Ej101;
 
 const xAOD::Jet* j100_obj(nullptr);
 const xAOD::Jet* j101_obj(nullptr);
 const xAOD::Jet* Ej100_obj(nullptr);
 const xAOD::Jet* Ej101_obj(nullptr);

 float Ej100_tau1;
 float Ej100_tau2;
 float Ej100_tau21;

 float Ej101_tau1;
 float Ej101_tau2;
 float Ej101_tau21;

 bool goodjet=false;

 cout<<"Check2"<<endl;
 //j3 finding leading and subleading Akt10 jets
 for (auto j : *jet10){
   jet10_pt.push_back(j->pt());
   jet10_m.push_back(j->pt());

  if (j->pt()>j101.Pt()){
    j101=j->p4();
    j101_obj=j;
  }
  if (j->pt()>j100.Pt()){
    j101=j100;
    j101_obj=j100_obj;
    j100=j->p4();
    j100_obj=j;

  }
}
cout<<"check2.5"<<endl;
float j100_tau1 = j100_obj->auxdata<float>("Tau1_wta");
float j100_tau2 = j100_obj->auxdata<float>("Tau2_wta");
float j100_tau21 = (j100_tau1>0) ? (j100_tau2/j100_tau1): -1;
cout<<"check2.7"<<endl;
float j101_tau1 = j101_obj->auxdata<float>("Tau1_wta");
float j101_tau2 = j101_obj->auxdata<float>("Tau2_wta");
float j101_tau21 = (j101_tau1>0) ? (j101_tau2/j101_tau1): -1;


 cout<<"Check3"<<endl;

//j3 finding the jets after the kt10 

if (j100_tau21>j101_tau21){
  Ej100=j101;
  Ej101=j100;

  Ej100_obj=j101_obj;
  Ej101_obj=j100_obj;
}
else{
  Ej101=j101;
  Ej100=j100;

  Ej101_obj=j101_obj;
  Ej101_obj=j100_obj;
}

jet100_pt.push_back(j100.Pt());
jet100_m.push_back(j100.M());
jet100_tau1.push_back(j100_tau1);
jet100_tau2.push_back(j100_tau2);
jet100_tau21.push_back(j100_tau21);

jet101_pt.push_back(j101.Pt());
jet101_m.push_back(j101.M());
jet101_tau1.push_back(j101_tau1);
jet101_tau2.push_back(j101_tau2);
jet101_tau21.push_back(j101_tau21);

Ejet100_pt.push_back(Ej100.Pt());
Ejet100_m.push_back(Ej100.M());
Ejet100_tau1.push_back(Ej100_tau1);
Ejet100_tau2.push_back(Ej100_tau2);
Ejet100_tau21.push_back(Ej100_tau21);

Ejet101_pt.push_back(Ej101.Pt());
Ejet101_m.push_back(Ej101.M());
Ejet101_tau1.push_back(Ej101_tau1);
Ejet101_tau2.push_back(Ej101_tau2);
Ejet101_tau21.push_back(Ej101_tau21);


 cout<<"Check4"<<endl;

//j3: AntiKT4 calculations for the event selection j2

bool Egoodjet=false;
TLorentzVector Ekt40;

for (auto j :*kt4s) {
  if ((Ej100.DeltaR(j->p4())>1.0) && (j->pt()>Ekt40.Pt())){
    Ekt40= j->p4();
    Egoodjet=true;
  }
}
  if (Egoodjet==false){
  cout<<"j2: Event skipped, there are no jet outside of the exclusion region."<<endl;
  continue;
  }

Ekt40_pt.push_back(Ekt40.Pt());
Ekt40_m.push_back(Ekt40.M());


 cout<<"Check5"<<endl;

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

 cout<<"Check6"<<endl;

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


 cout<<"Check7"<<endl;

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

 cout<<"Check8"<<endl;

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
