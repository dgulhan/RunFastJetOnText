#include "QR.h"
#include <TCanvas.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TF1.h>
#include <TFile.h>
// using namespace fastjet;
using namespace std;

// g++ `root-config --cflags` runQR.cc -o runQR `fastjet-install/bin/fastjet-config --cxxflags --libs --plugins` `root-config --libs`
      
int main(int argc, char *argv[]) { 
 
 QR * file = new QR(argv[1]);
 vector <Event> event_vector = file->get_event_vector();

 TFile* outf = new TFile(argv[2],"recreate");
 
 TNtuple* nt_dijet = new TNtuple("nt_dijet","","phi1:eta1:pt1:phi2:eta2:pt2:phi3:eta3:pt3:dphi:QG1:QG2:QG3");
 TNtuple* nt_track = new TNtuple("nt_track","","jtpt:jtphi:jteta:pt:eta:phi:QG");
 TNtuple* nt_jet = new TNtuple("nt_jet","","eta:phi:pt:ntrk:QG");

 for(vector <Event>::iterator it = event_vector.begin(); it != event_vector.end(); ++it){
  std::vector<std::pair<float, std::pair<float, std::pair <float, int> > > > jets;
  int njet=0;

  for(int i=0;i<it->get_number_of_jets();i++){
   Jet * Qjet = it->get_ieth_jet(i);
   if(fabs(Qjet->get_eta())>2) continue;
   if(Qjet->get_pt()<10) continue;
   jets.push_back(std::make_pair(Qjet->get_pt(),std::make_pair(Qjet->get_eta(), std::make_pair(Qjet->get_phi(),Qjet->get_QG()))));
   njet++;
   int ntrk=0;
   for(int j=0;j<it->get_number_of_particles();j++){
    Particle * particle = it->get_ieth_particle(j);
    nt_track->Fill(Qjet->get_pt(),Qjet->get_phi(),Qjet->get_eta(),particle->get_pt(),particle->get_eta(),particle->get_phi(),Qjet->get_QG());
    if(sqrt((particle->get_phi()-Qjet->get_phi())*(particle->get_phi()-Qjet->get_phi())+(particle->get_eta()-Qjet->get_eta())*(particle->get_eta()-Qjet->get_eta()))<0.3){
      ntrk++;
    }
   }
   nt_jet->Fill(Qjet->get_eta(),Qjet->get_phi(),Qjet->get_pt(),ntrk,Qjet->get_QG());
  }
  std::sort(jets.begin(),jets.end());
  
  float pt1=-99;
  float pt2=-99;
  float pt3=-99;
  float phi1=-99;
  float phi2=-99;
  float phi3=-99;
  float eta1=-99;
  float eta2=-99;
  float eta3=-99;
  float QG1=0;
  float QG2=0;
  float QG3=0;

  if(njet>0){
   pt1=jets[njet-1].first;
   eta1=jets[njet-1].second.first;
   phi1=jets[njet-1].second.second.first;
   QG1=jets[njet-1].second.second.second;
   if(njet>1){
    pt2=jets[njet-2].first;
    eta2=jets[njet-2].second.first;
    phi2=jets[njet-2].second.second.first;
    QG2=jets[njet-2].second.second.second;
    if(njet>2){
     pt3=jets[njet-3].first;
     eta3=jets[njet-3].second.first;
     phi3=jets[njet-3].second.second.first;
     QG3=jets[njet-3].second.second.second;
    }
   }
  }
  nt_dijet->Fill(phi1,eta1,pt1,phi2,eta2,pt2,phi3,eta3,pt3,acos(cos(phi1-phi2)),QG1,QG2,QG3);
 }
 outf->Write();
 outf->Close();
}
