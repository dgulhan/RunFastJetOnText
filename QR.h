#include "fastjet/ClusterSequence.hh"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
using namespace fastjet;
using namespace std;

struct Particle
{
 private:
  double phi;
  double pt;
  double eta;
  double pt_original;
 public:
  Particle(double pt, double eta, double phi,double pt_original)
  {
   this->eta = eta;
   this->phi = phi;
   this->pt = pt;
   this->pt = pt;
   this->pt_original=pt_original;
   // cout<<pt_original<<endl;
  }
  double get_eta(){return eta;}
  double get_phi(){return phi;}
  double get_pt(){return pt;}
  double get_pt_original(){return pt_original;}
};

struct Jet
{
 private:
  double pt;
  double eta;
  double phi;
  int QG;
  int number_of_incones;
  Particle **incones;
 public:
   void reset()
  { 
   pt=0;
   eta=-99;
   phi=-99;
   QG=0;
   incones=NULL;
   number_of_incones=0;
  }
  Jet(double pt,double phi,double eta, int QG)
  {
   reset();
   this->pt=pt;
   this->phi=phi;
   this->eta=eta;
   this->QG=QG;
  }
  double get_pt(){return pt;}
  double get_phi(){return phi;}
  double get_eta(){return eta;}
  double get_QG(){return QG;}
  int get_number_of_incones(){return number_of_incones;}
  void add_incone(double pt, double eta, double phi,double pt_original)
  {
   Particle *incone = new Particle(pt,eta,phi,pt_original);
   if(incones == NULL)
   {
    incones = (Particle **)calloc(2,sizeof(Particle *));
    incones[0] = incone;
  	incones[1] = NULL;
   }
   else
   {
    incones = (Particle **)realloc(incones, (number_of_incones+2)*sizeof(Particle *));
	  incones[number_of_incones] = incone;
	  incones[number_of_incones+1] = NULL;
   }
   number_of_incones++;
  }
  Particle * get_ieth_incone(int i){return incones[i];}
};

struct Event
{
 private:
  Jet **jets;//!SET BY READING FROM THE DATA FILE
  Jet **Qjets;//!SET BY READING FROM THE DATA FILE
  Particle **particles;
  Particle **Qparticles;
  int number_of_jets;
  int number_of_Qjets;
  int number_of_particles;//!SAME AS NF BUT IS AUTOMATICALLY SET WHEN ADDING THE FRAGMENTS (see add_particle)
  int number_of_Qparticles;//!SAME AS NF BUT IS AUTOMATICALLY SET WHEN ADDING THE FRAGMENTS (see add_particle)

 public:
  void reset()
  {
   particles = NULL;
   Qparticles = NULL;
   jets = NULL;
   Qjets = NULL;
   number_of_jets=0;
   number_of_Qjets=0;
   number_of_particles=0;
   number_of_Qparticles=0;
  }
  
  Event()
  {
   reset();
  }
  
  ~Event()
  {
   
  }
  
  void add_particle(double pt, double eta, double phi)
  { 
   Particle *particle = new Particle(pt,eta,phi,pt);
   if(particles == NULL)
   {
    particles = (Particle **)calloc(2,sizeof(Particle *));
    particles[0] = particle;
	  particles[1] = NULL;
   }
   else
   {
    particles = (Particle **)realloc(particles, (number_of_particles+2)*sizeof(Particle *));
	  particles[number_of_particles] = particle;
	  particles[number_of_particles+1] = NULL;
   }
   number_of_particles++;
  }
  
  void add_Qparticle(double pt, double eta, double phi,double pt_original)
  { 
   Particle *particle = new Particle(pt,eta,phi,pt_original);
   if(Qparticles == NULL)
   {
    Qparticles = (Particle **)calloc(2,sizeof(Particle *));
    Qparticles[0] = particle;
	  Qparticles[1] = NULL;
   }
   else
   {
    Qparticles = (Particle **)realloc(Qparticles, (number_of_Qparticles+2)*sizeof(Particle *));
	  Qparticles[number_of_Qparticles] = particle;
	  Qparticles[number_of_Qparticles+1] = NULL;
   }
   number_of_Qparticles++;
  }
  
  void add_jet(Jet *jet)
  {
   if(jets == NULL)
   {
    jets = (Jet **)calloc(2,sizeof(Jet *));
    jets[0] = jet;
	  jets[1] = NULL;
   }
   else
   {
    jets = (Jet **)realloc(jets, (number_of_jets+2)*sizeof(Jet *));
	  jets[number_of_jets] = jet;
  	jets[number_of_jets+1] = NULL;
   }
   number_of_jets++;
  }
  int get_number_of_particles(){return number_of_particles;}
  int get_number_of_Qparticles(){return number_of_particles;}
  int get_number_of_jets(){return number_of_jets;}
  Particle * get_ieth_particle(int i){return particles[i];}  
  Particle * get_ieth_Qparticle(int i){return Qparticles[i];}  
  Jet * get_ieth_jet(int i){return jets[i];}  
  
  void add_Qjet(Jet *jet)
  {
   if(Qjets == NULL)
   {
    Qjets = (Jet **)calloc(2,sizeof(Jet *));
    Qjets[0] = jet;
	  Qjets[1] = NULL;
   }
   else
   {
    Qjets = (Jet **)realloc(Qjets, (number_of_Qjets+2)*sizeof(Jet *));
	  Qjets[number_of_Qjets] = jet;
  	Qjets[number_of_Qjets+1] = NULL;
   }
   number_of_Qjets++;
  }
  int get_number_of_Qjets(){return number_of_Qjets;}
  Jet * get_ieth_Qjet(int i){return Qjets[i];}  
};

class QR
{
 private:
   vector <Event> events;
   
 public:
  vector <Event> get_event_vector()
  {
   return events;
  } 
  QR(char* file_path) {
   double R = 0.3;
   JetDefinition jet_def(antikt_algorithm, R);
   // cout<<"file: "<<file_path.data()<<endl;
   cout<<"file: "<<file_path<<endl;
   ifstream input_file(file_path);
   // ifstream input_file(file_path.data());
   // int state=0;
   string line;
   vector<PseudoJet> particles;
   vector<PseudoJet> Qparticles;
   Event current;
   // int iparticle=0;

   while(getline(input_file,line))
   {
    if(line.length()==0) continue;
    // if(state==0 and string::npos == line.find("END"))
    if(string::npos == line.find("END"))
    {
     int id,QG;
     double px, py, pz, E, Qpx, Qpy, Qpz, QE;
     stringstream line_s;
     line_s << line;
     line_s >> id >> px >> py >> pz >> E >> Qpx >> Qpy >> Qpz >> QE >> QG;
     fastjet::PseudoJet Particle(px,py,pz,E);
     Particle.set_user_index(QG);
     particles.push_back(Particle); 
     fastjet::PseudoJet QParticle(Qpx,Qpy,Qpz,QE);
     QParticle.set_user_index(QG);
     Qparticles.push_back(QParticle); 
     // cout<<"QG in assignment="<<Particle.user_index()<<" pT "<<sqrt(pow(px,2)+pow(py,2))<<endl;
    }
    // if(state==0 and string::npos !=line.find("END")){
    if(string::npos !=line.find("END")){
     ClusterSequence cs(particles, jet_def);
     ClusterSequence Qcs(Qparticles, jet_def);
     for(unsigned i = 0; i < particles.size(); i++){
      current.add_particle(particles[i].perp(),particles[i].eta(),particles[i].phi());
     }
     for(unsigned i = 0; i < Qparticles.size(); i++){
      // cout<<particles[i].perp();
      current.add_Qparticle(Qparticles[i].perp(),Qparticles[i].eta(),Qparticles[i].phi(),particles[i].perp());
     }
     vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
     vector<PseudoJet> Qjets = sorted_by_pt(Qcs.inclusive_jets());
     for(unsigned i = 0; i < Qjets.size(); i++){

      vector<PseudoJet> constituents = Qjets[i].constituents();
	    int QG=0;
      for (unsigned j = 0; j < constituents.size(); j++) {
       for(int k=0;k<Qparticles.size();k++){
        if(Qparticles[k].perp()==constituents[j].perp() && Qparticles[k].eta()==constituents[j].eta()) QG=Qparticles[k].user_index();
       }
      }
      // cout<<QG<<endl;
      Jet *jet = new Jet(Qjets[i].perp(),Qjets[i].phi(),Qjets[i].eta(),QG);
      for (unsigned j = 0; j < constituents.size(); j++) {
       jet->add_incone(constituents[j].perp(),constituents[j].eta(),constituents[j].phi(),constituents[j].perp());
      }
      current.add_Qjet(jet);
     }
     for(unsigned i = 0; i < jets.size(); i++){
      vector<PseudoJet> constituents = jets[i].constituents();
      int QG=0;
      for (unsigned j = 0; j < constituents.size(); j++) {
       for(int k=0;k<Qparticles.size();k++){
        if(particles[k].perp()==constituents[j].perp() && particles[k].eta()==constituents[j].eta()) QG=particles[k].user_index();
       }
      }
      Jet *jet = new Jet(jets[i].perp(),jets[i].phi(),jets[i].eta(),QG);
      for (unsigned j = 0; j < constituents.size(); j++) {
       jet->add_incone(constituents[j].perp(),constituents[j].eta(),constituents[j].phi(),constituents[j].perp());
      }
      current.add_jet(jet);
     }
	   events.push_back(current);
     current.reset();
     particles.clear();
	   Qparticles.clear();
    }
   }
   input_file.close();
 } 
};
