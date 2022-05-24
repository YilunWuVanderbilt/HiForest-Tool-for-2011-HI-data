#include <memory>
#include <TMath.h>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Ref.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "math.h"

//classes to extract PFJet information
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrectionUncertainty.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/JetReco/interface/Jet.h"
//#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
//#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TTree.h"
#include "TFile.h"
#include<vector>

#include "TRandom3.h"


//
// class declaration
//

using namespace std;

class JetAnalyzer : public edm::EDAnalyzer {
public:
  explicit JetAnalyzer(const edm::ParameterSet&);
  ~JetAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
//  std::vector<float> factorLookup(float eta);
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);  
      
  //declare the input tag for PFJetCollection
  edm::InputTag jetInput;
  
  // ----------member data ---------------------------    
  
  bool isData;

  int numjet; //number of jets in the event
  TTree *mtree;
  std::vector<float> jet_e;
  std::vector<float> jet_pt;
  std::vector<float> jet_px;
  std::vector<float> jet_py;
  std::vector<float> jet_pz;
  std::vector<float> jet_eta;
  std::vector<float> jet_phi;
  std::vector<float> jet_ch;
  std::vector<float> jet_mass;
  std::vector<double> jet_btag;
  std::vector<float> corr_jet_pt;
  std::vector<float> corr_jet_ptUp;
  std::vector<float> corr_jet_ptDown;
  std::vector<float> corr_jet_ptSmearUp;
  std::vector<float> corr_jet_ptSmearDown;
  float btagWeight;
  float btagWeightUp;
  float btagWeightDn;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

JetAnalyzer::JetAnalyzer(const edm::ParameterSet& iConfig)
{
//now do what ever initialization is needed
  jetInput = iConfig.getParameter<edm::InputTag>("InputCollection");
  edm::Service<TFileService> fs;
  mtree = fs->make<TTree>("Events", "Events");
  
  isData = iConfig.getParameter<bool>("isData");

	
  mtree->Branch("numberjet",&numjet);
  mtree->GetBranch("numberjet")->SetTitle("Number of Jets");
  mtree->Branch("jet_e",&jet_e);
  mtree->GetBranch("jet_e")->SetTitle("Uncorrected Jet Energy");
  mtree->Branch("jet_pt",&jet_pt);
  mtree->GetBranch("jet_pt")->SetTitle("Uncorrected Transverse Jet Momentum");
  mtree->Branch("jet_px",&jet_px);
  mtree->GetBranch("jet_px")->SetTitle("X-Component of Jet Momentum");
  mtree->Branch("jet_py",&jet_py); 
  mtree->GetBranch("jet_py")->SetTitle("Y-Component of Jet Momentum");
  mtree->Branch("jet_pz",&jet_pz);
  mtree->GetBranch("jet_pz")->SetTitle("Z-Component of Jet Momentum");
  mtree->Branch("jet_eta",&jet_eta);
  mtree->GetBranch("jet_eta")->SetTitle("Jet Eta");
  mtree->Branch("jet_phi",&jet_phi);
  mtree->GetBranch("jet_phi")->SetTitle("Jet Phi");
  mtree->Branch("jet_ch",&jet_ch);
  mtree->GetBranch("jet_ch")->SetTitle("Jet Charge");
  mtree->Branch("jet_mass",&jet_mass);
  mtree->GetBranch("jet_mass")->SetTitle("Jet Mass");
  mtree->Branch("jet_btag",&jet_btag);
  mtree->GetBranch("jet_btag")->SetTitle("Jet Btagging Discriminant (CSV)");
  mtree->Branch("corr_jet_pt",&corr_jet_pt);
  mtree->GetBranch("corr_jet_pt")->SetTitle("Corrected Transverse Jet Momentum");
  mtree->Branch("corr_jet_ptUp",&corr_jet_ptUp);
  mtree->GetBranch("corr_jet_ptUp")->SetTitle("Corrected Transverse Jet Momentum (JEC Shifted Up)");
  mtree->Branch("corr_jet_ptDown",&corr_jet_ptDown);
  mtree->GetBranch("corr_jet_ptDown")->SetTitle("Corrected Transverse Jet Momentum (JEC Shifted Down)");
  mtree->Branch("corr_jet_ptSmearUp",&corr_jet_ptSmearUp);
  mtree->GetBranch("corr_jet_ptSmearUp")->SetTitle("Corrected Transverse Jet Momentum (JER Shifted Up)");
  mtree->Branch("corr_jet_ptSmearDown",&corr_jet_ptSmearDown);	
  mtree->GetBranch("corr_jet_ptSmearDown")->SetTitle("Corrected Transverse Jet Momentum (JER Shifted Down)");
  mtree->Branch("btag_Weight", &btagWeight);
  mtree->GetBranch("btag_Weight")->SetTitle("B-Tag event weight");
  mtree->Branch("btag_WeightUp", &btagWeightUp);
  mtree->GetBranch("btag_WeightUp")->SetTitle("B-Tag Up event weight");
  mtree->Branch("btag_WeightDn", &btagWeightDn);
  mtree->GetBranch("btag_WeightDn")->SetTitle("B-Tag Down event weight");
}

JetAnalyzer::~JetAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}
       
void
JetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  using namespace std;

  Handle<reco::CaloJetCollection> myjets;
  iEvent.getByLabel(jetInput, myjets);
  Handle<reco::JetTagCollection> btags;
  iEvent.getByLabel(InputTag("combinedSecondaryVertexBJetTags"), btags);
  Handle<double> rhoHandle;
  iEvent.getByLabel(InputTag("fixedGridRhoAll"), rhoHandle);
  Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(InputTag("offlinePrimaryVertices"), vertices);
  Handle<reco::JetFlavourMatchingCollection> injets;
  if (isData){
    iEvent.getByLabel(InputTag("InputCollection"), injets);
  }

  numjet = 0;
  jet_e.clear();
  jet_pt.clear();
  jet_px.clear();
  jet_py.clear();
  jet_pz.clear();
  jet_eta.clear();
  jet_phi.clear();
  jet_ch.clear();
  jet_mass.clear();
  jet_btag.clear();
  corr_jet_pt.clear();
  corr_jet_ptUp.clear();
  corr_jet_ptDown.clear();
  corr_jet_ptSmearUp.clear();
  corr_jet_ptSmearDown.clear();
  
  if (myjets.isValid()){
    double corr, corrUp, corrDown;
    float ptscale, ptscale_down, ptscale_up, res;
    double MC = 1;
    btagWeight = 1;
    btagWeightUp = 1;
    btagWeightDn = 1;
    int min_pt = 0;
    cout << "Integer printf(\"LoopTest1\\n\") = "<< MC << endl;
    for (reco::CaloJetCollection::const_iterator itjet=myjets->begin(); itjet!=myjets->end(); ++itjet){
      reco::Candidate::LorentzVector uncorrJet = itjet->p4();
      
      corr = 1;
      cout << "Integer printf(\"LoopTest2\\n\") = "<< corr << endl;
      
      corrUp = 1.0;
      corrDown = 1.0;

      ptscale = 1;
      ptscale_down = 1;
      ptscale_up = 1;
      res = 1;
      
      if (ptscale*corr*uncorrJet.pt() >= min_pt){
	
	jet_e.push_back(itjet->energy());
	jet_pt.push_back(itjet->pt());
	jet_px.push_back(itjet->px());
	jet_py.push_back(itjet->py());
	jet_pz.push_back(itjet->pz());
	jet_eta.push_back(itjet->eta());
	jet_phi.push_back(itjet->phi());
	jet_ch.push_back(itjet->charge());
	jet_mass.push_back(itjet->mass());
	if(btags.isValid() && (itjet - myjets->begin()) < btags->size()) {
	  jet_btag.push_back(btags->operator[](itjet - myjets->begin()).second);
	}
	else jet_btag.push_back(-999);
	corr_jet_pt.push_back(ptscale*corr*uncorrJet.pt());
	corr_jet_ptUp.push_back(ptscale*corrUp*uncorrJet.pt());
	corr_jet_ptDown.push_back(ptscale*corrDown*uncorrJet.pt());
	corr_jet_ptSmearUp.push_back(ptscale_up*corr*uncorrJet.pt());
	corr_jet_ptSmearDown.push_back(ptscale_down*corr*uncorrJet.pt());
	
	
	++numjet;
      }
    }
    btagWeight = (btagWeight/MC);
    btagWeightUp = (btagWeightUp/MC);
    btagWeightDn = (btagWeightDn/MC);   
  }
  
  mtree->Fill();
  return;
  
}

// ------------ method called once each job just before starting event loop  ------------
void
JetAnalyzer::beginJob()
{}

// ------------ method called once each job just after ending the event loop  ------------
void
JetAnalyzer::endJob()
{}

// ------------ method called when starting to processes a run  ------------
void
JetAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a run  ------------
void
JetAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{}
// ------------ method called when starting to processes a luminosity block  ------------
void
JetAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a luminosity block  ------------
void
JetAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
 //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetAnalyzer);
