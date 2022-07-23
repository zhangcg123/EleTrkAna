// -*- C++ -*-
//
// Package:    EleTrkAna/DemoAnalyzer
// Class:      DemoAnalyzer
//
/**\class DemoAnalyzer DemoAnalyzer.cc EleTrkAna/DemoAnalyzer/plugins/DemoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Chenguang Zhang
//         Created:  Thu, 21 Jul 2022 09:22:11 GMT
//
//


// system include files
#include <memory>

#include "TTree.h"
#include "TLorentzVector.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class DemoAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit DemoAnalyzer(const edm::ParameterSet&);
      ~DemoAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      void bookPassedEventTree(TString treeName, TTree *tree);
      void assign_kinematics_to_leptons( std::vector<pat::Electron> &eles );


      TTree *passedEventsTree_All;


      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
      edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;
      edm::EDGetTokenT<edm::View<pat::Electron> > elecSrc_;
      
      edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;
      edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeometryToken_; 

      // ----------event variables -----------------------
      ULong64_t Run, Event, LumiSect;

      float ele1_correctedEcalEnergy, ele2_correctedEcalEnergy;
      float ele1_ecalTrkEnergyPreCorr, ele2_ecalTrkEnergyPreCorr;
      
      float pt1_ecaltrk, pt2_ecaltrk, eta1_ecaltrk, eta2_ecaltrk, phi1_ecaltrk, phi2_ecaltrk, m1_ecaltrk, m2_ecaltrk, massz_ecaltrk;
      float pt1_ecal, pt2_ecal, eta1_ecal, eta2_ecal, phi1_ecal, phi2_ecal, m1_ecal, m2_ecal, massz_ecal;
      float pt1_trkmean, pt2_trkmean, eta1_trkmean, eta2_trkmean, phi1_trkmean, phi2_trkmean, m1_trkmean, m2_trkmean, massz_trkmean;
      float pt1_trkmode, pt2_trkmode, eta1_trkmode, eta2_trkmode, phi1_trkmode, phi2_trkmode, m1_trkmode, m2_trkmode, massz_trkmode;
      float pt1_trkatpca, pt2_trkatpca, eta1_trkatpca, eta2_trkatpca, phi1_trkatpca, phi2_trkatpca, m1_trkatpca, m2_trkatpca, massz_trkatpca;
      float pt1_trkatbs, pt2_trkatbs, eta1_trkatbs, eta2_trkatbs, phi1_trkatbs, phi2_trkatbs, m1_trkatbs, m2_trkatbs, massz_trkatbs;
      
      int theVertex;
      bool passedZSelection;
      float corr1, corr2;
      TLorentzVector lep1, lep2;
      std::vector<pat::Electron> tightelectrons;

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
DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig)
 :
  vertexSrc_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc"))),
  beamSpotSrc_(consumes<reco::BeamSpot>(iConfig.getUntrackedParameter<edm::InputTag>("beamSpotSrc"))),
  elecSrc_(consumes<edm::View<pat::Electron> >(iConfig.getUntrackedParameter<edm::InputTag>("electronSrc"))),
  magneticFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
  trackerGeometryToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>())
{
   //now do what ever initialization is needed

}


DemoAnalyzer::~DemoAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace pat;

   // Get collections
   edm::Handle<reco::VertexCollection> vertex;
   iEvent.getByToken(vertexSrc_,vertex);

   edm::Handle<edm::View<pat::Electron> > electrons;
   iEvent.getByToken(elecSrc_,electrons);

   edm::Handle<reco::BeamSpot> beamSpot;
   iEvent.getByToken(beamSpotSrc_,beamSpot);
   //const reco::BeamSpot BS = *beamSpot;

   edm::ESHandle<MagneticField> theMGField = iSetup.getHandle(magneticFieldToken_);
   edm::ESHandle<TrackerGeometry> theTracerGeometry = iSetup.getHandle(trackerGeometryToken_);

   // initialize variables
   Run = -999; 
   Event = -999;
   LumiSect = -999;

   ele1_correctedEcalEnergy = -999.9; ele2_correctedEcalEnergy = -999.9;
   ele1_ecalTrkEnergyPreCorr = -999.9; ele2_ecalTrkEnergyPreCorr = -999.9;
   
   pt1_ecaltrk = -999.9; pt2_ecaltrk = -999.9; eta1_ecaltrk = -999.9; eta2_ecaltrk = -999.9; phi1_ecaltrk = -999.9; phi2_ecaltrk = -999.9; 
   m1_ecaltrk = -999.9; m2_ecaltrk = -999.9; massz_ecaltrk = -999.9;
   pt1_ecal = -999.9; pt1_ecal = -999.9; eta1_ecal = -999.9; eta2_ecal = -999.9; phi1_ecal = -999.9; phi2_ecal = -999.9;
   m1_ecal = -999.9; m2_ecal = -999.9; massz_ecal = 999.9;
   pt1_trkmean = -999.9; pt2_trkmean = -999.9; eta1_trkmean = -999.9; eta2_trkmean = -999.9; phi1_trkmean = -999.9; phi2_trkmean = -999.9; 
   m1_trkmean = -999.9; m2_trkmean = -999.9; massz_trkmean = -999.9;
   pt1_trkmode = -999.9; pt2_trkmode = -999.9; eta1_trkmode = -999.9; eta2_trkmode = -999.9; phi1_trkmode = -999.9; phi2_trkmode = -999.9; 
   m1_trkmode = -999.9; m2_trkmode = -999.9; massz_trkmode = -999.9;
   pt1_trkatpca = -999.9; pt2_trkatpca = -999.9; eta1_trkatpca = -999.9; eta2_trkatpca = -999.9; phi1_trkatpca = -999.9; phi2_trkatpca = -999.9; 
   m1_trkatpca = -999.9; m2_trkatpca = -999.9; massz_trkatpca = -999.9;
   pt1_trkatbs = -999.9; pt2_trkatbs = -999.9; eta1_trkatbs = -999.9; eta2_trkatbs = -999.9; phi1_trkatbs = -999.9; phi2_trkatbs = -999.9; 
   m1_trkatbs = -999.9; m2_trkatbs = -999.9; massz_trkatbs = -999.9;

   passedZSelection = false;
   corr1 = 1.0; corr2 = 1.0;
   lep1.SetPtEtaPhiM(0,0,0,0); lep2.SetPtEtaPhiM(0,0,0,0);
   tightelectrons.clear();
   theVertex = -1;

   for ( unsigned int i = 0; i < vertex->size(); i++){
	   const reco::Vertex *PV = &(vertex->at(i));
	   if (PV->isFake()) continue;
	   if (PV->ndof()<=4 || PV->position().Rho()>2.0 || fabs(PV->position().Z())>24.0) continue;
	   theVertex=(int)i; break;
   }
   
   if ( theVertex >=0 ){

   	for ( edm::View<pat::Electron>::const_iterator ele=electrons->begin(); ele != electrons->end(); ++ele){

		if ( ele->isElectronIDAvailable("cutBasedElectronID-Fall17-94X-V2-tight") || ele->electronID("cutBasedElectronID-Fall17-94X-V2-tight") ){

			tightelectrons.push_back( *ele );
		}
	}
	if ( tightelectrons.size() == 2 ){
	
		if ( tightelectrons[0].charge() * tightelectrons[1].charge() < 0 ){
	
			passedZSelection = true;
		}
	}

	if ( passedZSelection ){
		
		//sort as pt
		if ( tightelectrons[0].pt() < tightelectrons[1].pt() ){
			pat::Electron tmp = tightelectrons[0];
			tightelectrons.push_back( tmp );
			tightelectrons.erase( tightelectrons.begin() );
		}

		// Beamspot constrain and E-p combination
		//beamspot_constrain_and_epcombination( tightelectrons, theMGField, theTracerGeometry );

		// assign lepton kinematics
		assign_kinematics_to_leptons( tightelectrons );

		// add lepton stuff and probably ele constraint
		passedEventsTree_All->Fill();
		
	}// if ( passedZSelection )

   }// if ( theVertex >=0 )

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void
DemoAnalyzer::beginJob()
{
	using namespace edm;
	using namespace std;
	using namespace pat;

	bookPassedEventTree("passedEvents", passedEventsTree_All);
}

// ------------ method called once each job just after ending the event loop  ------------
void
DemoAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

void
DemoAnalyzer::assign_kinematics_to_leptons( std::vector<pat::Electron> &eles) {
		
		ele1_ecalTrkEnergyPreCorr = eles[0].userFloat("ecalTrkEnergyPreCorr");
		ele2_ecalTrkEnergyPreCorr = eles[1].userFloat("ecalTrkEnergyPreCorr");
		ele1_correctedEcalEnergy = eles[0].correctedEcalEnergy();
		ele2_correctedEcalEnergy = eles[1].correctedEcalEnergy();
		
		// lep.Energy = ecalTrkEnergyPreCorr
		pt1_ecaltrk = eles[0].pt();
		pt2_ecaltrk = eles[1].pt();
		eta1_ecaltrk = eles[0].eta();
		eta2_ecaltrk = eles[1].eta();
		phi1_ecaltrk = eles[0].phi();
		phi2_ecaltrk = eles[1].phi();
		m1_ecaltrk = 0.00051;//eles[0].mass();
		m2_ecaltrk = 0.00051;//eles[1].mass();
		lep1.SetPtEtaPhiM( pt1_ecaltrk, eta1_ecaltrk, phi1_ecaltrk, m1_ecaltrk );
		lep2.SetPtEtaPhiM( pt2_ecaltrk, eta2_ecaltrk, phi2_ecaltrk, m2_ecaltrk );
		massz_ecaltrk = ( lep1 + lep2 ).M();
		
		// lep.Energy = correctedEcalEnergy
		corr1 = eles[0].correctedEcalEnergy()/eles[0].userFloat("ecalTrkEnergyPreCorr");
		corr2 = eles[1].correctedEcalEnergy()/eles[1].userFloat("ecalTrkEnergyPreCorr");
		pt1_ecal = pt1_ecaltrk * corr1;
		pt2_ecal = pt2_ecaltrk * corr2;
		eta1_ecal = eta1_ecaltrk;
		eta2_ecal = eta2_ecaltrk;
		phi1_ecal = phi1_ecaltrk;
		phi2_ecal = phi2_ecaltrk;
		m1_ecal = m1_ecaltrk;
		m2_ecal = m1_ecaltrk;
		lep1.SetPtEtaPhiM( pt1_ecal, eta1_ecal, phi1_ecal, m1_ecal );
		lep2.SetPtEtaPhiM( pt2_ecal, eta2_ecal, phi2_ecal, m2_ecal );
		massz_ecal = ( lep1 + lep2 ).M();

		// lep.Energy = ele.gsfTrack.pMode
		auto gsfTrk1 = eles[0].gsfTrack();
		auto gsfTrk2 = eles[1].gsfTrack();
		
		pt1_trkmode = gsfTrk1->ptMode();
		pt2_trkmode = gsfTrk2->ptMode();
		eta1_trkmode = gsfTrk1->etaMode();
		eta2_trkmode = gsfTrk2->etaMode();
		phi1_trkmode = gsfTrk1->phiMode();
		phi2_trkmode = gsfTrk2->phiMode();
		m1_trkmode = m1_ecaltrk;
		m2_trkmode = m2_ecaltrk;
		lep1.SetPtEtaPhiM( pt1_trkmode, eta1_trkmode, phi1_trkmode, m1_trkmode );
		lep2.SetPtEtaPhiM( pt2_trkmode, eta2_trkmode, phi2_trkmode, m2_trkmode );
		massz_trkmode = ( lep1 + lep2 ).M();

		// lep.Energy = ele.gsfTrack.p
		pt1_trkmean = gsfTrk1->pt();
		pt2_trkmean = gsfTrk2->pt();
		eta1_trkmean = gsfTrk1->eta();
		eta2_trkmean = gsfTrk2->eta();
		phi1_trkmean = gsfTrk1->phi();
		phi2_trkmean = gsfTrk2->phi();
		m1_trkmean = m1_trkmode;
		m2_trkmean = m2_trkmode;
		lep1.SetPtEtaPhiM( pt1_trkmean, eta1_trkmean, phi1_trkmean, m1_trkmean );
		lep2.SetPtEtaPhiM( pt2_trkmean, eta2_trkmean, phi2_trkmean, m2_trkmean );
		massz_trkmean = ( lep1 + lep2 ).M();

		// lep.Energy = lep.trackMomentumAtVtx.R
		pt1_trkatpca = eles[0].trackMomentumAtVtx().Rho();
		pt2_trkatpca = eles[1].trackMomentumAtVtx().Rho();
		eta1_trkatpca = eles[0].trackMomentumAtVtx().Eta();
		eta2_trkatpca = eles[1].trackMomentumAtVtx().Eta();
		phi1_trkatpca = eles[0].trackMomentumAtVtx().Phi();
		phi2_trkatpca = eles[1].trackMomentumAtVtx().Phi();
		m1_trkatpca = m1_trkmean;
		m2_trkatpca = m2_trkmean;
		lep1.SetPtEtaPhiM( pt1_trkatpca, eta1_trkatpca, phi1_trkatpca, m1_trkatpca );
		lep2.SetPtEtaPhiM( pt2_trkatpca, eta2_trkatpca, phi2_trkatpca, m2_trkatpca );
		massz_trkatpca = ( lep1 + lep2 ).M();

		// lep.Energy = lep.trackMomentumAtVtxWithConstraint.R
		pt1_trkatbs = eles[0].trackMomentumAtVtxWithConstraint().Rho();
		pt2_trkatbs = eles[1].trackMomentumAtVtxWithConstraint().Rho();
		eta1_trkatbs = eles[0].trackMomentumAtVtxWithConstraint().Eta();
		eta2_trkatbs = eles[1].trackMomentumAtVtxWithConstraint().Eta();
		phi1_trkatbs = eles[0].trackMomentumAtVtxWithConstraint().Phi();
		phi1_trkatbs = eles[1].trackMomentumAtVtxWithConstraint().Phi();
		m1_trkatbs = m1_trkmean;
		m2_trkatbs = m2_trkmean;
		lep1.SetPtEtaPhiM( pt1_trkatbs, eta1_trkatbs, phi1_trkatbs, m1_trkatbs );
		lep2.SetPtEtaPhiM( pt2_trkatbs, eta2_trkatbs, phi2_trkatbs, m2_trkatbs );
		massz_trkatbs = ( lep1 + lep2 ).M();

}

void
DemoAnalyzer::bookPassedEventTree(TString treeName, TTree *tree) {

	using namespace edm;
	using namespace pat;
	using namespace std;

	tree->Branch("Run",&Run,"Run/l");
	tree->Branch("Event",&Event,"Event/l");
	tree->Branch("LumiSect",&LumiSect,"LumiSect/l");
	tree->Branch("ele1_ecalTrkEnergyPreCorr",&ele1_ecalTrkEnergyPreCorr,"ele1_ecalTrkEnergyPreCorr/F");
	tree->Branch("ele2_ecalTrkEnergyPreCorr",&ele2_ecalTrkEnergyPreCorr,"ele2_ecalTrkEnergyPreCorr/F");
	tree->Branch("ele1_correctedEcalEnergy",&ele1_correctedEcalEnergy,"ele1_correctedEcalEnergy/F");
	tree->Branch("ele2_correctedEcalEnergy",&ele2_correctedEcalEnergy,"ele2_correctedEcalEnergy/F");
	tree->Branch("pt1_ecaltrk",&pt1_ecaltrk,"pt1_ecaltrk/F");
	tree->Branch("pt2_ecaltrk",&pt2_ecaltrk,"pt2_ecaltrk/F");
	tree->Branch("eta1_ecaltrk",&eta1_ecaltrk,"eta1_ecaltrk/F");
	tree->Branch("eta1_ecaltrk",&eta2_ecaltrk,"eta1_ecaltrk/F");
	tree->Branch("phi1_ecaltrk",&phi1_ecaltrk,"phi1_ecaltrk/F");
	tree->Branch("phi2_ecaltrk",&phi2_ecaltrk,"phi2_ecaltrk/F");
	tree->Branch("m1_ecaltrk",&m1_ecaltrk,"m1_ecaltrk/F");
	tree->Branch("m2_ecaltrk",&m2_ecaltrk,"m2_ecaltrk/F");
	tree->Branch("massz_ecaltrk",&massz_ecaltrk,"massz_ecaltrk/F");

	tree->Branch("pt1_ecal",&pt1_ecal,"pt1_ecal/F");
	tree->Branch("pt2_ecal",&pt2_ecal,"pt2_ecal/F");
	tree->Branch("eta1_ecal",&eta1_ecal,"eta1_ecal/F");
	tree->Branch("eta1_ecal",&eta2_ecal,"eta1_ecal/F");
	tree->Branch("phi1_ecal",&phi1_ecal,"phi1_ecal/F");
	tree->Branch("phi2_ecal",&phi2_ecal,"phi2_ecal/F");
	tree->Branch("m1_ecal",&m1_ecal,"m1_ecal/F");
	tree->Branch("m2_ecal",&m2_ecal,"m2_ecal/F");
	tree->Branch("massz_ecal",&massz_ecal,"massz_ecal/F");

	tree->Branch("pt1_trkmean",&pt1_trkmean,"pt1_trkmean/F");
	tree->Branch("pt2_trkmean",&pt2_trkmean,"pt2_trkmean/F");
	tree->Branch("eta1_trkmean",&eta1_trkmean,"eta1_trkmean/F");
	tree->Branch("eta1_trkmean",&eta2_trkmean,"eta1_trkmean/F");
	tree->Branch("phi1_trkmean",&phi1_trkmean,"phi1_trkmean/F");
	tree->Branch("phi2_trkmean",&phi2_trkmean,"phi2_trkmean/F");
	tree->Branch("m1_trkmean",&m1_trkmean,"m1_trkmean/F");
	tree->Branch("m2_trkmean",&m2_trkmean,"m2_trkmean/F");
	tree->Branch("massz_trkmean",&massz_trkmean,"massz_trkmean/F");

	tree->Branch("pt1_trkmode",&pt1_trkmode,"pt1_trkmode/F");
	tree->Branch("pt2_trkmode",&pt2_trkmode,"pt2_trkmode/F");
	tree->Branch("eta1_trkmode",&eta1_trkmode,"eta1_trkmode/F");
	tree->Branch("eta1_trkmode",&eta2_trkmode,"eta1_trkmode/F");
	tree->Branch("phi1_trkmode",&phi1_trkmode,"phi1_trkmode/F");
	tree->Branch("phi2_trkmode",&phi2_trkmode,"phi2_trkmode/F");
	tree->Branch("m1_trkmode",&m1_trkmode,"m1_trkmode/F");
	tree->Branch("m2_trkmode",&m2_trkmode,"m2_trkmode/F");
	tree->Branch("massz_trkmode",&massz_trkmode,"massz_trkmode/F");

	tree->Branch("pt1_trkatpca",&pt1_trkatpca,"pt1_trkatpca/F");
	tree->Branch("pt2_trkatpca",&pt2_trkatpca,"pt2_trkatpca/F");
	tree->Branch("eta1_trkatpca",&eta1_trkatpca,"eta1_trkatpca/F");
	tree->Branch("eta1_trkatpca",&eta2_trkatpca,"eta1_trkatpca/F");
	tree->Branch("phi1_trkatpca",&phi1_trkatpca,"phi1_trkatpca/F");
	tree->Branch("phi2_trkatpca",&phi2_trkatpca,"phi2_trkatpca/F");
	tree->Branch("m1_trkatpca",&m1_trkatpca,"m1_trkatpca/F");
	tree->Branch("m2_trkatpca",&m2_trkatpca,"m2_trkatpca/F");
	tree->Branch("massz_trkatpca",&massz_trkatpca,"masz_trkatpca/F");
	
	tree->Branch("pt1_trkatbs",&pt1_trkatbs,"pt1_trkatbs/F");
	tree->Branch("pt2_trkatbs",&pt2_trkatbs,"pt2_trkatbs/F");
	tree->Branch("eta1_trkatbs",&eta1_trkatbs,"eta1_trkatbs/F");
	tree->Branch("eta1_trkatbs",&eta2_trkatbs,"eta1_trkatbs/F");
	tree->Branch("phi1_trkatbs",&phi1_trkatbs,"phi1_trkatbs/F");
	tree->Branch("phi2_trkatbs",&phi2_trkatbs,"phi2_trkatbs/F");
	tree->Branch("m1_trkatbs",&m1_trkatbs,"m1_trkatbs/F");
	tree->Branch("m2_trkatbs",&m2_trkatbs,"m2_trkatbs/F");
	tree->Branch("massz_trkatbs",&massz_trkatbs,"massz_trkbs/F");
}
//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);
