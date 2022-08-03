// -*- C++ -*-
//
// Package:    EleTrkAna/ProducerTest
// Class:      ProducerTest
// 
/**\class ProducerTest ProducerTest.cc EleTrkAna/ProducerTest/plugins/ProducerTest.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Chenguang Zhang
//         Created:  Thu, 21 Jul 2022 17:16:49 GMT
//
//


// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "TrackingTools/GsfTools/interface/MultiGaussianStateTransform.h"
#include "TrackingTools/GsfTools/interface/GaussianSumUtilities1D.h"
#include "TrackingTools/GsfTools/interface/MultiGaussianState1D.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateMode.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "TrackingTools/GsfTracking/interface/GsfConstraintAtVertex.h"

#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronUtilities.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
//
// class declaration
//

class ProducerTest : public edm::stream::EDProducer<> {
   public:
      explicit ProducerTest(const edm::ParameterSet&);
      ~ProducerTest();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;
      
      void ParamCombine( const float &p, const float &perr, const float &e, const float &eerr, const int &elc, float &f, float &ferr );
      
      void gsfTrkMode_ParamCombine( pat::Electron &ele );
      void momatPCA_ParamCombine( pat::Electron &ele, const reco::BeamSpot &bs );
      void momatBS_ParamCombine( pat::Electron &ele, const reco::BeamSpot &bs, const edm::EventSetup& setup );

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;
      edm::EDGetTokenT<std::vector<reco::GsfElectron> > elecSrc_;
      //edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;
      //edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeometryToken_;
      edm::ESHandle<MagneticField> theMGField;
      edm::ESHandle<TrackerGeometry> theTrackerGeometry;
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
ProducerTest::ProducerTest(const edm::ParameterSet& iConfig):
  beamSpotSrc_(consumes<reco::BeamSpot>(iConfig.getUntrackedParameter<edm::InputTag>("beamSpotSrc"))),
  elecSrc_(consumes<std::vector<reco::GsfElectron> >(iConfig.getUntrackedParameter<edm::InputTag>("elecSrc")))
  //magneticFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
  //trackerGeometryToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>())
{
	produces<std::vector<pat::Electron> >();
}


ProducerTest::~ProducerTest()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
ProducerTest::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   
   edm::Handle<reco::BeamSpot> beamSpot;
   iEvent.getByToken(beamSpotSrc_,beamSpot);
   const reco::BeamSpot bs = *beamSpot;
   
   //input electrons
   edm::Handle<std::vector<reco::GsfElectron> > elesCollection;
   iEvent.getByToken(elecSrc_ , elesCollection);

   iSetup.get<IdealMagneticFieldRecord>().get(theMGField);
   iSetup.get<TrackerDigiGeometryRecord>().get(theTrackerGeometry);
   //auto const& theMGField = iSetup.getData(magneticFieldToken_);
   //auto const& theTrackerGeometry = iSetup.getData(trackerGeometryToken_);

   // output eles
   std::vector<pat::Electron> *theEles  = new std::vector<pat::Electron>;
   theEles->clear();

   for ( unsigned int i = 0; i < elesCollection->size(); i++ ){

	   pat::Electron e = elesCollection->at(i);
	   
	   ProducerTest::gsfTrkMode_ParamCombine( e );

	   momatPCA_ParamCombine( e, bs );
	   
	   momatBS_ParamCombine( e, bs, iSetup );

	   theEles->push_back( e );

   }

   std::unique_ptr<std::vector<pat::Electron> > ptr( theEles );
   iEvent.put( std::move( ptr ) );

}


void
ProducerTest::gsfTrkMode_ParamCombine( pat::Electron &ele ){
	
	float finalE = ele.correctedEcalEnergy();
	float finalEerr = ele.ecalEnergyError();

	const int elClass = ele.classification();

	auto const gsfTrk = ele.gsfTrack();
	const float trkP = gsfTrk->pMode();
	const float trkPErr = std::abs(gsfTrk->qoverpModeError())*trkP*trkP;

	const float EcalEnergy = ele.correctedEcalEnergy();
	const float EcalEnergyErr = ele.ecalEnergyError();

	ParamCombine( trkP, trkPErr, EcalEnergy, EcalEnergyErr, elClass, finalE, finalEerr );
	
	ele.addUserFloat( "elClass", elClass );
	ele.addUserFloat( "trkModeParamCombinedEnergy", finalE );
	ele.addUserFloat( "trkModeParamCombinedEnergy", finalEerr );
	ele.addUserFloat( "trkPMode", trkP );
	ele.addUserFloat( "trkPModeErr", trkPErr );

}

void ProducerTest::momatPCA_ParamCombine( pat::Electron &ele, const reco::BeamSpot &bs ){

	float finalE = ele.correctedEcalEnergy();
	float finalEerr = ele.ecalEnergyError();
	
	const int elClass = ele.classification();
	
	const float EcalEnergy = ele.correctedEcalEnergy();
	const float EcalEnergyErr = ele.ecalEnergyError();
	
	MultiTrajectoryStateTransform mtsTransform( theTrackerGeometry.product(), theMGField.product() );
	auto const gsfTrk = ele.gsfTrack();	
	TrajectoryStateOnSurface innTSOS = mtsTransform.innerStateOnSurface(*gsfTrk);
	if ( !innTSOS.isValid() )return;
	
	GlobalPoint bsPos;
	ele_convert( bs.position(), bsPos );
	TrajectoryStateOnSurface pcaTSOS = mtsTransform.extrapolatedState(innTSOS, bsPos);
	if ( !pcaTSOS.isValid() )
		pcaTSOS = innTSOS;

	GlobalVector pcaMom;
	math::XYZPointF momatPAC;
	multiTrajectoryStateMode::momentumFromModeCartesian( pcaTSOS, pcaMom );
	ele_convert( pcaMom, momatPAC );
	const float trkP = momatPAC.R();

	MultiGaussianState1D qpState(MultiGaussianStateTransform::multiState1D(pcaTSOS, 0));
	GaussianSumUtilities1D qpUtils(qpState);
	const float trkPErr = trkP * trkP * sqrt(qpUtils.mode().variance());

	ParamCombine( trkP, trkPErr, EcalEnergy, EcalEnergyErr, elClass, finalE, finalEerr );
	

	if ( !ele.hasUserFloat( "elClass" ) ){
		ele.addUserFloat( "elClass", elClass );
	}
	ele.addUserFloat( "momatPCAParamCombinedEnergy", finalE );
	ele.addUserFloat( "momatPCAParamCombinedEnergyErr", finalEerr );
	ele.addUserFloat( "momatPCAErr", trkPErr );
	ele.addUserFloat( "momatPCA", trkP );

}

void ProducerTest::momatBS_ParamCombine( pat::Electron &ele, const reco::BeamSpot &bs, const edm::EventSetup& setup ){

	float finalE = ele.correctedEcalEnergy();
	float finalEerr = ele.ecalEnergyError();
	
	const int elClass = ele.classification();
	
	const float EcalEnergy = ele.correctedEcalEnergy();
	const float EcalEnergyErr = ele.ecalEnergyError();

	MultiTrajectoryStateTransform mtsTransform( theTrackerGeometry.product(), theMGField.product() );
	GsfConstraintAtVertex constraintAtBS( setup );
	auto const gsfTrk = ele.gsfTrack();
		
	TrajectoryStateOnSurface innTSOS = mtsTransform.innerStateOnSurface(*gsfTrk);
	if ( !innTSOS.isValid() )
		return;
		
	TrajectoryStateOnSurface bsTSOS = constraintAtBS.constrainAtBeamSpot( *gsfTrk, bs );

	GlobalVector bsMom;
	math::XYZPointF momatBS;
	multiTrajectoryStateMode::momentumFromModeCartesian( bsTSOS, bsMom );
	ele_convert( bsMom, momatBS );
	const float trkP = momatBS.R();
	
	MultiGaussianState1D qpState(MultiGaussianStateTransform::multiState1D( bsTSOS, 0 ));
	GaussianSumUtilities1D qpUtils(qpState);
	const float trkPErr = trkP * trkP * sqrt(qpUtils.mode().variance());

	ParamCombine( trkP, trkPErr, EcalEnergy, EcalEnergyErr, elClass, finalE, finalEerr );

	if ( !ele.hasUserFloat( "elClass" ) ){
		ele.addUserFloat( "elClass", elClass );
	}
	ele.addUserFloat( "momatBSParamCombinedEnergy", finalE );
	ele.addUserFloat( "momatBSParamCombinedEnergyErr", finalEerr );
	ele.addUserFloat( "momatBSErr", trkPErr );
	ele.addUserFloat( "momatBS", trkP );

}


void ProducerTest::ParamCombine( const float &trkP, const float &trkPErr, const float &ecalE, const float &ecalEerr, const int &elClass, float &finalE, float &finalEerr ){

	if ( trkPErr/trkP > 0.5 && ecalEerr/ecalE <=0.5 ){

		finalE = ecalE; finalEerr = ecalEerr;
	
	}else if ( trkPErr/trkP <= 0.5 && ecalEerr/ecalE > 0.5 ){

		finalE = trkP; finalEerr = trkPErr;

	}else if ( trkPErr/trkP > 0.5 && ecalEerr/ecalE > 0.5 ){

		if ( trkPErr/trkP < ecalEerr/ecalE ){

			finalE = trkP; finalEerr = trkPErr;

		} else {

			finalE = ecalE; finalEerr = ecalEerr;
		}
	
	} else {

		float eOverP = ecalE/trkP;
		float errEOverP = sqrt( (ecalEerr/trkP)*(ecalEerr/trkP) + (ecalE*trkPErr/trkP/trkP)*(ecalE*trkPErr/trkP/trkP) );
		bool eleIsNotInCombination = false;
		if ( (eOverP > 1.0 + 2.5 * errEOverP) || (eOverP < 1.0 - 2.5 * errEOverP) || (eOverP < 0.8) || (eOverP > 1.3) ){
			eleIsNotInCombination = true;
		}
		if ( eleIsNotInCombination ) {
			if ( eOverP > 1.0 ) {
				finalE = ecalE; finalEerr = ecalEerr;
			} else {
				if ( elClass == reco::GsfElectron::GOLDEN ){
					finalE = ecalE; finalEerr = ecalEerr;
				}
				if ( elClass == reco::GsfElectron::BIGBREM ){
					if ( ecalE < 36.0 ){
						finalE = trkP; finalEerr = trkPErr;
					} else {
						finalE = ecalE; finalEerr = ecalEerr;
					}
				}
				if ( elClass == reco::GsfElectron::BADTRACK ){
					finalE = ecalE; finalEerr = ecalEerr;
				}
				if ( elClass == reco::GsfElectron::SHOWERING ){
					if ( ecalE < 30.0 ){
						finalE = trkP; finalEerr = trkPErr;
					} else {
						finalE = ecalE; finalEerr = ecalEerr;
					}
				}
				if ( elClass == reco::GsfElectron::GAP ){
					if ( ecalE < 60.0 ){
						finalE = trkP; finalEerr = trkPErr;
					} else {
						finalE = ecalE; finalEerr = ecalEerr;
					}
				}
			}
		} else {
			finalE = ( ecalE/ecalEerr/ecalEerr + trkP/trkPErr/trkPErr ) / ( 1/ecalEerr/ecalEerr + 1/trkPErr/trkPErr );
			finalEerr = sqrt( 1/(1/ecalEerr/ecalEerr + 1/trkPErr/trkPErr) );
		}
	}
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
ProducerTest::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
ProducerTest::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
ProducerTest::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
ProducerTest::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
ProducerTest::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
ProducerTest::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ProducerTest::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ProducerTest);
