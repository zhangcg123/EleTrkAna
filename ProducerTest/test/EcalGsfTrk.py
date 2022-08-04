import FWCore.ParameterSet.Config as cms

process = cms.Process("ProducerTest")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag='106X_upgrade2018_realistic_v16_L1v1'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       #'file:/eos/user/c/chenguan/Tuples/ForEleTrck/ParentFilesOfTheMINIAODSIM/13631887-2D4C-2F4A-8575-59FBD01348E7.root'
       'file:/eos/user/c/chenguan/CondorOutputs/tmp/out_3865187_0.root',
    )
)

process.ProducerTest = cms.EDProducer('ProducerTest',
	elecSrc = cms.untracked.InputTag('gedGsfElectrons'),
	elecIdSrc = cms.untracked.InputTag('egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-tight'),
	beamSpotSrc  = cms.untracked.InputTag("offlineBeamSpot"),
)

from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
			isMiniAOD=False,
			runEnergyCorrections=True,
			runVID=True,
			eleIDModules=['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff'],
			#phoIDModules=['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V2_cff'],
			era='2018-UL')

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('ecalgsftrk.root'),
    outputCommands = cms.untracked.vstring(
    	'drop *',
      	'keep patElectrons_ProducerTest_*_*',
    	'keep recoBeamSpot_offlineBeamSpot_*_*',
	'keep *_calibratedElectrons_*_*',
	#'keep *_egmGsfElectronIDs_*_*',
	#'keep *_gedGsfElectrons_*_*',
    	)
)

process.p = cms.Path( process.egammaPostRecoSeq * process.ProducerTest )
process.e = cms.EndPath(process.out)
