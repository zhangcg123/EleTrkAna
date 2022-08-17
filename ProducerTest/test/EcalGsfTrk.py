import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

import glob
mylist = glob.glob('/eos/user/c/chenguan/CondorOutputs/TEST/out_4452341*root')
for i in range(len(mylist)):
	mylist[i] = 'file:'+mylist[i]

process = cms.Process("ProducerTest")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag='106X_upgrade2018_realistic_v16_L1v1'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(*mylist),
    #fileNames = cms.untracked.vstring('file:/eos/user/c/chenguan/CondorOutputs/TEST/out_3865280_39.root'),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
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
    fileName = cms.untracked.string('ecalgsftrk4.root'),
    outputCommands = cms.untracked.vstring(
    	'drop *',
      	'keep patElectrons_ProducerTest_*_*',
	'keep recoVertexs_offlinePrimaryVertices_*_*',
	#'keep recoBeamSpot_offlineBeamSpot_*_*',
	#'keep *_calibratedElectrons_*_*',
	#'keep *_egmGsfElectronIDs_*_*',
	#'keep *_gedGsfElectrons_*_*',
    	)
)

process.p = cms.Path( process.egammaPostRecoSeq * process.ProducerTest )
process.e = cms.EndPath(process.out)
