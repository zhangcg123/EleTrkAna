import FWCore.ParameterSet.Config as cms

process = cms.Process("ProducerTest")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
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
	beamSpotSrc  = cms.untracked.InputTag("offlineBeamSpot"),
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('ecalgsftrk.root'),
    outputCommands = cms.untracked.vstring(
      	"keep patElectrons_ProducerTest_*_*",
	"keep recoBeamSpot_offlineBeamSpot_*_*",
	)
)

process.p = cms.Path(process.ProducerTest)
process.e = cms.EndPath(process.out)
