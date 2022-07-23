import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
					'file:/eos/user/c/chenguan/Tuples/ForEleTrck/APieceOfDYMINIAODSIM/001C8DDF-599C-5E45-BF2C-76F887C9ADE9.root'
					)
                            )

process.demo = cms.EDAnalyzer('DemoAnalyzer',
	vertexSrc    = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"),
	beamSpotSrc  = cms.untracked.InputTag("offlineBeamSpot"),
	electronSrc  = cms.untracked.InputTag("slimmedElectrons"),
                              )

process.p = cms.Path(process.demo)
