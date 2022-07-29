#! /usr/bin/env python3

import ROOT
import sys
from DataFormats.FWLite import Events, Handle
import pandas
sys.path.append('/afs/cern.ch/work/c/chenguan/private/pycommontool/')
from FileSystemClass import *
from HTMLClass import *

dirs = DirTree()
dirs.mkrootdir('test_dummy')
path = dirs.root

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.parseArguments()

events = Events ("ecalgsftrk.root")
handle  = Handle ("std::vector<pat::Electron>")
label = ("ProducerTest")

massZ_standard = ROOT.TH1F('massz_standard','',100,70,110)
massZ_trkmode = ROOT.TH1F('massz_trkmode','',100,70,110)
massZ_trkmode_paramcombine = ROOT.TH1F('massz_trkmode_paramcombine','',100,70,110)
massZ_trkmode_replacepmodewithbsconstr_paramcombine = ROOT.TH1F('massZ_trkmode_replacepmodewithbsconstr_paramcombine','',100,70,110)
massZ_trkmode_replacepmodewithpca_paramcombine = ROOT.TH1F('massZ_trkmode_replacepmodewithpca_paramcombine','',100,70,110)

# loop over events
for i,event in enumerate(events):
    #if i > options.maxEvents: break
    event.getByLabel (label, handle)
    eles = handle.product()

    tighteles = []
    lep1 = ROOT.TLorentzVector()
    lep2 = ROOT.TLorentzVector()
    for e in eles:
	    #tight electron
	    if e.electronID('cutBasedElectronID-Fall17-94X-V2-tight') == 1:
		    tighteles.append( e )
	    
	    if len( tighteles ) == 2 and tighteles[0].charge() * tighteles[1].charge() < 0:

		    #build z
		    lep1.SetPxPyPzE( tighteles[0].px(), tighteles[0].py(), tighteles[0].pz(), tighteles[0].userFloat('ecalTrkEnergyPreCorr') )
		    lep2.SetPxPyPzE( tighteles[1].px(), tighteles[1].py(), tighteles[1].pz(), tighteles[1].userFloat('ecalTrkEnergyPreCorr') )
		    massz_standard = ( lep1 + lep2 ).M()
		    massZ_standard.Fill( massz_standard )

		    lep1.SetPxPyPzE( tighteles[0].gsfTrack().pxMode(), tighteles[0].gsfTrack().pyMode(), tighteles[0].gsfTrack().pzMode(), tighteles[0].gsfTrack().pMode() )
		    lep2.SetPxPyPzE( tighteles[1].gsfTrack().pxMode(), tighteles[1].gsfTrack().pyMode(), tighteles[1].gsfTrack().pzMode(), tighteles[1].gsfTrack().pMode() )
		    massz_trkmode_pxpypze = ( lep1 + lep2 ).M()
		    #same as above
		    #lep1.SetPtEtaPhiM( tighteles[0].gsfTrack().ptMode(), tighteles[0].gsfTrack().etaMode(), tighteles[0].gsfTrack().phiMode(), 0.00051 )
		    #lep2.SetPtEtaPhiM( tighteles[1].gsfTrack().ptMode(), tighteles[1].gsfTrack().etaMode(), tighteles[1].gsfTrack().phiMode(), 0.00051 )
		    #massz_trkmode_ptetaphim = ( lep1 + lep2 ).M()

		    corr1 = tighteles[0].userFloat('trkModeParamCombinedEnergy')/tighteles[0].gsfTrack().pMode()
		    corr2 = tighteles[1].userFloat('trkModeParamCombinedEnergy')/tighteles[1].gsfTrack().pMode()
		    lep1.SetPxPyPzE( corr1*tighteles[0].gsfTrack().pxMode(), corr1*tighteles[0].gsfTrack().pyMode(), corr1*tighteles[0].gsfTrack().pzMode(), tighteles[0].userFloat('trkModeParamCombinedEnergy') )
		    lep2.SetPxPyPzE( corr2*tighteles[1].gsfTrack().pxMode(), corr2*tighteles[1].gsfTrack().pyMode(), corr2*tighteles[1].gsfTrack().pzMode(), tighteles[1].userFloat('trkModeParamCombinedEnergy') ) 
		    massz_trkmode_paramcombine = ( lep1 + lep2 ).M()

		    massZ_trkmode.Fill( massz_trkmode_pxpypze )
		    massZ_trkmode_paramcombine.Fill( massz_trkmode_paramcombine )

		    corr1 = tighteles[0].userFloat('trkModeParamCombinedEnergy_1')/tighteles[0].gsfTrack().pMode()
		    corr2 = tighteles[1].userFloat('trkModeParamCombinedEnergy_1')/tighteles[1].gsfTrack().pMode()
		    lep1.SetPxPyPzE( corr1*tighteles[0].gsfTrack().pxMode(), corr1*tighteles[0].gsfTrack().pyMode(), corr1*tighteles[0].gsfTrack().pzMode(), tighteles[0].userFloat('trkModeParamCombinedEnergy_1') )
		    lep2.SetPxPyPzE( corr2*tighteles[1].gsfTrack().pxMode(), corr2*tighteles[1].gsfTrack().pyMode(), corr2*tighteles[1].gsfTrack().pzMode(), tighteles[1].userFloat('trkModeParamCombinedEnergy_1') )
		    massz_trkmode_replacepmodewithbs = ( lep1 + lep2 ).M()
		    massZ_trkmode_replacepmodewithbsconstr_paramcombine.Fill( massz_trkmode_replacepmodewithbs )

		    corr1 = tighteles[0].userFloat('trkModeParamCombinedEnergy_2')/tighteles[0].gsfTrack().pMode()
		    corr2 = tighteles[1].userFloat('trkModeParamCombinedEnergy_2')/tighteles[1].gsfTrack().pMode()
		    lep1.SetPxPyPzE( corr1*tighteles[0].gsfTrack().pxMode(), corr1*tighteles[0].gsfTrack().pyMode(), corr1*tighteles[0].gsfTrack().pzMode(), tighteles[0].userFloat('trkModeParamCombinedEnergy_2') )
		    lep2.SetPxPyPzE( corr2*tighteles[1].gsfTrack().pxMode(), corr2*tighteles[1].gsfTrack().pyMode(), corr2*tighteles[1].gsfTrack().pzMode(), tighteles[1].userFloat('trkModeParamCombinedEnergy_2') )
		    massz_trkmode_replacepmodepca = ( lep1 + lep2 ).M()
		    massZ_trkmode_replacepmodewithpca_paramcombine.Fill( massz_trkmode_replacepmodepca )


c = ROOT.TCanvas('c','',1400,1000)
massZ_standard.Draw()
c.Print(path + '/' + massZ_standard.GetName() + '.png')
c.Clear()
massZ_trkmode.Draw()
c.Print(path + '/' + massZ_trkmode.GetName() + '.png')
c.Clear()
massZ_trkmode_paramcombine.Draw()
c.Print(path + '/' + massZ_trkmode_paramcombine.GetName() + '.png')
c.Clear()
massZ_trkmode_replacepmodewithpca_paramcombine.Draw()
c.Print(path + '/' + massZ_trkmode_replacepmodewithpca_paramcombine.GetName() + '.png' )
c.Clear()
massZ_trkmode_replacepmodewithbsconstr_paramcombine.Draw()
c.Print(path + '/' + massZ_trkmode_replacepmodewithbsconstr_paramcombine.GetName() + '.png' )
