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

userfloatnames = ['elClass', 
	'trkModeParamCombinedEnergy', 'trkModeParamCombinedEnergy', 'trkPMode', 'trkPModeErr',
	'momatPCAParamCombinedEnergy', 'momatPCAParamCombinedEnergyErr', 'momatPCAErr', 'momatPCA',
	'momatBSParamCombinedEnergy', 'momatBSParamCombinedEnergyErr', 'momatBSErr', 'momatBS']

list = []
names = []
# loop over events
for i,event in enumerate(events):
    #if i > options.maxEvents: break
    event.getByLabel (label, handle)
    eles = handle.product()

    for e in eles:
	userfloats = []
	for name in userfloatnames:
		if e.hasUserFloat (name ): 
			userfloats.append( '{:.3f}'.format(e.userFloat(name)) )
		else:
			userfloats.append( -999.9 )
	list.append([ 
		'{:.3f}'.format(e.p()), 
		'{:.3f}'.format(e.userFloat('ecalTrkEnergyPreCorr')), 
		'{:.3f}'.format(e.correctedEcalEnergy()),
		'{:.3f}'.format(e.trackMomentumAtVtx().R()), 
		'{:.3f}'.format(e.trackMomentumAtVtxWithConstraint().R()),
		])

	list[-1] = list[-1] + userfloats

names = [ 'e.p','ecalTrkEnergyPreCorr','correctedEcalEnergy','e.trackMomentumAtVtx().R()','e.trackMomentumAtVtxWithConstraint().R()' ] + userfloatnames

df = pandas.DataFrame( list, columns = names )
print df[['e.trackMomentumAtVtx().R()','momatPCA','e.trackMomentumAtVtxWithConstraint().R()','momatBS','elClass']]

myhtml = HTMLClass('Z to ee info for electron track study')
myhtml.section('')
string = myhtml.check()
string = string + '\n' + df[['e.trackMomentumAtVtx().R()','momatPCA','e.trackMomentumAtVtxWithConstraint().R()','momatBS','elClass']].to_html()
fp = open(path + '/test.html','w')
fp.write(string)
fp.close()
