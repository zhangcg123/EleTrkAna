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

handle_el  = Handle ("std::vector<pat::Electron>")
label_el = ("ProducerTest")

handle_bs = Handle ("reco::BeamSpot")
label_bs = ("offlineBeamSpot")

massZ_trkmode_paramcombine = ROOT.TH1F('massz_trkmode_paramcombine','',100,70,110)
massZ_momatpca_paramcombine = ROOT.TH1F('massZ_momatpca_paramcombine','',100,70,110)
massZ_momatbs_paramcombine = ROOT.TH1F('massZ_momatbs_paramcombine','',100,70,110)

# loop over events
for i,event in enumerate(events):
    #if i > options.maxEvents: break
    event.getByLabel (label_el, handle_el)
    eles = handle_el.product()
    
    event.getByLabel (label_bs, handel_bs)
    bs = handel_el.product()

    lep1 = ROOT.TLorentzVector()
    lep2 = ROOT.TLorentzVector()
    for e in eles:
	    
	    if len( eles ) == 2 and eles[0].charge() * eles[1].charge() < 0:


		    lep1.SetPtEtaPhiM( eles[0].pt(), eles[0].eta(), eles[0].phi(), 0.0005 )
		    lep2.SetPtEtaPhiM( eles[1].pt(), eles[1].eta(), eles[1].phi(), 0.0005 ) 
		    massz_trkmode_paramcombine = ( lep1 + lep2 ).M()
		    massZ_trkmode_paramcombine.Fill( massz_trkmode_paramcombine )

		    if eles[0].hasUserFloat('momatPCAParamCombinedEnergy'):
		    	corr1 = eles[0].userFloat('momatPCAParamCombinedEnergy')/eles[0].userFloat('momatPCA')
			lep1.SetPxPyPzE(corr1*eles[0].trackMomentumAtVtx().X(), corr1*eles[0].trackMomentumAtVtx().Y(), corr1*eles[0].trackMomentumAtVtx().Z(), eles[0].userFloat('momatPCAParamCombinedEnergy') )
		    if eles[1].hasUserFloat('momatPCAParamCombinedEnergy'):
			corr2 = eles[1].userFloat('momatPCAParamCombinedEnergy')/eles[1].userFloat('momatPCA')
			lep2.SetPxPyPzE(corr2*eles[1].trackMomentumAtVtx().X(), corr2*eles[1].trackMomentumAtVtx().Y(), corr2*eles[1].trackMomentumAtVtx().Z(), eles[1].userFloat('momatPCAParamCombinedEnergy') )
		    massz_momatpca_paramcombine = ( lep1 + lep2 ).M()
		    massZ_momatpca_paramcombine.Fill( massz_momatpca_paramcombine )

		    if eles[0].hasUserFloat('momatBSParamCombinedEnergy'):
		    	corr1 = eles[0].userFloat('momatBSParamCombinedEnergy')/eles[0].userFloat('momatBS')
			lep1.SetPxPyPzE(corr1*eles[0].trackMomentumAtVtxWithConstraint().X(), corr1*eles[0].trackMomentumAtVtxWithConstraint().Y(), corr1*eles[0].trackMomentumAtVtxWithConstraint().Z(), eles[0].userFloat('momatBSParamCombinedEnergy') )
		    if eles[1].hasUserFloat('momatBSParamCombinedEnergy'):
			corr2 = eles[1].userFloat('momatBSParamCombinedEnergy')/eles[1].userFloat('momatBS')
			lep2.SetPxPyPzE(corr2*eles[1].trackMomentumAtVtxWithConstraint().X(), corr2*eles[1].trackMomentumAtVtxWithConstraint().Y(), corr2*eles[1].trackMomentumAtVtxWithConstraint().Z(), eles[1].userFloat('momatBSParamCombinedEnergy') )
		    massz_momatbs_paramcombine = ( lep1 + lep2 ).M()
		    massZ_momatbs_paramcombine.Fill( massz_momatbs_paramcombine )


x = ROOT.RooRealVar('x','',70,110)
d_norm = ROOT.RooDataHist('d_norm','',ROOT.RooArgList(x),ROOT.RooFit.Import(massZ_trkmode_paramcombine))
d_pca = ROOT.RooDataHist('d_pca','',ROOT.RooArgList(x),ROOT.RooFit.Import(massZ_momatpca_paramcombine))
d_bs = ROOT.RooDataHist('d_bs','',ROOT.RooArgList(x),ROOT.RooFit.Import(massZ_momatbs_paramcombine))

d_norm_plot = ROOT.RooDataHist('d_norm_plot','',ROOT.RooArgList(x),ROOT.RooFit.Import(massZ_trkmode_paramcombine.Rebin(2)))
d_pca_plot = ROOT.RooDataHist('d_pca_plot','',ROOT.RooArgList(x),ROOT.RooFit.Import(massZ_momatpca_paramcombine.Rebin(2)))  
d_bs_plot = ROOT.RooDataHist('d_bs_plot','',ROOT.RooArgList(x),ROOT.RooFit.Import(massZ_momatbs_paramcombine.Rebin(2)))   

d_to_fit = [d_norm, d_pca, d_bs]
d_to_plot = [d_norm_plot, d_pca_plot, d_bs_plot]

for i,d in enumerate(d_to_fit):
	if i == 0: continue

	m = ROOT.RooRealVar('m','',90,85,92)
	s = ROOT.RooRealVar('s','',1,0,10)
	n1 = ROOT.RooRealVar('n1','',10,0,70)
	a1 = ROOT.RooRealVar('a1','',1,0,10)
	n2 = ROOT.RooRealVar('n2','',10,0,70) 
	a2 = ROOT.RooRealVar('a2','',1,0,70)
	dcb = ROOT.RooDoubleCB('dcb','',x,m,s,a1,n1,a2,n2)

	dcb.fitTo(d)
	chi2_var = ROOT.RooChi2Var('chi2', '', dcb, d, ROOT.RooFit.DataError(ROOT.RooAbsData.Expected))
	frame = x.frame(ROOT.RooFit.Bins(30))
	frame.SetTitle('')
	d_to_plot[i].plotOn(frame)
	dcb.plotOn(frame)
	c = ROOT.TCanvas('c','',1400,1000)
	c.cd()
	frame.Draw()
	lx = ROOT.TLatex()
	lx.SetNDC()
	lx.SetTextSize(0.05)
	lx.SetTextFont(42)
	lx.SetTextAlign(23)
	lx.DrawLatex(0.7,0.9,"#chi^{2}/dof="+str(chi2_var.getVal()/float(d.numEntries()))[0:4])
	lx.DrawLatex(0.7,0.8,"mean="+str(m.getVal())[0:4]+"+/-"+str(m.getError())[0:4])
	lx.DrawLatex(0.7,0.7,"sigma="+str(s.getVal())[0:4]+"+/-"+str(s.getError())[0:4])
	c.Print(path + '/' + d.GetName() + '.png')
