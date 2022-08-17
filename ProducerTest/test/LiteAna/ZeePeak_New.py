#! /usr/bin/env python3

import ROOT
import sys
from DataFormats.FWLite import Events, Handle
import pandas
sys.path.append('/afs/cern.ch/work/c/chenguan/private/pycommontool/')
from FileSystemClass import *
from HTMLClass import *

dirs = DirTree()
dirs.mkrootdir('test_removedoublecout')
path = dirs.root

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.parseArguments()

rootlist = ["ecalgsftrk1.root","ecalgsftrk2.root","ecalgsftrk3.root","ecalgsftrk4.root","ecalgsftrk5.root","ecalgsftrk6.root"]
for idx in range(len(rootlist)):
	rootlist[idx] = '/eos/user/c/chenguan/Tuples/EleTrkStudies/PrivateFullSimDYJetsToLLM50/' + rootlist[idx]

events = Events (rootlist)

handle_el  = Handle ("std::vector<pat::Electron>")
label_el = ("ProducerTest")

handle_bs = Handle ("std::vector<reco::Vertex>")
label_bs = ("offlinePrimaryVertices")

massZ_trkmode_paramcombine = ROOT.TH1F('massz_trkmode_paramcombine','',800,70,110)
massZ_momatpca_paramcombine = ROOT.TH1F('massZ_momatpca_paramcombine','',800,70,110)
massZ_momatbs_paramcombine = ROOT.TH1F('massZ_momatbs_paramcombine','',800,70,110)

massZ_momatpca = ROOT.TH1F('massZ_momatpca','',800,70,110)
massZ_momatbs = ROOT.TH1F('massZ_momatbs','',800,70,110)
massZ_ecal = ROOT.TH1F('massZ_ecal','',800,70,110)

# loop over events
for i,event in enumerate(events):
    
	if i%1000 == 0: print i
	#if i > options.maxEvents: break

	event.getByLabel (label_el, handle_el)
	eles = handle_el.product()

	event.getByLabel (label_bs, handle_bs)
	bs = handle_bs.product()

	if len( eles ) != 2 or eles[0].charge() * eles[1].charge() > 0:
	    continue

	goodpv = False
	for vtx in bs:
	    if vtx.isFake() or vtx.ndof() <= 4 or vtx.position().Rho() > 2.0 or abs(vtx.position().Z()) > 24.0:
		    continue
	    else:
		    goodpv = True
		    break
	if not goodpv: continue
	    
	lep1 = ROOT.TLorentzVector()
	lep2 = ROOT.TLorentzVector()
		    
	# trkmode combined	--- as the base	    
	corr1 = eles[0].userFloat('trkModeParamCombinedEnergy')/eles[0].userFloat('trkPMode')
	corr2 = eles[1].userFloat('trkModeParamCombinedEnergy')/eles[1].userFloat('trkPMode')
	lep1.SetPxPyPzE( eles[0].userFloat('trkPxMode')*corr1, eles[0].userFloat('trkPyMode')*corr1, corr1*eles[0].userFloat('trkPzMode'), eles[0].userFloat('trkModeParamCombinedEnergy') )
	lep2.SetPxPyPzE( eles[1].userFloat('trkPxMode')*corr2, eles[1].userFloat('trkPyMode')*corr2, corr2*eles[1].userFloat('trkPzMode'), eles[1].userFloat('trkModeParamCombinedEnergy') ) 
	massz_trkmode_paramcombine = ( lep1 + lep2 ).M()
	massZ_trkmode_paramcombine.Fill( massz_trkmode_paramcombine )

	# momatpca combined
	if (eles[0].userFloat('momatPCAParamCombinedEnergyErr') < eles[0].userFloat('trkModeParamCombinedEnergyErr')):
		corr1 = eles[0].userFloat('momatPCAParamCombinedEnergy')/eles[0].userFloat('momatPCA')
		lep1.SetPxPyPzE(corr1*eles[0].userFloat('momxatPCA'), corr1*eles[0].userFloat('momyatPCA'), corr1*eles[0].userFloat('momzatPCA'), eles[0].userFloat('momatPCAParamCombinedEnergy') )
	if (eles[1].userFloat('momatPCAParamCombinedEnergyErr') < eles[1].userFloat('trkModeParamCombinedEnergyErr')):
		corr2 = eles[1].userFloat('momatPCAParamCombinedEnergy')/eles[1].userFloat('momatPCA')
		lep2.SetPxPyPzE(corr2*eles[1].userFloat('momxatPCA'), corr2*eles[1].userFloat('momyatPCA'), corr2*eles[1].userFloat('momzatPCA'), eles[1].userFloat('momatPCAParamCombinedEnergy') )
	massz_momatpca_paramcombine = ( lep1 + lep2 ).M()
	massZ_momatpca_paramcombine.Fill( massz_momatpca_paramcombine )

	# bs combined
	if ( eles[0].userFloat('momatBSParamCombinedEnergyErr') < eles[0].userFloat('trkModeParamCombinedEnergyErr') ):
		corr1 = eles[0].userFloat('momatBSParamCombinedEnergy')/eles[0].userFloat('momatBS')
		lep1.SetPxPyPzE(corr1*eles[0].userFloat('momxatBS'), corr1*eles[0].userFloat('momyatBS'), corr1*eles[0].userFloat('momzatBS'), eles[0].userFloat('momatBSParamCombinedEnergy') )
	if ( eles[1].userFloat('momatBSParamCombinedEnergyErr') < eles[1].userFloat('trkModeParamCombinedEnergyErr') ):
		corr2 = eles[1].userFloat('momatBSParamCombinedEnergy')/eles[1].userFloat('momatBS')
		lep2.SetPxPyPzE(corr2*eles[1].userFloat('momxatBS'), corr2*eles[1].userFloat('momyatBS'), corr2*eles[1].userFloat('momzatBS'), eles[1].userFloat('momatBSParamCombinedEnergy') )
	massz_momatbs_paramcombine = ( lep1 + lep2 ).M()
	massZ_momatbs_paramcombine.Fill( massz_momatbs_paramcombine )

	# ecal only 
	corr1 = eles[0].correctedEcalEnergy()/eles[0].userFloat('trkPMode')
	lep1.SetPxPyPzE(corr1*eles[0].userFloat('trkPxMode'),corr1*eles[0].userFloat('trkPyMode'),corr1*eles[0].userFloat('trkPzMode'),eles[0].correctedEcalEnergy() )
	corr2 = eles[1].correctedEcalEnergy()/eles[1].userFloat('trkPMode')
	lep2.SetPxPyPzE(corr2*eles[1].userFloat('trkPxMode'),corr2*eles[1].userFloat('trkPyMode'),corr2*eles[1].userFloat('trkPzMode'),eles[1].correctedEcalEnergy() )
	massz_ecal = (lep1 + lep2).M()
	massZ_ecal.Fill(massz_ecal)

	lep1.SetPxPyPzE( eles[0].userFloat('trkPxMode'), eles[0].userFloat('trkPyMode'), eles[0].userFloat('trkPzMode'),eles[0].userFloat('trkPMode') )
	lep2.SetPxPyPzE( eles[1].userFloat('trkPxMode'), eles[1].userFloat('trkPyMode'), eles[1].userFloat('trkPzMode'),eles[1].userFloat('trkPMode') )
	# pac only
	if eles[0].userFloat('momatPCAErr') < eles[0].userFloat('trkPModeErr'):
		lep1.SetPxPyPzE(eles[0].userFloat('momxatPCA'),eles[0].userFloat('momyatPCA'),eles[0].userFloat('momzatPCA'),eles[0].userFloat('momatPCA') )
	if eles[1].userFloat('momatPCAErr') < eles[1].userFloat('trkPModeErr'):
		lep2.SetPxPyPzE(eles[1].userFloat('momxatPCA'),eles[1].userFloat('momyatPCA'),eles[1].userFloat('momzatPCA'),eles[1].userFloat('momatPCA') )
	massz_momatpca = (lep1+lep2).M()
	massZ_momatpca.Fill(massz_momatpca)
	# bs only
	if eles[0].userFloat('momatBSErr') < eles[0].userFloat('trkPModeErr'):
		lep1.SetPxPyPzE(eles[0].userFloat('momxatBS'),eles[0].userFloat('momyatBS'),eles[0].userFloat('momzatBS'),eles[0].userFloat('momatBS') )
	if eles[1].userFloat('momatBSErr') < eles[1].userFloat('trkPModeErr'):
		lep2.SetPxPyPzE(eles[1].userFloat('momxatBS'),eles[1].userFloat('momyatBS'),eles[1].userFloat('momzatBS'),eles[1].userFloat('momatBS') )
	massz_momatbs = (lep1+lep2).M()
	massZ_momatbs.Fill(massz_momatbs)


x = ROOT.RooRealVar('x','',70,110)
d_norm = ROOT.RooDataHist('d_norm','',ROOT.RooArgList(x),ROOT.RooFit.Import(massZ_trkmode_paramcombine))
d_pca = ROOT.RooDataHist('d_pca','',ROOT.RooArgList(x),ROOT.RooFit.Import(massZ_momatpca_paramcombine))
d_bs = ROOT.RooDataHist('d_bs','',ROOT.RooArgList(x),ROOT.RooFit.Import(massZ_momatbs_paramcombine))
d_ecal = ROOT.RooDataHist('d_ecal','',ROOT.RooArgList(x),ROOT.RooFit.Import(massZ_ecal))

d_norm_plot = ROOT.RooDataHist('d_norm_plot','',ROOT.RooArgList(x),ROOT.RooFit.Import(massZ_trkmode_paramcombine.Rebin(8)))
d_pca_plot = ROOT.RooDataHist('d_pca_plot','',ROOT.RooArgList(x),ROOT.RooFit.Import(massZ_momatpca_paramcombine.Rebin(8)))  
d_bs_plot = ROOT.RooDataHist('d_bs_plot','',ROOT.RooArgList(x),ROOT.RooFit.Import(massZ_momatbs_paramcombine.Rebin(8)))   
d_ecal_plot = ROOT.RooDataHist('d_ecal_plot','',ROOT.RooArgList(x),ROOT.RooFit.Import(massZ_ecal.Rebin(8)))
d_trkpca_plot = ROOT.RooDataHist('d_trkpca_plot','',ROOT.RooArgList(x),ROOT.RooFit.Import(massZ_momatpca.Rebin(16)))
d_trkbs_plot = ROOT.RooDataHist('d_trkbs_plot','',ROOT.RooArgList(x),ROOT.RooFit.Import(massZ_momatbs.Rebin(16)))

d_to_fit = [d_norm, d_pca, d_bs, d_ecal]
d_to_plot = [d_norm_plot, d_pca_plot, d_bs_plot, d_ecal_plot]

dcbs = []

for i,d in enumerate(d_to_fit):
	#if i == 0: continue
	x.setBins(1000)
	zm = ROOT.RooRealVar('zm','',91.153)
	zw = ROOT.RooRealVar('zw','',2.493)
	zshape = ROOT.RooBreitWigner('zshape','',x,zm,zw)

	m = ROOT.RooRealVar('m','',-0.3,-3,3)
	s = ROOT.RooRealVar('s','',1,0,10)
	n1 = ROOT.RooRealVar('n1','',10,0,70)
	a1 = ROOT.RooRealVar('a1','',1,0,10)
	n2 = ROOT.RooRealVar('n2','',10,0,70) 
	a2 = ROOT.RooRealVar('a2','',1,0,70)
	dcb = ROOT.RooDoubleCB('dcb','',x,m,s,a1,n1,a2,n2)
	fft = ROOT.RooFFTConvPdf('fft','',x,zshape,dcb)
	fft.setBufferFraction(0.2)

	fft.fitTo(d)
	chi2_var = ROOT.RooChi2Var('chi2', '', fft, d, ROOT.RooFit.DataError(ROOT.RooAbsData.Expected))
	frame = x.frame(ROOT.RooFit.Bins(30))
	frame.SetTitle('')
	d_to_plot[i].plotOn(frame)
	fft.plotOn(frame)
	c = ROOT.TCanvas('c','',1400,1000)
	c.SetTicks(1,1)
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

ctot = ROOT.TCanvas('ctot','',1400,1000)
ctot.cd()
ctot.SetTicks(1,1)
frametot = x.frame()
frametot.SetTitle('')
d_bs_plot.plotOn(frametot,ROOT.RooFit.MarkerColor(2))
d_pca_plot.plotOn(frametot,ROOT.RooFit.MarkerColor(1))
frametot.Draw()
ctot.Print(path + '/pcainblack_bsinread.png' )

frametrk = x.frame()
frametrk.SetTitle('')
d_trkbs_plot.plotOn(frametrk,ROOT.RooFit.MarkerColor(2))
d_trkpca_plot.plotOn(frametrk,ROOT.RooFit.MarkerColor(1))
ctot.Clear()
frametrk.Draw()
ctot.Print(path + '/trkpcabs.png' )
