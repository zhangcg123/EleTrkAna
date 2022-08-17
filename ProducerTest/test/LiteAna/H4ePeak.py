#! /usr/bin/env python3

import ROOT
import sys
from DataFormats.FWLite import Events, Handle
import pandas
sys.path.append('/afs/cern.ch/work/c/chenguan/private/pycommontool/')
from FileSystemClass import *
from HTMLClass import *

dirs = DirTree()
dirs.mkrootdir('test_higgspeak')
path = dirs.root

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.parseArguments()


########

def findhiggscandidate( electrons ):
	
	
	# 1st, make all sfos electron pairs
	n_zs = 0
	z_lepindex1 = []
	z_lepindex2 = []
	Zs = []
	Zmass = []
	Zeta = []
	Zpt = []
	Zphi = []
	for i,e1 in enumerate( electrons ):
		for j,e2 in enumerate( electrons ):
			if i == j:
				continue
			if e1.charge() + e2.charge() != 0:
				continue
			lepi = ROOT.TLorentzVector()
			lepj = ROOT.TLorentzVector()
			lepi.SetPtEtaPhiM( e1.pt(), e1.eta(), e1.phi(), e1.mass() )
			lepj.SetPtEtaPhiM( e2.pt(), e2.eta(), e2.phi(), e2.mass() )
			Z = lepi + lepj
			Zs.append( Z )
			Zpt.append( Z.Pt() )
			Zeta.append( Z.Eta() )
			Zphi.append( Z.Phi() )
			Zmass.append( Z.M() )
			z_lepindex1.append( i )
			z_lepindex2.append( j )
	# 2nd, make all possible z pairs
	for i,z1 in enumerate( Zs ):
		higgsleps = []
		for j,z2 in enumerate( Zs ):
			if i == j:
				continue
			i1 = z_lepindex1[i]
			i2 = z_lepindex2[i]
			j1 = z_lepindex1[j]
			j2 = z_lepindex2[j]
			
			if i1 == j1 or i1 == j2 or i2 == j1 or i2 == j2:
				continue

			lep_i1 = ROOT.TLorentzVector()
			lep_i2 = ROOT.TLorentzVector()
			lep_j1 = ROOT.TLorentzVector()
			lep_j2 = ROOT.TLorentzVector()

			lep_i1.SetPtEtaPhiM( electrons[i1].pt(), electrons[i1].eta(), electrons[i1].phi(), electrons[i1].mass() )
			lep_i2.SetPtEtaPhiM( electrons[i2].pt(), electrons[i2].eta(), electrons[i2].phi(), electrons[i2].mass() )
			lep_j1.SetPtEtaPhiM( electrons[j1].pt(), electrons[j1].eta(), electrons[j1].phi(), electrons[j1].mass() )
			lep_j2.SetPtEtaPhiM( electrons[j2].pt(), electrons[j2].eta(), electrons[j2].phi(), electrons[j2].mass() )

			zi = ROOT.TLorentzVector()
			zj = ROOT.TLorentzVector()
			zi.SetPtEtaPhiM( Zpt[i], Zeta[i], Zphi[i], Zmass[i] )
			zj.SetPtEtaPhiM( Zpt[j], Zeta[j], Zphi[j], Zmass[j] )

			Z1 = ROOT.TLorentzVector()
			Z2 = ROOT.TLorentzVector()
			Z1_index = 0
			Z2_index = 0
			Z1_lepindex = [0,0]
			Z2_lepindex = [0,0]

			# z1 and z2
			genzmass = 91.19
			if ( zi.M() - genzmass ) <= ( zj.M() - genzmass ):
				Z1_index = i
				Z2_index = j
				Z1 = zi
				Z2 = zj
			else:
				Z1_index = j
				Z2_index = i
				Z2 = zi
				Z1 = zj
			if Z1.M() < 40:
				continue
			# check pt
			leadingptcut = 20
			subleadingptcut = 10
			allpt = [ lep_i1.Pt(), lep_i2.Pt(), lep_j1.Pt(), lep_j2.Pt() ]
			allpt.sort()
			if allpt[-1] < leadingptcut or allpt[-2] < subleadingptcut:
				continue
			# all dr between leptons
			passedDeltaR = True
			for i,lep1 in enumerate([lep_i1,lep_i2,lep_j1,lep_j2]):
				for j,lep2 in enumerate([lep_i1,lep_i2,lep_j1,lep_j2] ):
					if i == j: continue
					if lep1.DeltaR( lep2 ) < 0.02:
						passedDeltaR = False
			if not passedDeltaR:
				continue
			# all os lepton pair mass
			passedAllOSPairMass = True
			lep1_tmp = ROOT.TLorentzVector()
			lep2_tmp = ROOT.TLorentzVector()
			for i,e1 in enumerate([electrons[i1],electrons[i2],electrons[j1],electrons[j2]]):
				for j,e2 in enumerate([electrons[i1],electrons[i2],electrons[j1],electrons[j2]], start = i + 1):
					if e1.charge() * e2.charge() < 0:
						lep1_tmp.SetPtEtaPhiM(e1.pt(), e1.eta(), e1.phi(), e1.mass() )
						lep2_tmp.SetPtEtaPhiM(e2.pt(), e2.eta(), e2.phi(), e2.mass() )
						m_tmp = ( lep1_tmp + lep2_tmp ).M()
						if m_tmp < 4.0:
							passedAllOSPairMass = False
			if not passedAllOSPairMass:
				continue
			higgsleps.append( electrons[i1] )
			higgsleps.append( electrons[i2] )
			higgsleps.append( electrons[j1] )
			higgsleps.append( electrons[j2] )
			
			break
		
		else:
			
			continue
		
		break

	return higgsleps








			# smart cut

events = Events (["ecalgsftrk6.root","ecalgsftrk5.root","ecalgsftrk4.root","ecalgsftrk3.root","ecalgsftrk1.root","ecalgsftrk0.root"])

handle_el  = Handle ("std::vector<pat::Electron>")
label_el = ("ProducerTest")

handle_bs = Handle ("std::vector<reco::Vertex>")
label_bs = ("offlinePrimaryVertices")

massZ_trkmode_paramcombine = ROOT.TH1F('massz_trkmode_paramcombine','',300,110,140)
massZ_momatpca_paramcombine = ROOT.TH1F('massZ_momatpca_paramcombine','',300,110,140)
massZ_momatbs_paramcombine = ROOT.TH1F('massZ_momatbs_paramcombine','',300,110,140)

massZ_momatpca = ROOT.TH1F('massZ_momatpca','',300,110,140)
massZ_momatbs = ROOT.TH1F('massZ_momatbs','',300,110,140)
massZ_ecal = ROOT.TH1F('massZ_ecal','',300,110,140)

# loop over events
for i,event in enumerate(events):
    
    if i%1000 == 0: print i
    #if i > options.maxEvents: break
    event.getByLabel (label_el, handle_el)
    eles = handle_el.product()
    
    event.getByLabel (label_bs, handle_bs)
    bs = handle_bs.product()

    goodpv = False
    for vtx in bs:
	    if vtx.isFake() or vtx.ndof() <= 4 or vtx.position().Rho() > 2.0 or abs(vtx.position().Z()) > 24.0:
		    continue
	    else:
		    goodpv = True
		    break
    if not goodpv: continue
	    
	    
    #build 4e events
    if len( eles ) < 4: continue
    
    higgsleps = findhiggscandidate( eles )	    
     
    if len( higgsleps ) != 4:
	    continue

    lep1 = ROOT.TLorentzVector()
    lep2 = ROOT.TLorentzVector()
    lep3 = ROOT.TLorentzVector()
    lep4 = ROOT.TLorentzVector()
    
    # trkmode combined	--- as the base	    
    corr1 = higgsleps[0].userFloat('trkModeParamCombinedEnergy')/higgsleps[0].userFloat('trkPMode')
    corr2 = higgsleps[1].userFloat('trkModeParamCombinedEnergy')/higgsleps[1].userFloat('trkPMode')
    corr3 = higgsleps[2].userFloat('trkModeParamCombinedEnergy')/higgsleps[2].userFloat('trkPMode')
    corr4 = higgsleps[3].userFloat('trkModeParamCombinedEnergy')/higgsleps[3].userFloat('trkPMode')
    lep1.SetPxPyPzE( higgsleps[0].userFloat('trkPxMode')*corr1, higgsleps[0].userFloat('trkPyMode')*corr1, corr1*higgsleps[0].userFloat('trkPzMode'), higgsleps[0].userFloat('trkModeParamCombinedEnergy') )
    lep2.SetPxPyPzE( higgsleps[1].userFloat('trkPxMode')*corr2, higgsleps[1].userFloat('trkPyMode')*corr2, corr2*higgsleps[1].userFloat('trkPzMode'), higgsleps[1].userFloat('trkModeParamCombinedEnergy') ) 
    lep3.SetPxPyPzE( higgsleps[2].userFloat('trkPxMode')*corr3, higgsleps[2].userFloat('trkPyMode')*corr3, corr3*higgsleps[2].userFloat('trkPzMode'), higgsleps[2].userFloat('trkModeParamCombinedEnergy') )
    lep4.SetPxPyPzE( higgsleps[3].userFloat('trkPxMode')*corr4, higgsleps[3].userFloat('trkPyMode')*corr4, corr4*higgsleps[3].userFloat('trkPzMode'), higgsleps[3].userFloat('trkModeParamCombinedEnergy') ) 
    massz_trkmode_paramcombine = ( lep1 + lep2 + lep3 + lep4 ).M()
    massZ_trkmode_paramcombine.Fill( massz_trkmode_paramcombine )
    
    # momatpca combined
    if (higgsleps[0].userFloat('momatPCAParamCombinedEnergyErr') < higgsleps[0].userFloat('trkModeParamCombinedEnergyErr')):
	corr1 = higgsleps[0].userFloat('momatPCAParamCombinedEnergy')/higgsleps[0].userFloat('momatPCA')
	lep1.SetPxPyPzE(corr1*higgsleps[0].userFloat('momxatPCA'), corr1*higgsleps[0].userFloat('momyatPCA'), corr1*higgsleps[0].userFloat('momzatPCA'), higgsleps[0].userFloat('momatPCAParamCombinedEnergy') )
    if (higgsleps[1].userFloat('momatPCAParamCombinedEnergyErr') < higgsleps[1].userFloat('trkModeParamCombinedEnergyErr')):
	corr2 = higgsleps[1].userFloat('momatPCAParamCombinedEnergy')/higgsleps[1].userFloat('momatPCA')
	lep2.SetPxPyPzE(corr2*higgsleps[1].userFloat('momxatPCA'), corr2*higgsleps[1].userFloat('momyatPCA'), corr2*higgsleps[1].userFloat('momzatPCA'), higgsleps[1].userFloat('momatPCAParamCombinedEnergy') ) 
    if (higgsleps[2].userFloat('momatPCAParamCombinedEnergyErr') < higgsleps[2].userFloat('trkModeParamCombinedEnergyErr')):
	corr3 = higgsleps[2].userFloat('momatPCAParamCombinedEnergy')/higgsleps[2].userFloat('momatPCA')
	lep3.SetPxPyPzE(corr3*higgsleps[2].userFloat('momxatPCA'), corr3*higgsleps[2].userFloat('momyatPCA'), corr3*higgsleps[2].userFloat('momzatPCA'), higgsleps[2].userFloat('momatPCAParamCombinedEnergy') )
    if (higgsleps[3].userFloat('momatPCAParamCombinedEnergyErr') < higgsleps[3].userFloat('trkModeParamCombinedEnergyErr')):
	corr4 = higgsleps[3].userFloat('momatPCAParamCombinedEnergy')/higgsleps[3].userFloat('momatPCA')
	lep4.SetPxPyPzE( corr4*higgsleps[3].userFloat('momxatPCA'), corr4*higgsleps[3].userFloat('momyatPCA'), corr4*higgsleps[3].userFloat('momzatPCA'), higgsleps[3].userFloat('momatPCAParamCombinedEnergy') )
    massz_momatpca_paramcombine = ( lep1 + lep2 + lep3 + lep4 ).M()
    massZ_momatpca_paramcombine.Fill( massz_momatpca_paramcombine )

    # bs combined
    if ( higgsleps[0].userFloat('momatBSParamCombinedEnergyErr') < higgsleps[0].userFloat('trkModeParamCombinedEnergyErr') ):
	corr1 = higgsleps[0].userFloat('momatBSParamCombinedEnergy')/higgsleps[0].userFloat('momatBS')
	lep1.SetPxPyPzE(corr1*higgsleps[0].userFloat('momxatBS'), corr1*higgsleps[0].userFloat('momyatBS'), corr1*higgsleps[0].userFloat('momzatBS'), higgsleps[0].userFloat('momatBSParamCombinedEnergy') )
    if ( higgsleps[1].userFloat('momatBSParamCombinedEnergyErr') < higgsleps[1].userFloat('trkModeParamCombinedEnergyErr') ):
	corr2 = higgsleps[1].userFloat('momatBSParamCombinedEnergy')/higgsleps[1].userFloat('momatBS')
	lep2.SetPxPyPzE(corr2*higgsleps[1].userFloat('momxatBS'), corr2*higgsleps[1].userFloat('momyatBS'), corr2*higgsleps[1].userFloat('momzatBS'), higgsleps[1].userFloat('momatBSParamCombinedEnergy') )
    if ( higgsleps[2].userFloat('momatBSParamCombinedEnergyErr') < higgsleps[2].userFloat('trkModeParamCombinedEnergyErr') ):
	corr3 = higgsleps[2].userFloat('momatBSParamCombinedEnergy')/higgsleps[2].userFloat('momatBS')
	lep3.SetPxPyPzE(corr3*higgsleps[2].userFloat('momxatBS'), corr3*higgsleps[2].userFloat('momyatBS'), corr3*higgsleps[2].userFloat('momzatBS'), higgsleps[2].userFloat('momatBSParamCombinedEnergy') )
    if ( higgsleps[3].userFloat('momatBSParamCombinedEnergyErr') < higgsleps[3].userFloat('trkModeParamCombinedEnergyErr') ):
	corr4 = higgsleps[3].userFloat('momatBSParamCombinedEnergy')/higgsleps[3].userFloat('momatBS')
	lep4.SetPxPyPzE(corr4*higgsleps[3].userFloat('momxatBS'), corr4*higgsleps[3].userFloat('momyatBS'), corr4*higgsleps[3].userFloat('momzatBS'), higgsleps[3].userFloat('momatBSParamCombinedEnergy') )
    massz_momatbs_paramcombine = ( lep1 + lep2 + lep3 + lep4 ).M()
    massZ_momatbs_paramcombine.Fill( massz_momatbs_paramcombine )
    
    # ecal only 
    corr1 = higgsleps[0].correctedEcalEnergy()/higgsleps[0].userFloat('trkPMode')
    lep1.SetPxPyPzE(corr1*higgsleps[0].userFloat('trkPxMode'),corr1*higgsleps[0].userFloat('trkPyMode'),corr1*higgsleps[0].userFloat('trkPzMode'),higgsleps[0].correctedEcalEnergy() )
    corr2 = higgsleps[1].correctedEcalEnergy()/higgsleps[1].userFloat('trkPMode')
    lep2.SetPxPyPzE(corr2*higgsleps[1].userFloat('trkPxMode'),corr2*higgsleps[1].userFloat('trkPyMode'),corr2*higgsleps[1].userFloat('trkPzMode'),higgsleps[1].correctedEcalEnergy() )
    corr3 = higgsleps[2].correctedEcalEnergy()/higgsleps[2].userFloat('trkPMode')
    lep3.SetPxPyPzE(corr3*higgsleps[2].userFloat('trkPxMode'),corr3*higgsleps[2].userFloat('trkPyMode'),corr3*higgsleps[2].userFloat('trkPzMode'),higgsleps[2].correctedEcalEnergy() )
    corr4 = higgsleps[3].correctedEcalEnergy()/higgsleps[3].userFloat('trkPMode')
    lep4.SetPxPyPzE(corr4*higgsleps[3].userFloat('trkPxMode'),corr4*higgsleps[3].userFloat('trkPyMode'),corr4*higgsleps[3].userFloat('trkPzMode'),higgsleps[3].correctedEcalEnergy() )
    massz_ecal = (lep1 + lep2 + lep3 + lep4 ).M()
    massZ_ecal.Fill(massz_ecal)

    '''
    lep1.SetPxPyPzE( higgsleps[0].userFloat('trkPxMode'), higgsleps[0].userFloat('trkPyMode'), higgsleps[0].userFloat('trkPzMode'),higgsleps[0].userFloat('trkPMode') )
    lep2.SetPxPyPzE( higgsleps[1].userFloat('trkPxMode'), higgsleps[1].userFloat('trkPyMode'), higgsleps[1].userFloat('trkPzMode'),higgsleps[1].userFloat('trkPMode') )
    # pac only
    if higgsleps[0].userFloat('momatPCAErr') < higgsleps[0].userFloat('trkPModeErr'):
	    lep1.SetPxPyPzE(higgsleps[0].userFloat('momxatPCA'),higgsleps[0].userFloat('momyatPCA'),higgsleps[0].userFloat('momzatPCA'),higgsleps[0].userFloat('momatPCA') )
    if higgsleps[1].userFloat('momatPCAErr') < higgsleps[1].userFloat('trkPModeErr'):
	    lep2.SetPxPyPzE(higgsleps[1].userFloat('momxatPCA'),higgsleps[1].userFloat('momyatPCA'),higgsleps[1].userFloat('momzatPCA'),higgsleps[1].userFloat('momatPCA') )
    massz_momatpca = (lep1+lep2).M()
    massZ_momatpca.Fill(massz_momatpca)
    # bs only
    if higgsleps[0].userFloat('momatBSErr') < higgsleps[0].userFloat('trkPModeErr'):
	    lep1.SetPxPyPzE(higgsleps[0].userFloat('momxatBS'),higgsleps[0].userFloat('momyatBS'),higgsleps[0].userFloat('momzatBS'),higgsleps[0].userFloat('momatBS') )
    if higgsleps[1].userFloat('momatBSErr') < higgsleps[1].userFloat('trkPModeErr'):
	    lep2.SetPxPyPzE(higgsleps[1].userFloat('momxatBS'),higgsleps[1].userFloat('momyatBS'),higgsleps[1].userFloat('momzatBS'),higgsleps[1].userFloat('momatBS') )
    massz_momatbs = (lep1+lep2).M()
    massZ_momatbs.Fill(massz_momatbs)
    '''

x = ROOT.RooRealVar('x','',110,140)
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
	#x.setBins(1000)
	#zm = ROOT.RooRealVar('zm','',91.153)
	#zw = ROOT.RooRealVar('zw','',2.493)
	#zshape = ROOT.RooBreitWigner('zshape','',x,zm,zw)

	m = ROOT.RooRealVar('m','',124,122,126)
	s = ROOT.RooRealVar('s','',2,1,3)
	n1 = ROOT.RooRealVar('n1','',10,0,70)
	a1 = ROOT.RooRealVar('a1','',1,0,10)
	n2 = ROOT.RooRealVar('n2','',10,0,70) 
	a2 = ROOT.RooRealVar('a2','',1,0,70)
	dcb = ROOT.RooDoubleCB('dcb','',x,m,s,a1,n1,a2,n2)
	#fft = ROOT.RooFFTConvPdf('fft','',x,zshape,dcb)
	#fft.setBufferFraction(0.2)

	dcb.fitTo(d)
	chi2_var = ROOT.RooChi2Var('chi2', '', dcb, d, ROOT.RooFit.DataError(ROOT.RooAbsData.Expected))
	frame = x.frame(ROOT.RooFit.Bins(30))
	frame.SetTitle('')
	d_to_plot[i].plotOn(frame)
	dcb.plotOn(frame)
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
	lx.DrawLatex(0.7,0.8,"mean="+str(m.getVal())[0:6]+"+/-"+str(m.getError())[0:4])
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
'''
frametrk = x.frame()
frametrk.SetTitle('')
d_trkbs_plot.plotOn(frametrk,ROOT.RooFit.MarkerColor(2))
d_trkpca_plot.plotOn(frametrk,ROOT.RooFit.MarkerColor(1))
ctot.Clear()
frametrk.Draw()
ctot.Print(path + '/trkpcabs.png' )
'''
