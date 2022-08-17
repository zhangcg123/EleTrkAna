#! /usr/bin/env python3

import ROOT
import sys
from DataFormats.FWLite import Events, Handle
import pandas
sys.path.append('/afs/cern.ch/work/c/chenguan/private/pycommontool/')
from FileSystemClass import *
from HTMLClass import *

dirs = DirTree()
dirs.mkrootdir('test_electron_perr_vs_p_trkpmodeatpca_endcap')
path = dirs.root

def histo_to_graph( h ):
	g = ROOT.TGraph()
	for i in range(h.GetNbinsX()):
		g.SetPoint( i, h.GetBinCenter(i+1), h.GetBinContent(i+1) )
	return g

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.parseArguments()

events = Events (["ecalgsftrk1.root","ecalgsftrk2.root","ecalgsftrk3.root","ecalgsftrk4.root","ecalgsftrk5.root","ecalgsftrk6.root"])

handle_el  = Handle ("std::vector<pat::Electron>")
label_el = ("ProducerTest")

handle_bs = Handle ("std::vector<reco::Vertex>")
label_bs = ("offlinePrimaryVertices")

ecal_only = ROOT.TH2F('ecal_only','ecal',200,0,200,100,0,0.2)

trkmode = ROOT.TH2F('trkmode','',200,0,200,100,0,0.2)
trkmode_paramcombine = ROOT.TH2F('trkmode_paramcombine','',200,0,200,100,0,0.2)

momatpca = ROOT.TH2F('momatpca','trk_mom_at_pca',200,0,200,100,0,0.2)
momatpca_paramcombine = ROOT.TH2F('momatpca_paramcombine','comb_mom_at_pca',200,0,200,100,0,0.2)

momatbs = ROOT.TH2F('momatbs','trk_mom_at_pca',200,0,200,100,0,0.2)
momatbs_paramcombine = ROOT.TH2F('momatbs_paramcombine','comb_mom_at_pca',200,0,200,100,0,0.2)

# loop over events
for i,event in enumerate(events):
    if i%1000==0: print i
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
	    	     
    for e in eles:
	    if abs(e.eta()) < 1.4 : continue
	    trkmode.Fill( e.userFloat('trkPMode'), e.userFloat('trkPModeErr')/e.userFloat('trkPMode') )
	    trkmode_paramcombine.Fill( e.userFloat('trkModeParamCombinedEnergy'), e.userFloat('trkModeParamCombinedEnergyErr')/e.userFloat('trkModeParamCombinedEnergy') )
	    ecal_only.Fill( e.correctedEcalEnergy(), e.correctedEcalEnergyError()/e.correctedEcalEnergy() )
	    momatpca.Fill( e.userFloat('momatPCA'), e.userFloat('momatPCAErr')/e.userFloat('momatPCA') )
	    momatbs.Fill( e.userFloat('momatBS'), e.userFloat('momatBSErr')/e.userFloat('momatBS') )
	    momatpca_paramcombine.Fill( e.userFloat('momatPCAParamCombinedEnergy'), e.userFloat('momatPCAParamCombinedEnergyErr')/e.userFloat('momatPCAParamCombinedEnergy') )
	    momatbs_paramcombine.Fill( e.userFloat('momatBSParamCombinedEnergy'), e.userFloat('momatBSParamCombinedEnergyErr')/e.userFloat('momatBSParamCombinedEnergy') )

profile_ecal = ecal_only.ProfileX()
project_ecal = ecal_only.ProjectionX()

profile_trkmode = trkmode.ProfileX()
profile_trkmode_combine = trkmode_paramcombine.ProfileX()

project_trkmode = trkmode.ProjectionX()
project_trkmode_combine = trkmode_paramcombine.ProjectionX()

profile_momatpca = momatpca.ProfileX()
project_momatpca = momatpca.ProjectionX()
profile_momatpca_combine = momatpca_paramcombine.ProfileX()
project_momatpca_combine = momatpca_paramcombine.ProjectionX()

profile_momatbs = momatbs.ProfileX()
project_momatbs = momatbs.ProjectionX()
profile_momatbs_combine = momatbs_paramcombine.ProfileX()
project_momatbs_combine = momatbs_paramcombine.ProjectionX()

gs = []

for h in [profile_ecal, profile_momatpca, profile_momatbs, profile_momatbs_combine, profile_momatpca_combine, profile_trkmode_combine, profile_trkmode]:
	h.Rebin(4)	
	g = histo_to_graph( h )
	g.SetMarkerStyle(20)
	g.SetMarkerSize(1)
	gs.append( g )

	c = ROOT.TCanvas('c','',1400,1000)
	g.Draw('ap')
	c.Print(path + '/' + h.GetName() + '.png')
for h in [project_ecal, project_momatpca, project_momatbs, project_momatbs_combine, project_momatpca_combine, project_trkmode, project_trkmode_combine]:
	c.Clear()
	h.SetStats(0)
	h.Draw('hist')
	c.Print(path + '/' + h.GetName() + '.png')


# ratio

numerator = profile_momatpca.Clone() 
denominator = profile_momatbs.Clone()
numerator.Divide( denominator )
g_mompcaoverbs = histo_to_graph( numerator )
g_mompcaoverbs.SetMarkerStyle(24)
g_mompcaoverbs.SetMarkerColor(4)
g_mompcaoverbs.SetMarkerSize(1)

numerator = profile_momatpca_combine.Clone()
denominotor = profile_momatbs_combine.Clone()
numerator.Divide( denominotor )
g_combpcaoverbs = histo_to_graph( numerator )
g_combpcaoverbs.SetMarkerStyle(34)
g_combpcaoverbs.SetMarkerSize(1)
g_combpcaoverbs.SetMarkerColor(4)


# for a pretty plot
legend = ROOT.TLegend(.26,.4,.45,.58)
legend.SetTextSize(0.04)
legend.SetLineWidth(0)
legend.SetFillColor(0)
legend.SetBorderSize()
legend.AddEntry(gs[0],'ecal')
legend.AddEntry(gs[1],'trk. pca')
legend.AddEntry(gs[2],'trk. bs')
legend.AddEntry(gs[4],'comb. pca')
legend.AddEntry(gs[3],'comb. bs')

ctot = ROOT.TCanvas('ctot','',1000,1000)
ctot.cd()
c1 = ROOT.TPad('p1','',0, 0.5, 1, 1.0)
c1.Draw()
c1.cd()
c1.SetTopMargin(0.4)
c1.SetBottomMargin(0.01)
c1.SetRightMargin(0.2)
c1.SetLeftMargin(0.2)
c1.SetTicks(1,1)
gs[0].GetYaxis().SetTitle('elec. mom. resol.')
gs[0].GetYaxis().SetTitleSize(0.06)
gs[0].GetYaxis().SetTitleOffset(1.2)
gs[0].GetHistogram().SetMaximum( 0.15 )
gs[0].GetHistogram().SetMinimum( 0.0 )
gs[0].SetMarkerStyle(25)
gs[0].Draw('ap')
gs[0].GetXaxis().SetLabelSize(0)
gs[0].GetYaxis().SetLabelSize(0.05)
gs[0].GetYaxis().SetNdivisions(510)
gs[1].SetMarkerStyle(24)
gs[1].Draw('p')
gs[2].SetMarkerStyle(24)
gs[2].SetMarkerColor(2)
gs[2].Draw('p')
gs[3].SetMarkerStyle(34)
gs[3].Draw('p')
gs[4].SetMarkerStyle(34)
gs[4].SetMarkerColor(2)
gs[4].Draw('p')
gs[5].SetMarkerStyle(32)
gs[5].SetMarkerColor(2)
#gs[5].Draw('p')
gs[6].SetMarkerStyle(32)
gs[6].SetMarkerColor(1)
#gs[6].Draw('p')
legend.Draw('same')
ctot.cd()
c2 = ROOT.TPad('p2','', 0, 0.35, 1, 0.5)
c2.Draw()
c2.cd()
c2.SetTicks(1,1)
c2.SetTopMargin(0.01)
c2.SetBottomMargin(0.01)
c2.SetRightMargin(0.2)
c2.SetLeftMargin(0.2)
g_mompcaoverbs.Draw('ap')
g_mompcaoverbs.GetXaxis().SetLabelSize(0)
g_mompcaoverbs.GetYaxis().SetLabelSize(0.18)
g_mompcaoverbs.GetYaxis().SetNdivisions(503)
g_mompcaoverbs.GetHistogram().SetMaximum(1.6)
g_mompcaoverbs.GetHistogram().SetMinimum(0.5)
g_mompcaoverbs.GetYaxis().SetTitle('trk. pca/bs')
g_mompcaoverbs.GetYaxis().SetTitleSize(0.19)
g_mompcaoverbs.GetYaxis().SetTitleOffset(0.37)
l1 = ROOT.TLine(0,1,g_mompcaoverbs.GetXaxis().GetXmax(),1)
l2 = ROOT.TLine(0,1.3,g_mompcaoverbs.GetXaxis().GetXmax(),1.3)
l1.SetLineStyle(1)
l1.SetLineColor(2)
l1.SetLineWidth(2)
l2.SetLineStyle(2)
l2.SetLineColor(2)
l2.SetLineWidth(2)
l1.Draw('same')
l2.Draw('same')
ctot.cd()
c3 = ROOT.TPad('p3','', 0, 0.0, 1, 0.35)
c3.Draw()
c3.cd()
c3.SetTicks(1,1)
c3.SetTopMargin(0.01)
c3.SetBottomMargin(0.6)
c3.SetRightMargin(0.2) 
c3.SetLeftMargin(0.2)
g_combpcaoverbs.Draw('ap')
g_combpcaoverbs.GetYaxis().SetLabelSize(0.08)
g_combpcaoverbs.GetXaxis().SetLabelSize(0.08)
g_combpcaoverbs.GetYaxis().SetNdivisions(503)
g_combpcaoverbs.GetHistogram().SetMaximum(1.2)
g_combpcaoverbs.GetHistogram().SetMinimum(0.8)
g_combpcaoverbs.GetXaxis().CenterTitle(True)
g_combpcaoverbs.GetXaxis().SetTitleSize(0.08)
g_combpcaoverbs.GetYaxis().SetTitleSize(0.08)
g_combpcaoverbs.GetXaxis().SetTitle('P')
g_combpcaoverbs.GetYaxis().SetTitle('comb. pca/bs')
l1.Draw('same')
l3 = l2.Clone()
l3.SetY1(1.05)
l3.SetY2(1.05)
l3.Draw('same')
ctot.Print(path + '/mom_bsvspca.png')
ctot.Print(path + '/mom_bsvspca.pdf')
