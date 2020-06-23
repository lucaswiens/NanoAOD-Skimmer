# Copyright (C) 2014 Colin Bernet
# https://github.com/cbernet/heppy/blob/master/LICENSE
import ROOT
import math, os
import array
import operator

ROOT.PyConfig.IgnoreCommandLineOptions = True
# for met object
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import matchObjectCollection, matchObjectCollectionMultiple

import itertools
from ROOT import TLorentzVector, TVector2, std
mt2obj = ROOT.heppy.Davismt2.Davismt2()
from deltar import bestMatch



def leptonSel( event):
	event.allelectrons = Collection(event, "Electron")
	allmuons = Collection(event, "Muon")
	event.inclusiveLeptons = []
	event.selectedLeptons = []
	event.selectedMuons = []
	event.selectedElectrons = []
	event.otherLeptons = []
	inclusiveMuons = []
	inclusiveElectrons = []
	for mu in allmuons:
		if (mu.pt>10 and abs(mu.eta)<2.4 and
				abs(mu.dxy)<0.5 and abs(mu.dz)<1.):
			inclusiveMuons.append(mu)
	for ele in event.allelectrons:
		if ( ele.cutBased >=1 and
				ele.pt>10 and abs(ele.eta)<2.4 ):#and abs(ele.dxy)<0.5 and abs(ele.dz)<1. and ele.lostHits <=1.0) :
			inclusiveElectrons.append(ele)
	event.inclusiveLeptons = inclusiveMuons + inclusiveElectrons
	# make loose leptons (basic selection)
	for mu in inclusiveMuons :
			if (mu.pt > 10 and abs(mu.eta) < 2.4 and mu.miniPFRelIso_all < 0.4 and mu.isPFcand and abs(mu.dxy)<0.05 and abs(mu.dz)<1):
				event.selectedLeptons.append(mu)
				event.selectedMuons.append(mu)
			else:
				event.otherLeptons.append(mu)
	looseMuons = event.selectedLeptons[:]
	for ele in inclusiveElectrons :
		ele.looseIdOnly = ele.cutBased >=1
		if (ele.looseIdOnly and
					ele.pt>10 and abs(ele.eta)<2.4 and ele.miniPFRelIso_all < 0.4 and ele.isPFcand and #ele.convVeto and abs(ele.dxy)<0.05 and abs(ele.dz)<0.5  and ele.lostHits <=1.0 and
					(bestMatch(ele, looseMuons)[1] > (0.05**2))):
				event.selectedLeptons.append(ele)
				event.selectedElectrons.append(ele)
		else:
				event.otherLeptons.append(ele)
	event.otherLeptons.sort(key = lambda l : l.pt, reverse = True)
	event.selectedLeptons.sort(key = lambda l : l.pt, reverse = True)
	event.selectedMuons.sort(key = lambda l : l.pt, reverse = True)
	event.selectedElectrons.sort(key = lambda l : l.pt, reverse = True)
	event.inclusiveLeptons.sort(key = lambda l : l.pt, reverse = True)
	return event.otherLeptons , event.selectedLeptons, event.selectedMuons , event.selectedElectrons ,event.inclusiveLeptons , event.allelectrons
