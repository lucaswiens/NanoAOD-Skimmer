import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from leptonselector import leptonSel
import susyCutter

class susy1LeptonBaseProducer(Module):
	def __init__(self):
		pass
	def beginJob(self):
		pass
	def endJob(self):
		pass
	def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree, isMC=True):
		self.isMC = isMC
		self.cut = susyCutter.susyCutter()
		self.out = wrappedOutputTree
		self.out.branch("nGoodLepton", "I")
		self.out.branch("nVetoLepton", "I")
		self.out.branch("nAntiIsolatedGoodLepton", "I")
		self.out.branch("leadingLeptonSelected", "I")
		self.out.branch("secondLeptonSelected", "I")
		self.out.branch("nGoodElectron", "I")
		self.out.branch("nVetoElectron", "I")
		self.out.branch("nGoodMuon", "I")
		self.out.branch("nVetoMuon", "I")
		self.out.branch("leptonPt", "F")
		self.out.branch("leptonEta", "F")
		self.out.branch("leptonPhi", "F")
		self.out.branch("leptonPdgId", "I")
		self.out.branch("leptonRelIso", "F")
		self.out.branch("leptonMiniIso", "F")
		self.out.branch("leptonHOverE", "F")
		self.out.branch("leptonSF", "F")
		self.out.branch("MuonEffSF", "F")
		self.out.branch("ElectronEffSF", "F")
		self.out.branch("Lep2_pt", "F")


	def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
		pass
	def analyze(self, event):
		"""process event, return True (go to next module) or False (fail, go to next event)"""
		leptonSel(event)
		nLooseLepton = len(event.selectedLeptons)
		#Select good leptons (stricter than medium but not as strict as tight)
		selectedGoodLeptons, selectedVetoLeptons, antiIsolatedGoodLeptons = self.selectGoodLeptons(event.selectedLeptons, event.otherLeptons)
		nGoodLepton = len(selectedGoodLeptons)
		nVetoLepton = len(selectedVetoLeptons)
		nAntiIsolatedGoodLepton = len(antiIsolatedGoodLeptons)
		self.out.fillBranch("nGoodLepton", nGoodLepton)
		self.out.fillBranch("nVetoLepton", nVetoLepton)
		self.out.fillBranch("nAntiIsolatedGoodLepton", nAntiIsolatedGoodLepton)
		if nGoodLepton > 0:
			self.out.fillBranch("leadingLeptonSelected", 1)
			if nGoodLepton > 1:
				self.out.fillBranch("secondLeptonSelected", 1)
			else:
				self.out.fillBranch("secondLeptonSelected", 0)
		elif nAntiIsolatedGoodLepton == 0:
			self.out.fillBranch("leadingLeptonSelected", 0)
		else:
			self.out.fillBranch("leadingLeptonSelected", -1)

		# get number of tight el and mu
		goodElectron = [lep for lep in selectedGoodLeptons if abs(lep.pdgId) == 11]
		vetoElectron = [lep for lep in selectedVetoLeptons if abs(lep.pdgId) == 11]
		goodMuon = [lep for lep in selectedGoodLeptons if abs(lep.pdgId) == 13]
		vetoMuon = [lep for lep in selectedVetoLeptons if abs(lep.pdgId) == 13]
		nGoodElectron = len(goodElectron)
		nVetoElectron = len(vetoElectron)
		nGoodMuon = len(goodMuon)
		nVetoMuon = len(vetoMuon)
		self.out.fillBranch("nGoodElectron", nGoodElectron)
		self.out.fillBranch("nVetoElectron", nVetoElectron)
		self.out.fillBranch("nGoodMuon", nGoodMuon)
		self.out.fillBranch("nVetoMuon", nVetoMuon)
		""" TODO investigate electron energy calibration """
		## for El in tightEl:
		## 	# this branch is for investigating the electron energy calibraition
		## 	self.out.fillBranch("GoodEl_eCorr", El.eCorr )
		## 	self.h_eCorr_vs_eta_tight.Fill(El.eCorr, El.eta)
		## 	self.h_eCorr_vs_phi_tight.Fill(El.eCorr, El.phi)
		## 	self.h_eCorr_vs_pt_tight.Fill(El.eCorr, El.pt)
		## for El in VetoEl:
		## 	self.h_eCorr_vs_eta_veto.Fill(El.eCorr, El.eta)
		## 	self.h_eCorr_vs_phi_veto.Fill(El.eCorr, El.phi)
		## 	self.h_eCorr_vs_pt_veto.Fill(El.eCorr, El.pt)

		# Leading Lepton Variables
		if nGoodLepton > 0:
			self.out.fillBranch("leptonPt",selectedGoodLeptons[0].pt)
			self.out.fillBranch("leptonEta", selectedGoodLeptons[0].eta)
			self.out.fillBranch("leptonPhi", selectedGoodLeptons[0].phi)
			self.out.fillBranch("leptonPdgId", selectedGoodLeptons[0].pdgId)

			self.out.fillBranch("leptonRelIso", selectedGoodLeptons[0].pfRelIso03_all)
			self.out.fillBranch("leptonMiniIso", selectedGoodLeptons[0].miniPFRelIso_all)
			if hasattr(selectedGoodLeptons[0],"hoe"):
				self.out.fillBranch("leptonHOverE", selectedGoodLeptons[0].hoe)
			if self.isMC:
				for tlep in selectedGoodLeptons:
					if abs(tlep.pdgId) == 13:
						#sf_mu = self._worker_mu.getSF(tlep.pdgId,tlep.pt,tlep.eta)
						sf_mu = 1.0
						self.out.fillBranch("leptonSF", sf_mu)
						self.out.fillBranch("MuonEffSF", sf_mu)
					elif abs(tlep.pdgId) == 11:
						#sf_el = self._worker_el.getSF(tlep.pdgId,tlep.pt,tlep.eta)
						sf_el = 1.0
						self.out.fillBranch("leptonSF", sf_el)
						self.out.fillBranch("ElectronEffSF", sf_el)
					else:
						self.out.fillBranch("MuonEffSF", 1.0)
						self.out.fillBranch("ElectronEffSF", 1.0)
						self.out.fillBranch("leptonSF", 1.0)
			else:
				self.out.fillBranch("MuonEffSF", 1.0)
				self.out.fillBranch("ElectronEffSF", 1.0)
				self.out.fillBranch("leptonSF", 1.0)

		# save second leading lepton vars
		if len(selectedGoodLeptons) > 1:# 2nd tight lep
			self.out.fillBranch("Lep2_pt", selectedGoodLeptons[1].pt)

		# jets = Collection(event, "Jet")
		# jetPt = sum([jet.pt for jet in jets if jet.pt > 30])

		# MET = Object(event, "MET")
		# genMET = Object(event, "GenMET")
		# self.out.fillBranch("nLepton", nElectrons + nMuons)
		# self.out.fillBranch("LT", electronPt + muonPt)
		# self.out.fillBranch("HT", jetPt)
		return True
	def selectGoodLeptons(self, selectedLeptons, otherLeptons):
		selectedGoodLeptons = [] #Not really tight but tighter than medium WP
		selectedVetoLeptons = []
		antiIsolatedGoodLeptons = []
		for lepton in selectedLeptons:
			leptonEta = abs(lepton.eta)
			# Pt cut
			if lepton.pt < 10: continue

			# Iso cut -- to be compatible with the trigger
			if lepton.miniPFRelIso_all > self.cut.trig_miniIsoCut: continue
			# MUONS
			if(abs(lepton.pdgId) == 13):
				if leptonEta > self.cut.muonEta: continue
				#passId = False
				passId = lepton.mediumId
				passIso = lepton.miniPFRelIso_all < self.cut.muo_miniIsoCut
				passIP = lepton.sip3d < self.cut.goodMu_sip3d

				# selected muons
				if passId and passIP and passIso:
					selectedGoodLeptons.append(lepton)
				else:
					selectedVetoLeptons.append(lepton)
				if not passIso:
					antiIsolatedGoodLeptons.append(lepton)

			# ELECTRONS
			elif(abs(lepton.pdgId) == 11):

				if leptonEta > self.cut.electronEta: continue

				passIso = False
				passGoodID = lepton.cutBased == 4
				passMediumID = lepton.cutBased >= 3
				passVetoID = lepton.cutBased >= 1
				# selected
				if passGoodID:
					passIso = lepton.miniPFRelIso_all < self.cut.ele_miniIsoCut
					if passIso:
						selectedGoodLeptons.append(lepton)
					else:
						selectedVetoLeptons.append(lepton)

				# anti-selected
				elif not passMediumID:
					# all anti leptons are veto for selected
					selectedVetoLeptons.append(lepton)

					# Iso check
					passIso = lepton.miniPFRelIso_all < self.cut.Lep_miniIsoCut
					# other checks
					passOther = False
					if hasattr(lepton,"hoe"):
						passOther = lepton.hoe > 0.01

					# fill
					if passIso and passOther:
						antiIsolatedGoodLeptons.append(lepton)
				# Veto leptons
				elif passVetoID:
					selectedVetoLeptons.append(lepton)

		for lepton in otherLeptons:
			# check acceptance
			leptonEta = abs(lepton.eta)
			if leptonEta > self.cut.electronEta: continue
			# Pt cut
			if lepton.pt < 10: continue #TODO define this later in cut
			# Iso cut -- to be compatible with the trigger
			if lepton.miniPFRelIso_all > self.cut.trig_miniIsoCut: continue

			# Muons
			if(abs(lepton.pdgId) == 13):
				passIso = lepton.miniPFRelIso_all > self.cut.muo_miniIsoCut
				if passIso:
					antiIsolatedGoodLeptons.append(lepton)

			# Electrons
			elif(abs(lepton.pdgId) == 11):
				if leptonEta > self.cut.electronEta: continue
				if lepton.miniPFRelIso_all > self.cut.Lep_miniIsoCut: continue

				passMediumID = lepton.cutBased >= 3
				passVetoID = lepton.cutBased >= 1
				if not passMediumID:
					passOther = False
					if hasattr(lepton,"hoe"):
						passOther = lepton.hoe > 0.01
					if passOther:
						antiIsolatedGoodLeptons.append(lepton)
		return selectedGoodLeptons, selectedVetoLeptons, antiIsolatedGoodLeptons


# define modules using the syntax 'name = lambda: constructor' to avoid having them loaded when not needed
susy1LeptonProducer = lambda : susy1LeptonBaseProducer()

