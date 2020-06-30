import ROOT
import array
import os

ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import matchObjectCollection, matchObjectCollectionMultiple

from leptonselector import leptonSel
import susyCutter

class susy1LeptonBaseProducer(Module):
	def __init__(self, isMC=True , isSignal=False, era=2016, muonSelectionTag = "MediumWP_2016", electronSelectionTag = "TightWP_2016"):
		self.isMC = isMC
		self.isSignal = isSignal
		self.era = era
		self.mt2obj = ROOT.heppy.Davismt2.Davismt2()
		self.muonSelectionTag = muonSelectionTag
		self.electronSelectionTag = electronSelectionTag
		muonFile, electronFile = [], [] #For some reason python complains about electronFile being referenced before its assignment unless it is declared here..

		if isMC:
			if self.era == 2016:
				if self.muonSelectionTag=="LooseWP_2016":
					muonFile = ["Mu_Trg.root","Mu_ID.root","Mu_Iso.root"]
					muonHistogram = ["IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio",
							"MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio",
							"LooseISO_LooseID_pt_eta/pt_abseta_ratio"]
				if self.electronSelectionTag=="GPMVA90_2016":
					electronFile = ["EGM2D_eleGSF.root","EGM2D_eleMVA90.root"]
					electronHistogram = ["EGamma_SF2D", "EGamma_SF2D"]
				if self.muonSelectionTag=="MediumWP_2016":
					muonFile = ["Mu_Trg.root","Mu_ID.root","Mu_Iso.root"]
					muonHistogram = ["IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio",
							"MC_NUM_MediumID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio",
							"LooseISO_MediumID_pt_eta/pt_abseta_ratio"]
				if self.electronSelectionTag=="Tight_2016":
					electronFile = ["EGM2D_eleGSF.root","EGM2D_eleMVA90.root"]
					electronHistogram = ["EGamma_SF2D", "EGamma_SF2D"]
			elif self.era == 2017:
				if self.muonSelectionTag=="LooseWP_2017":
					muonFile=["Mu_Trg17.root","Mu_ID17.root","Mu_Iso17.root"]
					muonHistogram = ["IsoMu27_PtEtaBins/pt_abseta_ratio",
							"NUM_LooseID_DEN_genTracks_pt_abseta",
							"NUM_LooseRelIso_DEN_LooseID_pt_abseta"]
				if self.electronSelectionTag=="GPMVA90_2017":
					electronFile = ["EGM2D_eleGSF17.root","EGM2D_eleMVA90_17.root"]
					electronHistogram = ["EGamma_SF2D", "EGamma_SF2D"]
				if self.muonSelectionTag=="MediumWP_2017":
					muonFile=["Mu_Trg17.root","Mu_ID17.root","Mu_Iso17.root"]
					muonHistogram = ["IsoMu27_PtEtaBins/pt_abseta_ratio",
							"NUM_MediumID_DEN_genTracks_pt_abseta",
							"NUM_LooseRelIso_DEN_MediumID_pt_abseta"]
				if self.electronSelectionTag=="Tight_2017":
					electronFile = ["EGM2D_eleGSF17.root","EGM2D_eleMVA90_17.root"]
					electronHistogram = ["EGamma_SF2D", "EGamma_SF2D"]
			else :
				raise ValueError("ERROR: Invalid era = '%i'!" % self.era)

		if self.isMC:
			muonFile = ["%s/src/PhysicsTools/NanoAODTools/python/postprocessing/data/leptonSF/" % os.environ['CMSSW_BASE'] + f for f in muonFile]
			electronFile = ["%s/src/PhysicsTools/NanoAODTools/python/postprocessing/data/leptonSF/" % os.environ['CMSSW_BASE'] + f for f in electronFile]
			self.muonFile = ROOT.std.vector(str)(len(muonFile))
			self.muonHistogram = ROOT.std.vector(str)(len(muonFile))
			for i in range(len(muonFile)):
				self.muonFile[i] = muonFile[i]
				self.muonHistogram[i] = muonHistogram[i];
			self.electronFile = ROOT.std.vector(str)(len(electronFile))
			self.electronHistogram = ROOT.std.vector(str)(len(electronFile))
			for i in range(len(electronFile)):
				self.electronFile[i] = electronFile[i]
				self.electronHistogram[i] = electronHistogram[i];

		if "/LeptonEfficiencyCorrector_cc.so" not in ROOT.gSystem.GetLibraries():
			print "Load C++ Worker"
			ROOT.gROOT.ProcessLine(".L %s/src/PhysicsTools/NanoAODTools/python/postprocessing/helpers/LeptonEfficiencyCorrector.cc+" % os.environ['CMSSW_BASE'])
		pass

	def beginJob(self):
		if self.isMC :
			self.workerMuon = ROOT.LeptonEfficiencyCorrector(self.muonFile, self.muonHistogram)
			self.workerElectron = ROOT.LeptonEfficiencyCorrector(self.electronFile, self.electronHistogram)
		pass

	def endJob(self):
		pass

	def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree, isMC=True):
		self.isMC = isMC
		self.cut = susyCutter.susyCutter()
		self.out = wrappedOutputTree
		self.filename = inputFile.GetName()

		self.out.branch("nGoodLepton", "I")
		self.out.branch("nVetoLepton", "I")
		self.out.branch("nAntiIsolatedGoodLepton", "I")
		self.out.branch("leadingLeptonSelected", "I")
		self.out.branch("secondLeptonSelected", "I")
		self.out.branch("nGoodElectron", "I")
		self.out.branch("nVetoElectron", "I")
		self.out.branch("nGoodMuon", "I")
		self.out.branch("nVetoMuon", "I")
		self.out.branch("nJet30", "I")
		self.out.branch("nCleanJet30", "I")
		self.out.branch("nJet40", "I")
		self.out.branch("leptonPdgId", "I")
		self.out.branch("isDPhiSignal", "I")
		self.out.branch("Flag_fastSimCorridorJetCleaning", "I")
		self.out.branch("nJet", "I")
		self.out.branch("LSLjetptGT80", "I")
		self.out.branch("nBJet", "I")
		self.out.branch("nBJets30", "I")
		self.out.branch("nBJetDeep", "I")
		self.out.branch("nBJets40", "I")
		self.out.branch("LSLjetptGT80", "I")

		self.out.branch("leptonPt", "F")
		self.out.branch("leptonEta", "F")
		self.out.branch("leptonPhi", "F")
		self.out.branch("leptonRelIso", "F")
		self.out.branch("leptonMiniIso", "F")
		self.out.branch("leptonHOverE", "F")
		self.out.branch("leptonSF", "F")
		self.out.branch("MuonEffSF", "F")
		self.out.branch("ElectronEffSF", "F")
		self.out.branch("secondLeptonPt", "F")
		self.out.branch("goodElectronECorr", "F")
		self.out.branch("vetoElectronECorr", "F")
		self.out.branch("jetPt", "F")
		self.out.branch("jetEta", "F")
		self.out.branch("jetPt_2", "F")
		self.out.branch("jetEta_2", "F")
		self.out.branch("HT", "F")
		self.out.branch("jetHt30", "F")
		self.out.branch("cleanJetHt30", "F")
		self.out.branch("centralJetHt40", "F")
		self.out.branch("cleanJetHt30", "F")

		self.out.branch("DeltaPhiLepW", "F")
		self.out.branch("dPhi", "F")
		self.out.branch("LT", "F")
		self.out.branch("Lp", "F")
		self.out.branch("MT", "F")
		self.out.branch("GendPhi", "F")
		self.out.branch("GenLT", "F")
		self.out.branch("GenMET", "F")
		self.out.branch("isSignalRegion", "I")
		self.out.branch("RA2_muJetFilter", "B") # normally true for now # don't know how to get the Muon energy fraction from EMEF
		self.out.branch("dileptonMass", "F")
		self.out.branch("minMWjj", "F")
		self.out.branch("minMWjjPt", "F")
		self.out.branch("bestMWjj", "F")
		self.out.branch("bestMWjjPt", "F")
		self.out.branch("bestMTopHad", "F")
		self.out.branch("bestMTopHadPt", "F")
		self.out.branch("iso_MT2", "F")
		self.out.branch("iso_pt", "F")
		self.out.branch("iso_had", "I")
		self.out.branch("iso_had", "I")
		self.out.branch("iso_Veto", "B")
		self.out.branch("PD_JetHT", "B")
		self.out.branch("PD_SingleElectron", "B")
		self.out.branch("PD_SingleMuon", "B")
		self.out.branch("PD_MET", "B")



	def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
		pass

	def analyze(self, event):
		"""process event, return True (go to next module) or False (fail, go to next event)"""
		# Leptons
		leptonSel(event)
		nLooseLepton = len(event.selectedLeptons)
		#Select good leptons (stricter than medium but not as strict as tight)
		selectedGoodLeptons = []
		selectedVetoLeptons = []
		antiIsolatedGoodLeptons = []
		for lepton in event.selectedLeptons:
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
					if hasattr(lepton, "hoe"):
						passOther = lepton.hoe > 0.01

					# fill
					if passIso and passOther:
						antiIsolatedGoodLeptons.append(lepton)
				# Veto leptons
				elif passVetoID:
					selectedVetoLeptons.append(lepton)

		for lepton in event.otherLeptons:
			# check acceptance
			leptonEta = abs(lepton.eta)
			if leptonEta > self.cut.electronEta: continue
			# Pt cut
			if lepton.pt < self.cut.minLeptonPt: continue
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
					if hasattr(lepton, "hoe"):
						passOther = lepton.hoe > 0.01
					if passOther:
						antiIsolatedGoodLeptons.append(lepton)
		#selectedGoodLeptons, selectedVetoLeptons, antiIsolatedGoodLeptons = self.selectGoodLeptons(event.selectedLeptons, event.otherLeptons)
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
		for electron in goodElectron:
			self.out.fillBranch("goodElectronECorr", electron.eCorr)
		for electron in vetoElectron:
			self.out.fillBranch("vetoElectronECorr", electron.eCorr)

		if nGoodLepton > 0:
			self.out.fillBranch("leptonPt", selectedGoodLeptons[0].pt)
			self.out.fillBranch("leptonEta", selectedGoodLeptons[0].eta)
			self.out.fillBranch("leptonPhi", selectedGoodLeptons[0].phi)
			self.out.fillBranch("leptonPdgId", selectedGoodLeptons[0].pdgId)

			self.out.fillBranch("leptonRelIso", selectedGoodLeptons[0].pfRelIso03_all)
			self.out.fillBranch("leptonMiniIso", selectedGoodLeptons[0].miniPFRelIso_all)
			if hasattr(selectedGoodLeptons[0], "hoe"):
				self.out.fillBranch("leptonHOverE", selectedGoodLeptons[0].hoe)
			if self.isMC:
				for lepton in selectedGoodLeptons:
					if abs(lepton.pdgId) == 13:
						muonScaleFactor = self.workerMuon.getSF(lepton.pdgId, lepton.pt, lepton.eta)
						self.out.fillBranch("leptonSF", muonScaleFactor)
						self.out.fillBranch("MuonEffSF", muonScaleFactor)
					elif abs(lepton.pdgId) == 11:
						electronScaleFactor = self.workerElectron.getSF(lepton.pdgId, lepton.pt, lepton.eta)
						self.out.fillBranch("leptonSF", electronScaleFactor)
						self.out.fillBranch("ElectronEffSF", electronScaleFactor)
					else:
						self.out.fillBranch("MuonEffSF", 1.0)
						self.out.fillBranch("ElectronEffSF", 1.0)
						self.out.fillBranch("leptonSF", 1.0)
			else:
				self.out.fillBranch("MuonEffSF", 1.0)
				self.out.fillBranch("ElectronEffSF", 1.0)
				self.out.fillBranch("leptonSF", 1.0)

		if len(selectedGoodLeptons) > 1:
			self.out.fillBranch("secondLeptonPt", selectedGoodLeptons[1].pt)

		# Jets
		Jets = Collection(event, "Jet")
		jets = [j for j in Jets if j.pt > 20 and abs(j.eta) < 2.4]

		centralJets30 = []
		centralJets40 = []
		cleanJets25 = []
		cleanJets = []
		# fill this flag but defaults to 1 and then change it after the proper selection
		self.out.fillBranch("Flag_fastSimCorridorJetCleaning", 1)

		for index, jet in enumerate(jets):
			jet.pt = event.Jet_pt_nom[index]
			if self.isSignal: #only check for signals (see condition check above)
				self.out.fillBranch("isDPhiSignal",1)
				genJets = Collection(event, "GenJet" )
				pairs = matchObjectCollection(Jets, genJets)
				genJet = pairs[jet]
				if genJet is not None :
					if jet.pt > 20 and abs(jet.eta) < self.cut.jetEta and (genJet.pt == 0) and jet.chHEF < 0.1:
						self.out.fillBranch("Flag_fastSimCorridorJetCleaning", 0 )
			if jet.pt > 25 :
				cleanJets25.append(jet)
			if jet.pt  >  30 and abs(jet.eta) < self.cut.centralEta:
				centralJets30.append(jet)
			if jet.pt > 40 and abs(jet.eta) < self.cut.centralEta:
				centralJets40.append(jet)

		nJet = len(jets)
		nCentralJets30 = len(centralJets30)
		nCentralJets40 = len(centralJets40)
		self.out.fillBranch("nJet", nJet)
		self.out.fillBranch("nJet30", nCentralJets30)
		self.out.fillBranch("nJet40", nCentralJets40)

		##############################
		## Local cleaning from leptons
		##############################
		dRminCut = 0.4 #TODO susy cutter

		# Do cleaning a la CMG: clean max 1 jet for each lepton (the nearest)
		cleanJets30 = centralJets30
		if len(selectedGoodLeptons) > 0:
			tightLeptons = selectedGoodLeptons
		elif len(antiIsolatedGoodLeptons) > 0:
			tightLeptons = antiIsolatedGoodLeptons
		else:
			tightLeptons = []
		for lepton in tightLeptons:
			closestJet = None
			dRmin = 99
			for jet in centralJets30:
				dR = jet.p4().DeltaR(lepton.p4())
				if dR < dRmin:
					closestJet = jet
					dRmin = dR
			if dRmin < dRminCut:
				cleanJets30.remove(closestJet)
			dR = 99 #TODO think if this makes sense, maybe I just remove ones that are closer thant he closest
			for jet25 in cleanJets25:
				dR = jet25.p4().DeltaR(lepton.p4())
				if dR < dRmin:
					cleanJets.append(jet25)
		# cleaned jets
		nCleanJet30 = len(cleanJets30)

		self.out.fillBranch("nCleanJet30", nCleanJet30)

		if nCleanJet30 > 0:
			self.out.fillBranch("jetPt", cleanJets30[0].pt)
			self.out.fillBranch("jetEta", cleanJets30[0].eta)
		if nCleanJet30 > 1:
			self.out.fillBranch("jetPt_2", cleanJets30[1].pt)
			self.out.fillBranch("jetEta_2", cleanJets30[1].eta)

		# imho, use Jet2_pt > 80 instead TODO why?
		self.out.fillBranch("LSLjetptGT80", 1 if sum([jet.pt > 80 for jet in cleanJets30]) >= 2 else 0)

		self.out.fillBranch("HT", sum([j.pt for j in cleanJets30]))
		self.out.fillBranch("cleanJetHt30", sum([j.pt for j in cleanJets30]))
		self.out.fillBranch("jetHt30", sum([j.pt for j in jets if j.pt>30]))
		#centralJetHt40 = sum([jet.pt for jet in centralJets40])
		self.out.fillBranch("centralJetHt40", sum([j.pt for j in centralJets40]))
		#self.out.fillBranch("centralJetHt40", centralJetHt40)

		## B tagging WPs for CSVv2 (CSV-IVF)
		## from: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagging#Preliminary_working_or_operating
		# WP defined on top
		#TODO define in init function, possibly use susyCutter
		btagWP = self.cut.btagMediumWP
		btagDeepMediumWP = self.cut.btagDeepMediumWP
		btagLooseWP = self.cut.btagLooseWP

		BJetMedium30 = []
		BJetMedium40 = []

		nBJetDeep = 0

		for jet in cleanJets30:
			if jet.btagCSVV2 > btagWP:
				BJetMedium30.append(jet)
			if jet.btagDeepB > btagDeepMediumWP:
				nBJetDeep += 1

		for jet in centralJets40:
			if jet.btagCSVV2 > btagWP:
				BJetMedium40.append(jet)

		# using cleaned collection!
		self.out.fillBranch("nBJet", len(BJetMedium30))
		self.out.fillBranch("nBJets30", len(BJetMedium30))
		self.out.fillBranch("nBJetDeep", nBJetDeep)
		# using normal collection
		self.out.fillBranch("nBJets40", len(BJetMedium40))
		#TODO store btag maybe?

		# High Level Variables
		dPhiLeptonW = -999
		GendPhiLeptonW = -999
		LT = -999
		GenLT = -999
		Lp = -999
		MT = -999
		metP4 = ROOT.TLorentzVector(0,0,0,0)
		metP4.SetPtEtaPhiM(event.MET_pt_nom, 0., event.MET_phi_nom, 0)

		if len(tightLeptons) > 0:
			#GenrecoWp4 =  tightLeptons[0].p4() + GenmetP4 #this does not make sense, use genmet of leptons
			#GendPhiLeptonW = tightLeptons[0].p4().DeltaPhi(GenrecoWp4)
			#GenLT = tightLeptons[0].pt + GenmetP4.Pt()
			recoWp4 =  tightLeptons[0].p4() + metP4
			dPhiLeptonW = tightLeptons[0].p4().DeltaPhi(recoWp4)
			LT = tightLeptons[0].pt + metP4.Pt()
			Lp = tightLeptons[0].pt / recoWp4.Pt() * ROOT.TMath.Cos(dPhiLeptonW)
			MT = ROOT.TMath.Sqrt(2*metP4.Pt()*tightLeptons[0].pt * (1-ROOT.TMath.Cos(dPhiLeptonW)))

		dPhi = abs(dPhiLeptonW)
		self.out.fillBranch("DeltaPhiLepW", dPhiLeptonW)
		self.out.fillBranch("dPhi", dPhi)
		#self.out.fillBranch("ST", LT) #ST? What?
		self.out.fillBranch("LT", LT)
		self.out.fillBranch("Lp", Lp)
		self.out.fillBranch("MT", MT)
		#self.out.fillBranch("GendPhi", abs(GendPhiLeptonW))
		#self.out.fillBranch("GenLT", GenLT)
		#self.out.fillBranch("GenMET", GenmetP4.Pt())

		# SIGNAL REGION FLAG
		# isSignalRegion SR vs CR flag
		isSignalRegion = 0

		# 0-B SRs -- simplified dPhi
		if len(BJetMedium30) == 0:# check the no. of Bjets
			if LT < 250:   isSignalRegion = 0
			elif LT > 250: isSignalRegion = dPhi > 0.75
			# BLIND data
			if (not self.isMC)  and nCleanJet30 >= 5:
				isSignalRegion = - isSignalRegion
		# Multi-B SRs
		elif nCleanJet30 < 99:
			if LT < 250:   isSignalRegion = 0
			elif LT < 350: isSignalRegion = dPhi > 1.0
			elif LT < 600: isSignalRegion = dPhi > 0.75
			elif LT > 600: isSignalRegion = dPhi > 0.5
			# BLIND data
			if (not self.isMC) and nCleanJet30 >= 6:
				isSignalRegion = - isSignalRegion

		self.out.fillBranch("isSignalRegion", isSignalRegion)

		#############
		## Playground
		#############

		# di-lepton mass: opposite-sign, same flavour
		dileptonMass = 0

		if len(tightLeptons) > 1:

			leadingLepton = tightLeptons[0]
			id1 = leadingLepton.pdgId

			for trailingLepton in event.selectedLeptons[1:]:
				# TODO check if this could be tgrue for more than one pair of leptons
				if trailingLepton.pdgId + leadingLepton.pdgId == 0:
					dilepP4 = leadingLepton.p4() + trailingLepton.p4()
					dileptonMass = dilepP4.M()

		self.out.fillBranch("dileptonMass", dileptonMass)

		# RA2 proposed filter
		# TODO Learn what this is
		self.out.fillBranch("RA2_muJetFilter", True) # normally true for now # don't know how to get the Muon energy fraction from EMEF
		#for j in cJet30Clean:
		#	if j.pt > 200 and j.chEmEF > 0.5 and abs(ROOT.TMath.ACos(ROOT.TMath.Cos(j.phi-metP4.Phi()))) > (ROOT.TMath.pi - 0.4):
		#		self.out.fillBranch("RA2_muJetFilter", False)


		## MET FILTERS for data looks like the met filters are applied already for nanoAOD #TODO check how this is done, in which module
		# Top Tagging
		lightJets = [ j for j in cleanJets if not j.btagCSVV2 == btagWP ]
		bjetsLoose  = [ j for j in cleanJets if j.btagCSVV2== self.cut.btagLooseWP]
		minMWjj   = -999
		minMWjjPt = -999
		bestMWjj   = -999
		bestMWjjPt = -999
		bestMTopHad   = -999
		bestMTopHadPt = -999
		for i1, j1 in enumerate(lightJets):
			for i2 in xrange(i1+1, len(lightJets)):
				j2 = lightJets[i2]
				jjp4 = j1.p4() + j2.p4()
				mjj  = jjp4.M()
				if mjj > 30 and mjj < minMWjj:
					minMWjj = mjj
					minMWjjPt = jjp4.Pt()
					#self.out.fillBranch("minMWjj", minMWjj) #probably do not fill in here
					#self.out.fillBranch("minMWjjPt", minMWjjPt)
				if abs(mjj-80.4) < abs(bestMWjj-80.4):
					bestMWjj = mjj
					bestMWjjPt = jjp4.Pt()
					#self.out.fillBranch("bestMWjj", bestMWjj)
					#self.out.fillBranch("bestMWjjPt", bestMWjjPt)
					for bj in bjetsLoose:
						if deltaR(bj.eta(), bj.phi(), j1.eta(), j1.phi()) < 0.1 or deltaR(bj.eta(), bj.phi(), j2.eta(), j2.phi()) < 0.1: continue
						tp4 = jjp4 + bj.p4()
						mtop = tp4.M()
						if abs(mtop-172) < abs(bestMTopHad - 172):
							bestMTopHad = mtop
							bestMTopHadPt = tp4.Pt()
							#self.out.fillBranch("bestMTopHad", bestMTopHad)
							#self.out.fillBranch("bestMTopHadPt", bestMTopHadPt)

		self.out.fillBranch("minMWjj", minMWjj)
		self.out.fillBranch("minMWjjPt", minMWjjPt)
		self.out.fillBranch("bestMWjj", bestMWjj)
		self.out.fillBranch("bestMWjjPt", bestMWjjPt)
		self.out.fillBranch("bestMTopHad", bestMTopHad)
		self.out.fillBranch("bestMTopHadPt", bestMTopHadPt)
		# isolated tracks after basic selection (((pt>5 && (abs(pdgId) == 11 || abs(pdgId) == 13)) || pt > 10) && (abs(pdgId) < 15 || abs(eta) < self.cut.jetEta) && abs(dxy) < 0.2 && abs(dz) < 0.1 && ((pfIsolationDR03().chargedHadronIso < 5 && pt < 25) || pfIsolationDR03().chargedHadronIso/pt < 0.2)) and lepton veto
		# First check is the event has the IsoTrack or not
		if hasattr(event, "nIsoTrack"):
			tracks = [j for j in Collection(event, "IsoTrack", "nIsoTrack")]
			trackp4 = ROOT.TLorentzVector(0, 0, 0, 0)

			# min dR between good lep and iso track
			minDR = 0.1
			# MT2 cuts for hadronic and leptonic veto tracks
			#hadMT2cut = 60
			#lepMT2cut = 80
			if (len(tightLeptons) > 0) and len(tracks) > 0:
				#for i, t in enumerate(tracks):
				for track in tracks:
					# looking for opposite charged tracks
					#if tightLeptons[0].charge == track.charge: continue # not track charge is founded replace with the next copule of lines
					if track.isHighPurityTrack == False : continue
					#print track.miniPFRelIso_chg
					# not track mass is founded
					trackp4.SetPtEtaPhiM(track.pt, track.eta, track.phi,0.)
					if trackp4.DeltaR(tightLeptons[0].p4()) < minDR:
						p1 = tightLeptons[0].p4()
						p2 = trackp4
						a = array.array('d', [p1.M(), p1.Px(), p1.Py()])
						b = array.array('d', [p2.M(), p2.Px(), p2.Py()])
						c = array.array('d', [metP4.M(), metP4.Px(), metP4.Py()])
						self.mt2obj.set_momenta(a, b, c)
						self.mt2obj.set_mn(0)
						self.out.fillBranch("iso_MT2", self.mt2obj.get_mt2())
						self.out.fillBranch("iso_pt", p2.Pt())
					# cuts on MT2 as defined above
					if abs(track.pdgId)>10 and abs(track.pdgId)<14:
						self.out.fillBranch("iso_had", 0)  #leptonic
						cut = self.cut.lepMT2cut
					else:
						self.out.fillBranch("iso_had", 1)  #hadronic track
						cut = self.cut.hadMT2cut
					if self.mt2obj.get_mt2() <= cut:
						self.out.fillBranch("iso_Veto", True)

		#self.out.fillBranch("Xsec", self.xs)
		if 'JetHT' in self.filename:
			self.out.fillBranch("PD_JetHT", True)
		else:
			self.out.fillBranch("PD_JetHT", False)
		if 'SingleEle' in self.filename:
			self.out.fillBranch("PD_SingleElectron", True)
		else:
			self.out.fillBranch("PD_SingleElectron", False)
		if 'SingleMu' in self.filename:
			self.out.fillBranch("PD_SingleMuon", True)
		else:
			self.out.fillBranch("PD_SingleMuon", False)
		if 'MET' in self.filename:
			self.out.fillBranch("PD_MET", True)
		else:
			self.out.fillBranch("PD_MET", False)
		return True

# define modules using the syntax 'name = lambda: constructor' to avoid having them loaded when not needed
susy1LeptonBase2016 = lambda : susy1LeptonBaseProducer(isMC = True, isSignal = False, era = 2016, muonSelectionTag = "MediumWP_2016", electronSelectionTag = "TightWP_2016")
susy1LeptonBase2017 = lambda : susy1LeptonBaseProducer(isMC = True, isSignal = False, era = 2017, muonSelectionTag = "MediumWP_2017", electronSelectionTag = "TightWP_2017")
susy1LeptonBaseSignal2016 = lambda : susy1LeptonBaseProducer(isMC = True, isSignal = True, era = 2016, muonSelectionTag = "MediumWP_2016", electronSelectionTag = "TightWP_2016")
susy1LeptonBaseSignal2017 = lambda : susy1LeptonBaseProducer(isMC = True, isSignal = True, era = 2017, muonSelectionTag = "MediumWP_2017", electronSelectionTag = "TightWP_2017")
susy1LeptonBaseData2016 = lambda : susy1LeptonBaseProducer(isMC = False, isSignal = False, era = 2016, muonSelectionTag = "MediumWP_2016", electronSelectionTag = "TightWP_2016")
susy1LeptonBaseData2017 = lambda : susy1LeptonBaseProducer(isMC = False, isSignal = False, era = 2017, muonSelectionTag = "MediumWP_2017", electronSelectionTag = "TightWP_2017")
