class susyCutter():
	def __init__(self, era = 2016):
		""" Eta requirement """
		self.centralEta = 2.4
		self.electronEta = 2.4

		""" Jets """
		self.jetEta = 2.5

		""" MUONS """
		self.muID = 'medium' # 'medium'(2015) or 'ICHEPmediumMuonId' (2016)
		self.muonEta = 2.4

		""" Electrons """
		self.eleID = 'CB'

		""" Isolation """
		self.ele_miniIsoCut = 0.1
		self.muo_miniIsoCut = 0.2
		self.Lep_miniIsoCut = 0.4
		self.trig_miniIsoCut = 0.8

		""" Lepton cuts (for MVAID) """
		self.minLeptonPt = 10.
		self.goodEl_lostHits = 0
		self.goodEl_sip3d = 4
		self.goodMu_sip3d = 4


		self.hadMT2cut = 60
		self.lepMT2cut = 80

		if era == 2016:
			self.btagLooseWP = 0.5426
			self.btagMediumWP = 0.8484
			self.btagTightWP = 0.9535
			#https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
			self.btagDeepLooseWP = 0.2219
			self.btagDeepMediumWP = 0.6324
			self.btagDeepTightWP = 0.8958
		elif era == 2017:
			self.btagLooseWP = 0.5803
			self.btagMediumWP = 0.8838
			self.btagTightWP = 0.9693
			# DeepCSV (new Deep Flavour tagger)
			#https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
			self.btagDeepLooseWP = 0.1522
			self.btagDeepMediumWP = 0.4941
			self.btagDeepTightWP = 0.8001
