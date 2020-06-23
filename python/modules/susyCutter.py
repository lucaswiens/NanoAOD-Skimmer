class susyCutter():
	def __init__(self):
		""" Eta requirement """
		self.centralEta = 2.4
		self.electronEta = 2.4

		""" Jets """
		self.smearJER = True

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
		self.goodEl_lostHits = 0
		self.goodEl_sip3d = 4
		self.goodMu_sip3d = 4
