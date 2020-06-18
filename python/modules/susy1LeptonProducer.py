import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class susy1LeptonBaseProducer(Module):
	def __init__(self):
		pass
	def beginJob(self):
		pass
	def endJob(self):
		pass
	def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
		self.out = wrappedOutputTree
		self.out.branch("nMuons",  "F");
		self.out.branch("nElectrons",  "F");
		self.out.branch("nLeptons",  "F");
	def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
		pass
	def analyze(self, event):
		"""process event, return True (go to next module) or False (fail, go to next event)"""
		electrons = Collection(event, "Electron")
		muons = Collection(event, "Muon")
		#jets = Collection(event, "Jet")

		nMuons = len(muons)
		nElectrons = len(electrons)
		nLeptons = nMuons + nElectrons

		self.out.fillBranch("nMuons",nMuons)
		self.out.fillBranch("nElectrons",nElectrons)
		self.out.fillBranch("nLeptons",nLeptons)

		return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
susy1LeptonBase = lambda : susy1LeptonBaseProducer()

