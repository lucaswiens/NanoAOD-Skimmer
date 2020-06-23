# Selection Modules
from Susy1LeptonAnalysis.NanoAODSkimmer.modules.susy1LeptonProducer import susy1LeptonProducer
from Susy1LeptonAnalysis.NanoAODSkimmer.modules.susy1LeptonBaseSelector import susy1LeptonSelector
## from PhysicsTools.NanoAODTools.postprocessing.modules.selection.selectionProducer import selection
##
## #from PhysicsTools.NanoAODTools.postprocessing.modules.susy.genParticleProducer import genAll
## #jme
## #from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import createJMECorrector
## #common
## from PhysicsTools.NanoAODTools.postprocessing.modules.common.PrefireCorr import *
## from PhysicsTools.NanoAODTools.postprocessing.modules.common.lepSFProducer import lepSF
##
## #btv
## from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import btagSF2016, btagSF2017
##
## preFireCorr2016 = lambda : PrefCorr(jetroot="L1prefiring_jetpt_2016BtoH.root", jetmapname="L1prefiring_jetpt_2016BtoH", photonroot="L1prefiring_photonpt_2016BtoH.root", photonmapname="L1prefiring_photonpt_2016BtoH", branchnames=["PrefireWeight","PrefireWeight_Up", "PrefireWeight_Down"])
## preFireCorr2017 = lambda : PrefCorr(jetroot="L1prefiring_jetpt_2017BtoF.root", jetmapname="L1prefiring_jetpt_2017BtoF", photonroot="L1prefiring_photonpt_2017BtoF.root", photonmapname="L1prefiring_photonpt_2017BtoF", branchnames=["PrefireWeight","PrefireWeight_Up", "PrefireWeight_Down"])
