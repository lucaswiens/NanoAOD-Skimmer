#!/usr/bin/env python2
"""
author:	Lucas Wiens
mail:	lucas.wiens@desy.de
"""

import sys, os
import argparse
import subprocess
import shutil
from re import findall

def createFileAndModuleList(sampleFileName):
	fileList, moduleList = [], []
	sampleFile = open(sampleFileName, "r")
	for sample in sampleFile:
		sample = sample.strip()
		if not (sample.startswith("#") or sample in ["", "\n", "\r\n"]):
			fileList.append(createFileList(sample))
			moduleList.append(prepareModuleListAndArguments(sample))
	sampleFile.close()
	return fileList, moduleList

def createFileList(sample):
	return subprocess.check_output("dasgoclient -query=\"file dataset=" + str(sample) + "\"", shell=True).split()

def prepareModuleListAndArguments(sample):
	isMC = False
	isSig = False
	isUSER = False
	year = "unknown"
	runPeriod = "unkown"
	jesUncert = "unknown"
	redoJEC = "unknown"
	isFastSim = False
	module = "unknown"
	if "RunIISummer16" in str(sample) or "Run2016" in str(sample):
		year = 2016
		if "PUSummer16v3Fast" in str(sample):
			isFastSim == True
		else:
			isFastSim == False
	elif "RunIIFall17" in str(sample) or "Run2017" in str(sample):
		year = 2017
		if "TuneCP2" in str(sample):
			isFastSim == True
		else:
			isFastSim == False
	if "/NANOAODSIM" in sample:
		isMC = True
	if "/SMS-T1tttt" in sample:
		isSig = True
	if isMC and not isSig:
		#module = "susy_1l_FiltersMC,jecUncert,susy1lepTOPMC,susy_1l_gen"#,xsec,genpartsusymod
		module = ""#,xsec,genpartsusymod
		if year == 2016:
			#module +="preFireCorr2016"
			module +="susy1LeptonProducer"
		if year == 2017:
			# Temporarly use the jetmet uncertainty for 2016
			module +="preFireCorr2017"
		# if "TTJets" in str(sample) and year == 2016: module +=",susy_1l_nISR16,susy1lepTT_syst"
		# if "TTJets" in str(sample) and year == 2017: module +=",susy_1l_nISR17,susy1lepTT_syst"
		# if "WJets"  in str(sample): module +=",susy1lepWJets_syst"
	elif isMC and isSig:
		module = ""#,xsec,genpartsusymod
			# Temporarly use the jetmet uncertainty for 2016
		if year == 2016: module +="susy1LeptonProducer"
		if year == 2017: module +="susy1LeptonProducer"
		#if year == 2016: module +="preFireCorr2016,lepSF,btagSF2016"
		#if year == 2017: module +="preFireCorr2017,lepSF,btagSF2017"
	else:
		if year == 2016:
			module = ""
		elif year == 2017:
			module = ""

	if not isMC:
		runPeriod = findall("Run201..", sample)[0][-1]
	return module, isMC, isSig, year, runPeriod, isFastSim#, jesUncert, redoJEC

def getOSVariable(Var):
	try:
		variable = os.environ[Var]
	except KeyError:
		print "Please set the environment variable " + Var
		sys.exit(1)
	return variable



if __name__=="__main__":
	date = subprocess.check_output("date +\"%Y_%m_%d\"", shell=True).replace("\n", "")
	parser = argparse.ArgumentParser(description="Runs a NAF batch system for nanoAOD", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-i", "--input-file", required=True, help="Path to the file containing a list of samples.")
	parser.add_argument("-o", "--output", help="Path to the output directory", default = "batch/" + date)

	args = parser.parse_args()

	cmsswBase = getOSVariable("CMSSW_BASE")
	workarea = "%s/src/Susy1LeptonAnalysis/NanoAODSkimmer/batch" % cmsswBase
	#X509 = getOSVariable("X509_USER_PROXY")
	X509 = "/tmp/x509up_u33974"

	condTEMP = "./templates/submit.condor"
	wrapTEMP = "./templates/wrapnanoPost.sh"

	if  os.path.exists(args.output):
		keepDirectory = raw_input("Output directory already exists: " + str(args.output) + " Do you want to remove it [y/n]: ")
		if ( "y" in keepDirectory or "Y" in keepDirectory or "Yes" in keepDirectory):
			shutil.rmtree(str(args.output))
			os.makedirs(str(args.output))
			os.makedirs(str(args.output) + "/samples")
			os.makedirs(str(args.output) + "/condor")
			os.makedirs(str(args.output) + "/wrapper")
			os.makedirs(str(args.output) + "/logs")
			os.makedirs(str(args.output) + "/tree")
			os.makedirs(str(args.output) + "/output")
		elif ( "N" in keepDirectory or  "n" in keepDirectory or  "No" in keepDirectory ): print str(args.output) , "will be ovewritten by the job output -- take care"
		else:
			raise ValueError( "invalid input, answer with \"Yes\" or \"No\"")
	else:
		os.makedirs(str(args.output))
		os.makedirs(str(args.output) + "/samples")
		os.makedirs(str(args.output) + "/condor")
		os.makedirs(str(args.output) + "/wrapper")
		os.makedirs(str(args.output) + "/logs")
		os.makedirs(str(args.output) + "/tree")
		os.makedirs(str(args.output) + "/output")

	sampleFile = open(args.input_file, "r")
	for sample in sampleFile:
		sample = sample.strip()
		if not (sample.startswith("#") or sample in ["", "\n", "\r\n"]):
			fileList = createFileList(sample)
			moduleList, isMC, isSig, year, runPeriod, isFastSim = prepareModuleListAndArguments(sample)
			sampleName = sample.replace("/", "_")[1:]

			file = open(args.output + "/samples/" + sampleName + ".txt", "w+")
			for filename in fileList:
				file.write("root://cms-xrd-global.cern.ch/" + str(filename) + "\n")
			file.close()

			i = 1
			logDirectory = args.output + "/logs"
			for filename in fileList:
				os.system("cp " + condTEMP + " " + args.output + "/condor/" + sampleName + str(i) + ".submit")
				submitFileContent = open(args.output + "/condor/" + sampleName + str(i) + ".submit").read()
				submitFileContent = submitFileContent.replace("@EXECUTABLE", args.output + "/wrapper/" + sampleName + str(i))
				submitFileContent = submitFileContent.replace("@LOGS", logDirectory)
				submitFileContent = submitFileContent.replace("@X509", X509)
				submitFileContent = submitFileContent.replace("@TIME", "60*60*48")

				submitFile = open(args.output + "/condor/" + sampleName + str(i) + ".submit", "w")
				submitFile.write(submitFileContent)
				submitFile.close()

				skimTreeOutput = "_".join(filename.split("/")[3:5])
				os.system("cp " + wrapTEMP + " " + args.output + "/wrapper/" + sampleName + str(i))
				wrapperFileContent = open(args.output + "/wrapper/" + sampleName + str(i)).read()
				wrapperFileContent = wrapperFileContent.replace("@WORKDIR", workarea)
				wrapperFileContent = wrapperFileContent.replace("@EXEDIR", workarea)
				wrapperFileContent = wrapperFileContent.replace("@MODULES", moduleList)
				wrapperFileContent = wrapperFileContent.replace("@ISMC", "--is-mc"*isMC)
				wrapperFileContent = wrapperFileContent.replace("@ISFASTSIM", "--is-fastsim"*isFastSim)
				wrapperFileContent = wrapperFileContent.replace("@YEAR", str(year))
				wrapperFileContent = wrapperFileContent.replace("@RUNPERIOD", runPeriod)
				wrapperFileContent = wrapperFileContent.replace("@OUTPUT", "SkimmingOutput/" + date + "/" + skimTreeOutput)
				wrapperFileContent = wrapperFileContent.replace("@SKIMTREELOCATION", args.output)
				wrapperFileContent = wrapperFileContent.replace("@INPUTFILE", "root://cms-xrd-global.cern.ch/" + filename)
				wrapperFileContent = wrapperFileContent.replace("@X509", X509)

				wrapperFile = open(args.output + "/wrapper/" + sampleName + str(i), "w")
				wrapperFile.write(wrapperFileContent)
				wrapperFile.close()
				file = open(args.output + "/submitAllViaHTC", "a")
				file.write("condor_submit -name s02 " + args.output + "/condor/" + sampleName + str(i) + ".submit\n")
				file.close()
				i +=1
	os.system("chmod +744 " + args.output + "/submitAllViaHTC")
	print "submitAllViaHTC created in " + args.output + " to submit all jobs"
	sampleFile.close()
