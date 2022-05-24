#Loading necessary libraries
import FWCore.ParameterSet.Config as cms
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
import PhysicsTools.PythonAnalysis.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

import os 
import sys

relBase = os.environ['CMSSW_BASE']

#---- sys.argv takes the parameters given as input cmsRun PhysObjectExtractor/python/poet_cfg.py <isData (default=False)> <doPat (default=False)> 
#----  e.g: cmsRun PhysObjectExtractor/python/poet_cfg.py True True
#---- NB the first two parameters are always "cmsRun" and the config file name
#---- Work with data (if False, assumed MC simulations)
#---- This needs to be in agreement with the input files/datasets below.
if len(sys.argv) > 2:
    isData = eval(sys.argv[2])
else:
    isData = False

#---- Flag for using the Physics Analysis Toolkit for jets and MET
if len(sys.argv) > 3:
    doPat = eval(sys.argv[3])
else:
    doPat = False

process = cms.Process('HiForest')
process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))

#Number of events: put '-1' unless testing
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

#HiForest script init
process.load("HiForest_cff")
process.HiForest.inputLines = cms.vstring("HiForest V3",)
version = 'no git info'
process.HiForest.HiForestVersion = cms.string(version)

goodJSON = 'Cert_181530-183126_HI7TeV_PromptReco_Collisions11_JSON.txt'
myLumis = LumiList.LumiList(filename = goodJSON).getCMSSWString().split(',')
import FWCore.Utilities.FileUtils as FileUtils
files2011data = FileUtils.loadListFromFile ('CMS_HIRun2011_HIHighPt_RECO_15Apr2013-v1_root_file_index.txt')
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(*files2011data    
    )
)
process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
process.source.lumisToProcess.extend(myLumis)

#Global Tag: change the name according to the instructions
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/GR_R_44_V15.db')
process.GlobalTag.globaltag = 'GR_R_44_V15::All'
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.HiForest.GlobalTagLabel = process.GlobalTag.globaltag

#Define the output root file (change each run not to overwrite previous output)
process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string("HiForestAOD_DATAtest2011.root"))

#Init Trigger Analyzer
process.hltanalysis = cms.EDAnalyzer('TriggerInfoAnalyzer',
                              processName = cms.string("HLT"),
                              triggerName = cms.string("@"),         
                              datasetName = cms.string("HIDiMuon"),  #'HICorePhysics' to look at Core Physics only
                              triggerResults = cms.InputTag("TriggerResults","","HLT"),
                              triggerEvent   = cms.InputTag("hltTriggerSummaryAOD","","HLT")                             
                              )

#Collect event data
process.demo = cms.EDAnalyzer('Analyzer') #present analyzer is for muons - see details in Analyzer.cc for possible modifications

process.myevents = cms.EDAnalyzer('EventAnalyzer')


process.mytracks= cms.EDAnalyzer('TrackAnalyzer')

#---- Jet correction paths -- these correspond to the Global Tag. **Run jec_cfg.py first to get .txt files!!**
JecString = 'START53_V27_'
if isData: JecString = 'FT53_V21A_AN6_'

#---- Jets are simpler to work with in "Physics Analysis Toolkit" format. See more at https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookPAT
if doPat:
	#---- Load PAT configs and build some light sequences to process jets and MET
	process.load('PhysicsTools.PatAlgos.producersLayer1.jetProducer_cff')
	process.load('PhysicsTools.PatAlgos.producersLayer1.metProducer_cff')
	process.load('PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi')
	process.patCandidates = cms.Sequence(process.makePatJets+process.makePatMETs)
	process.selectedPatCandidates = cms.Sequence(process.selectedPatJets)
	process.patDefaultSequence = cms.Sequence(process.patCandidates * process.selectedPatCandidates)
	process.load('RecoJets.Configuration.RecoPFJets_cff')
	from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection, runBTagging
	from PhysicsTools.PatAlgos.tools.coreTools import runOnData

	#---- Choose which jet correction levels to apply
	jetcorrlabels = ['L1FastJet','L2Relative','L3Absolute']
	if isData:
		#---- For data we need to remove generator-level matching processes
		runOnData(process, ['Jets','METs'], "", None, [])
		jetcorrlabels.append('L2L3Residual')

	#---- Configure the addJetCollection tool
	#---- This process will make corrected jets with b-tagging included, and will make Type1-corrected MET
	process.ak5PFJets.doAreaFastjet = True
	addJetCollection(process,cms.InputTag('ak5PFJets'),
			 'AK5', 'PFCorr',
			 doJTA        = True,
			 #doBTagging   = True, 
			 jetCorrLabel = ('AK5PF', cms.vstring(jetcorrlabels)),
			 doType1MET   = True,
			 doL1Cleaning = False,
			 doL1Counters = False,
			 doJetID      = True,
			 jetIdLabel   = "ak5",
			 ) 
 
	#---- Configure the POET jet analyzer
	#---- Don't forget to run jec_cfg.py to get these .txt files!
	process.myjets= cms.EDAnalyzer('PatJetAnalyzer',
				       InputCollection = cms.InputTag("selectedPatJetsAK5PFCorr"),
				       isData = cms.bool(isData),
				       jecUncName = cms.FileInPath('HiForest/HiForestProducer/JEC/'+JecString+'Uncertainty_AK5PF.txt'), 
				       jerResName = cms.FileInPath('HiForest/HiForestProducer/JEC/JetResolutionInputAK5PF.txt')         
				       )
else:
	if not isData:
		#---- Get non-PAT access to the jet flavour information
		#from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
		from PhysicsTools.JetMCAlgos.SelectPartons_cff import myPartons
                #process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone()
		process.selectedHadronsAndPartons = myPartons.clone()
                #from PhysicsTools.JetMCAlgos.AK5PFJetsMCFlavourInfos_cfi import ak5JetFlavourInfos
		from PhysicsTools.JetMCAlgos.AK5CaloJetsMCFlavour_cff import AK5byValAlgo
                #process.jetFlavourInfosAK5PFJets = AK5byValAlgo.clone()

	#---- Configure the POET jet analyzer
	#---- Don't forget to run jec_cfg.py to get these .txt files!
	process.myjets= cms.EDAnalyzer('JetAnalyzer',
				       InputCollection = cms.InputTag("ak5PFJets"),
				       isData = cms.bool(isData),
				       jecL1Name = cms.FileInPath('HiForest/HiForestProducer/JEC/'+JecString+'L1FastJet_AK5PF.txt'), 
				       jecL2Name = cms.FileInPath('HiForest/HiForestProducer/JEC/'+JecString+'L2Relative_AK5PF.txt'),
				       jecL3Name = cms.FileInPath('HiForest/HiForestProducer/JEC/'+JecString+'L3Absolute_AK5PF.txt'),
				       jecResName = cms.FileInPath('HiForest/HiForestProducer/JEC/'+JecString+'L2L3Residual_AK5PF.txt'),
				       jecUncName = cms.FileInPath('HiForest/HiForestProducer/JEC/'+JecString+'Uncertainty_AK5PF.txt'),
				       jerResName = cms.FileInPath('HiForest/HiForestProducer/JEC/JetResolutionInputAK5PF.txt')
				       )

process.dump=cms.EDAnalyzer('EventContentAnalyzer') #easy check of Event structure and names without using the TBrowser

process.ana_step = cms.Path(process.hltanalysis+
		  	    #process.dump+  #uncomment if necessary to check the name. Do not forget to change the number of events to '1'
			    process.demo+
			    process.mytracks+
			    process.myjets+
                            process.HiForest 
)

