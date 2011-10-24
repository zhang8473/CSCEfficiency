HLTProc = "HLT"
inputfiles="rfio:/castor/cern.ch/cms/store/data/Run2011B/SingleMu/RAW/v1/000/175/888/02370460-2DDA-E011-9F88-003048F11C58.root"
outputfileName="MuTrkCand.root"

import sys
if sys.argv[0] == "cmsRun":
    for option in sys.argv:
        if "file:" in option or "rfio:" in option:
            inputfiles = [option]
        if "out:" in option:
            outputfilename = option.replace("out:","")
        if "hltproc:" in option:
            HLTProc = option.replace("hltproc:","")

import FWCore.ParameterSet.Config as cms
process = cms.Process("MuonSysPriEff")

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#from Configuration.PyReleaseValidation.autoCond import autoCond
#from Configuration.AlCa.autoCond import autoCond#after 4_3_0
#process.GlobalTag.globaltag = autoCond['startup']
#process.GlobalTag.globaltag = "GR_R_44_V1::All"
process.GlobalTag.globaltag = "GR_P_V25::All"

process.load("Configuration.StandardSequences.Services_cff")

process.source = cms.Source('PoolSource',
fileNames = cms.untracked.vstring(inputfiles)
#,eventsToProcess = cms.untracked.VEventRange('1:24575813-1:24575813')
#,skipEvents=cms.untracked.uint32(28)
)
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.MuTrkCand = cms.EDAnalyzer('MuonSysPriEff',
                                  FileName = cms.string(outputfileName),
                                  MuonPtCut = cms.double(3.),
                                  MuonEtaMin = cms.double(0.9),
                                  MuonEtaMax = cms.double(2.4),
                                  ThrowBadEvents = cms.bool(True),
                                  MinNumberOfMuons = cms.untracked.uint32(1),
                                  CSCDigisTag = cms.InputTag("csctfDigis"),
                                  TriggerResultsTag = cms.InputTag('TriggerResults','',HLTProc),
                                  triggerEventTag = cms.untracked.InputTag('hltTriggerSummaryAOD','',HLTProc),
                                  HLTObj_HLTNames = cms.vstring("HLT_Mu5_v10"),
                                  hltFilterNames = cms.VInputTag(cms.InputTag('hltSingleMu24L3Filtered24','',HLTProc)),
                                  StandardMuonCuts = cms.untracked.vstring("GlobalMuonPromptTight","TMLastStationLoose","TMLastStationTight","TMLastStationAngLoose","TMLastStationAngTight"),
                                  minTrackHits = cms.untracked.uint32(3)
)

########### Event cleaning and  Trigger selection ###########
# Select events based on the HLTtriggers....singleJet and BTag triggers
# Use the instructions provided at:
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/TriggerResultsFilter
# This eases the trigger selection for different HLT menus and also takes care of wildcard and trigger versioning
#############################################################
process.MuTrkCand.VertexFilterPSet = cms.PSet ( vertexCollection = cms.InputTag("offlinePrimaryVertices"),
                                                   minimumNDOF = cms.uint32(4),
                                                   maxAbsZ = cms.double(24),
                                                   maxd0 = cms.double(2.0)
                                                 )

process.MuTrkCand.noscrapingPSet = cms.PSet( numtrack = cms.untracked.uint32(10),
                                thresh = cms.untracked.double(0.2)
                                )

process.MuTrkCand.HLTFilterPSet = cms.PSet ( triggerConditions = cms.vstring(
                                                #SingleMuStream Triggers 
                                                "HLT_IsoMu*",
                                                "HLT_L1SingleMu*",
                                                "HLT_L1DoubleMu*",
                                                "HLT_L2Mu*",
                                                "HLT_Mu??_v*",
                                                "HLT_Mu?_v*",
                                                #DoubleMuStream Triggers
                                                "HLT_DoubleMu*",
                                                "HLT_L1DoubleMu*",
                                                "HLT_L2DoubleMu*",
                                                "HLT_Mu*_Jet*_v*"
                                                ),
                                              hltResults = cms.InputTag('TriggerResults','',HLTProc),
                                              l1tResults = cms.InputTag( "" ),
                                              l1tIgnoreMask = cms.bool(False),
                                              l1techIgnorePrescales   = cms.bool(False),
                                              daqPartitions           = cms.uint32(0x01),
                                              throw = cms.bool( False ) #set to false to deal with missing triggers while running over different trigger menus
                                            )
process.MuonMC = cms.Path(process.RawToDigi*process.reconstruction*process.MuTrkCand)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",moduleMemorySummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 100   #quench the message logger (optional)
process.MessageLogger.cerr_stats.threshold = cms.untracked.string('WARNING')
