
# Configuration file for CSCTriggerPrimitive study.
# Chaouki Boulahouache June-2010.
# Jinzhong Zhang June-2012.
#
# Import the right stuff
#
datatype = "RAW"
#Candidates are data: "RAW" "RAW-RECO" "FEVT"
#mc: in order of suggestions: "GEN-RAWDEBUG"(mc) "GEN-SIM-RAW"(mc) "GEN-RAW"(mc) "GEN-SIM".....
# if no RAW exists, you need to choose a mixing module manually
# pileup tracks cannot be identified in "GEN-RAW" sample  

HLTProcessName='HLT'

import sys,os

if "GEN" in datatype or "SIM" in datatype:
  runOnMC = True
else:
  runOnMC = False

CMSSW_Version_=os.getenv("CMSSW_VERSION")
Version_Number_=CMSSW_Version_.split('_')

if ( int(Version_Number_[1])<5 )  or ( int(Version_Number_[1])==5 and int(Version_Number_[2])<2 ):
  After520=False
else:
  After520=True

if ( int(Version_Number_[1])<4 ) or ( int(Version_Number_[1])==4 and int(Version_Number_[2])<3 ):
  After430=False
else:
  After430=True

if After430:
  from Configuration.AlCa.autoCond import autoCond # for >=CMSSW_4_3_0
else:
  from Configuration.PyReleaseValidation.autoCond import autoCond

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
#
##################################################################################################
process = cms.Process("ANALYSISNTUPLE")
###################################################################################################

# import of standard configurations
process.load("FWCore/MessageService/MessageLogger_cfi")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
# CSC TF
process.load("EventFilter.CSCTFRawToDigi.csctfunpacker_cfi")
process.csctfunpacker.producer = cms.InputTag("rawDataCollector")

process.schedule=cms.Schedule()

if runOnMC:
  process.load("Configuration.StandardSequences.RawToDigi_cff")
  process.load("Configuration.StandardSequences.Reconstruction_cff")
  process.GlobalTag.globaltag = autoCond['startup']
  print "The global tag is ", process.GlobalTag.globaltag
##------------------  IMPORT MODULES FOR PRECISE MC MATCHING  ------------
  if "DEBUG" not in datatype:
    process.load('Configuration/StandardSequences/Services_cff')
    process.load("SimGeneral.MixingModule.mixNoPU_cfi")
    process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")
    process.load("Configuration.StandardSequences.Digi_cff")
#process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
#process.TrackAssociatorByChi2ESProducer.chi2cut = cms.double(1e9)
    if "SIM" not in datatype and "FEVT" not in datatype:
      process.load('SimG4Core/Application/g4SimHits_cfi')
      process.trackingTruth = cms.Path(process.g4SimHits*cms.SequencePlaceholder("mix")*process.trDigi*process.trackingParticles)
    else:
      process.trackingTruth = cms.Path(cms.SequencePlaceholder("mix")*process.trDigi*process.trackingParticles)

  if "RAW" not in datatype and "FEVT" not in datatype:
    if After520:
      process.load("SimGeneral.MixingModule.mix_2012_Startup_50ns_PoissonOOTPU_cfi")
    else:
      process.load("SimGeneral.MixingModule.mix_E7TeV_Summer_2011_50ns_PoissonOOT")
    process.load('Configuration.StandardSequences.Digi_cff')
    process.load('Configuration.StandardSequences.SimL1Emulator_cff')
    process.load('Configuration.StandardSequences.DigiToRaw_cff')
    process.load('HLTrigger.Configuration.HLT_GRun_cff')
    HLTProcessName=process._Process__name
    process.CSCHaloData.HLTResultLabel = cms.InputTag("TriggerResults","",HLTProcessName)
    process.simtoraw=cms.Path(process.pdigi*process.SimL1Emulator*process.DigiToRaw)
    process.schedule.extend( [process.simtoraw] )
    process.schedule.extend( process.HLTSchedule )
  else:
    if "DEBUG" not in datatype:# Playback
      del process.RandomNumberGeneratorService.generator
      process.RandomNumberGeneratorService.restoreStateLabel = cms.untracked.string('randomEngineStateProducer')
      process.mix.playback = cms.untracked.bool(True)

  process.schedule.extend( [process.trackingTruth] )

  process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
  process.TrackAssociatorByHits.Cut_RecoToSim = cms.double(0.5)
else:
  process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
  if After520:
    process.load("Configuration.StandardSequences.Reconstruction_Data_cff")
  else:
    process.load("Configuration.StandardSequences.Reconstruction_cff")
  process.GlobalTag.globaltag = autoCond['com10']
  #Load bad chambers according to the run numbers
  process.load("CalibMuon.Configuration.getCSCConditions_frontier_cff")
  process.cscConditions.connect='frontier://FrontierProd/CMS_COND_31X_CSC'
  process.cscConditions.toGet = cms.VPSet(cms.PSet(
      record = cms.string('CSCBadChambersRcd'),
      tag = cms.string('CSCBadChambers_RunDependent_offline')
      ))
  process.cscConditions.loadAll = cms.bool(True)
  process.cscConditions.timetype = cms.string('runnumber')
  process.cscConditions.DBParameters = cms.PSet(
    authenticationPath = cms.untracked.string('/nfshome0/popcondev/conddb/'),
    authenticationMethod = cms.untracked.uint32(1)
    )
  process.es_prefer_cscConditions = cms.ESPrefer("PoolDBESSource","cscConditions")
  #end of loading bad chambers according to the run numbers

if "RECO" in datatype or "FEVT" in datatype:
  process.raw2digiandreconstruction = cms.Path(process.RawToDigi)
else:
  process.raw2digiandreconstruction = cms.Path(process.RawToDigi*process.reconstruction)

process.schedule.extend( [process.raw2digiandreconstruction] )

# ecal mapping
process.load("Geometry.EcalMapping.EcalMapping_cfi")

##------------------  IMPORT MODULES NEEDED FOR JPT  ------------
if After520:
   process.load("Configuration.Geometry.GeometryIdeal_cff")
else:
  process.load("Configuration.StandardSequences.Geometry_cff")

from Configuration.StandardSequences.MagneticField_cff import *

process.source = cms.Source("PoolSource",
#                            fileNames = cms.untracked.vstring("dcap:://pnfs/cms/WAX/11/store/mc/Fall11/DYJetsToEE-MuMu_PtZ-100To300_TuneZ2_7TeV-calchep-pythia/GEN-RAW/PU_S6_START42_V14B-v1/0000/2E0052CC-62B8-E111-ACBE-00261894388F.root")
#                            fileNames = cms.untracked.vstring("file:/uscms/home/zhangjin/nobackup/Test_DYToMuMu_M-20_CT10_GEN-SIM-RAW-HLTDEBUG-RECO.root")
#                            fileNames = cms.untracked.vstring("file:/uscms/home/zhangjin/nobackup/Test_DYMuMu200_42_RAW.root")
#                            fileNames = cms.untracked.vstring("dcap:://pnfs/cms/WAX/11/store/mc/Summer12/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/GEN-SIM-RAW/PU_S7_START52_V9-v5/0000/A6A3E090-B498-E111-A4AB-003048C68AA6.root")
                            #fileNames = cms.untracked.vstring("dcap:://pnfs/cms/WAX/11/store/mc/Summer12_DR53X/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM-RECODEBUG/DEBUG_PU_S10_START53_V7A-v1/0000/A8B78296-D6E0-E111-AD1F-003048678B1A.root")
                            #fileNames = cms.untracked.vstring("dcap:://pnfs/cms/WAX/11/store/mc/Summer12/DYToMuMu_M_20_TuneZ2star_8TeV_pythia6/GEN-SIM/START50_V13-v1/0000/001E08EA-EC56-E111-A0C6-002618943915.root")#/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12-START50_V13-v1/GEN-SIM
                            #fileNames = cms.untracked.vstring("dcap:://pnfs/cms/WAX/11/store/mc/Summer12/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/GEN-SIM/START52_V9-v3/0008/FEF6E3F9-FCD3-E111-8322-0025901D0C52.root")#/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/Summer12-START52_V9-v3/GEN-SIM
                            #fileNames = cms.untracked.vstring("dcap:://pnfs/cms/WAX/11/store/mc/Summer12/JPsiToMuMu_2MuPtEtaFilter_tuneD6T_8TeV-pythia6-evtgen/GEN-SIM/START52_V9-v3/0001/FEDFAB36-3BF8-E111-9DF6-68B599B94F60.root")#/JPsiToMuMu_2MuPtEtaFilter_tuneD6T_8TeV-pythia6-evtgen/Summer12-START52_V9-v3/GEN-SIM
                            fileNames = cms.untracked.vstring("dcap://cmsdca.fnal.gov:24137/pnfs/fnal.gov/usr/cms/WAX/11/store/data/Run2012C/SingleMu/RAW/v1/000/202/016/98D875A0-10F3-E111-B2F2-001D09F2AD84.root")#/SingleMu/Run2012C-v1/RAW
                            ,skipEvents=cms.untracked.uint32(0)
)
#
#  Number of events to be processed...
#
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(
  SkipEvent = cms.untracked.vstring( "Error: uninitialized ProxyBase used" ),
  IgnoreCompletely = cms.untracked.vstring( "ProductNotFound" )
  )

from RecoMuon.TrackingTools.MuonSegmentMatcher_cff import *
process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff") 
from TrackingTools.TrackAssociator.default_cfi import *
from RecoMuon.MuonIsolationProducers.trackExtractorBlocks_cff import MIsoTrackExtractorBlock

process.aodDump = cms.EDAnalyzer('TPTrackMuonSys',
                                 TrackAssociatorParameterBlock,
                                 MuonSegmentMatcher,
                                 TrackExtractor=cms.PSet(MIsoTrackExtractorBlock),
                                 rootFileName   = cms.untracked.string('CSCPFG_Ineff_DATA.root'),
                                 CSCUseTimingCorrections = cms.bool( True ),
                                 CSCUseGasGainCorrections = cms.bool( True ),
                                 isMC            = cms.untracked.bool(runOnMC),
                                 mcTag           = cms.untracked.InputTag('genParticles'),
                                 vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                 gTracksTag      = cms.untracked.InputTag('generalTracks'),
                                 trackProducer   = cms.InputTag('csctfunpacker:'),
                                 readBadChannels = cms.bool(True),
                                 readBadChambers = cms.bool(True),
                                 hltTag      = cms.untracked.InputTag("TriggerResults","",HLTProcessName),
                                 hltEvTag    = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
                                 HLTMuTrgNames =  cms.vstring("HLT_Mu?_v*","HLT_Mu??_v*","HLT_Mu???_v*","HLT_IsoMu?_v*","HLT_IsoMu??_v*","HLT_IsoMu???_v*","HLT_L2Mu?_v*","HLT_L2Mu??_v*","HLT_L2Mu???_v*","HLT_SingleMu*","HLT_L1SingleMu*"),
                                 HLTDiMuTrgName =  cms.string("HLT_DoubleMu?_v*"),
                                 #                      hltEvTag    = cms.untracked.InputTag("hltTriggerSummaryAOD","","REDIGI36X"),
                                 L1extraTag   = cms.untracked.InputTag("l1extraParticles"), 
                                 dedxTag         =  cms.untracked.InputTag('dedxHarmonic2')
                                 )

########################################
process.primaryVertexFilter = cms.EDFilter("VertexSelector",
   src = cms.InputTag("offlinePrimaryVertices"),
   cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"), # tracksSize() > 3 for the older cut
   filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
)

process.noscraping = cms.EDFilter("FilterOutScraping",
   applyfilter = cms.untracked.bool(True),
   debugOn = cms.untracked.bool(False),
   numtrack = cms.untracked.uint32(10),
   thresh = cms.untracked.double(0.25)
)

process.CSCTFunpacker=cms.Path(process.csctfunpacker)
process.filtersteps=cms.Path(process.primaryVertexFilter*process.noscraping)
process.outputstep = cms.EndPath(process.aodDump)

process.schedule.extend( [process.CSCTFunpacker,process.filtersteps,process.outputstep] )

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1   #quench the message logger (optional)
process.MessageLogger.cerr_stats.threshold = cms.untracked.string('WARNING')
"""
#For DEBUG:
process.load('Configuration.EventContent.EventContent_cff')
process.out = cms.OutputModule("PoolOutputModule",
                               splitLevel = cms.untracked.int32(0),
                               outputCommands = cms.untracked.vstring(
                                           'drop *_*_*_*',
                                           'keep *_mix_*_*',
                                           'keep *_generator_*_*',
                                           'keep *_randomEngineStateProducer_*_*',
                                           'keep *_mergedtruth__*',
                                           'keep *_mergedtruth_*_*',
                                           'keep *_muons_*_*'
                                           ),
                               fileName = cms.untracked.string('MYDEBUG.root'),
                               )
#process.outputstep=cms.EndPath(process.out)
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",moduleMemorySummary = cms.untracked.bool(True))
  if After520:
    process.highlevelreco=cms.Sequence(process.egammaHighLevelRecoPrePF*
                                       process.particleFlowReco*
                                       process.egammaHighLevelRecoPostPF*
                                       process.regionalCosmicTracksSeq*
                                       process.muoncosmichighlevelreco*
                                       process.muonshighlevelreco*
                                       process.particleFlowLinks)
  else:
    process.highlevelreco=cms.Sequence(process.egammaHighLevelRecoPrePF*
                                       process.particleFlowReco*
                                       process.egammaHighLevelRecoPostPF)
  process.reconstruction = cms.Sequence(process.localreco*process.globalreco*process.highlevelreco)
"""
