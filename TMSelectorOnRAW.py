import FWCore.ParameterSet.Config as cms
process = cms.Process("TrackerMu")

#process.load("SimTracker.TrackHistory.PlaybackWithReco_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.StandardSequences.Digi_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START38_V12::All')
process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")

process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.TrackAssociatorByHits.Cut_RecoToSim = cms.double(0.5)
process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
process.TrackAssociatorByChi2ESProducer.chi2cut = cms.double(1e9)

#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
#process.options = cms.untracked.PSet( SkipEvent = cms.untracked.vstring('ProductNotFound') )

# Playback
process.load("Configuration.StandardSequences.Services_cff")
del process.RandomNumberGeneratorService.generator
process.RandomNumberGeneratorService.restoreStateLabel = cms.untracked.string('randomEngineStateProducer')
process.mix.playback = cms.untracked.bool(True)

process.source = cms.Source('PoolSource',
fileNames = cms.untracked.vstring("rfio:/castor/cern.ch/user/z/zhangjin/ppMuX_START37_V5_S09-v1_RAW.root")
#,eventsToProcess = cms.untracked.VEventRange("1:5820068-1:5820068")
)

process.TMSelector = cms.EDFilter('CSCPriEff',
                     FileName = cms.untracked.string("TM.root"),
                     StandardMuonCuts = cms.untracked.vstring("GlobalMuonPromptTight","TMLastStationLoose","TMLastStationTight","TMLastStationAngLoose","TMLastStationAngTight"),
                     maxChamberDist = cms.untracked.double(-3.),
                     maxChamberDistPull = cms.untracked.double(-3.),
                     minTrackHits = cms.untracked.uint32(3)
)

process.trackingTruth = cms.Sequence(process.mix*process.doAllDigi*process.mergedtruth)

process.MuonMC = cms.Path(process.RawToDigi*process.trackingTruth*process.reconstruction*process.TMSelector)


from FWCore.MessageLogger.MessageLogger_cfi import *
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10   #quench the message logger (optional)
