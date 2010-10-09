import FWCore.ParameterSet.Config as cms
process = cms.Process("TrackerMu")

#process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
#process.TrackAssociatorByHits.Cut_RecoToSim = cms.double(0.5)
#process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
#process.TrackAssociatorByChi2ESProducer.chi2cut = cms.double(1e9)

#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
#process.options = cms.untracked.PSet( SkipEvent = cms.untracked.vstring('ProductNotFound') )
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('START38_V9::All')

process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("SimGeneral.TrackingAnalysis.trackingParticlesNoSimHits_cfi")
process.load("SimMuon.MCTruth.MuonAssociatorByHitsESProducer_NoSimHits_cfi")

process.source = cms.Source('PoolSource',
fileNames = cms.untracked.vstring("file:/home/zhangjin/ppMuX_START37_V5_S09-v1_RECO.root")
#fileNames = cms.untracked.vstring("rfio:/castor/cern.ch/user/z/zhangjin/testCSCEff.root")
#,eventsToProcess = cms.untracked.VEventRange("1:5820068-1:5820068")
#fileNames = cms.untracked.vstring("file:/afs/cern.ch/user/z/zhangjin/scratch0/ProblemEvt.root")
)

process.MuonSelector = cms.EDFilter('CSCPriEff',
                     FileName = cms.untracked.string("TM.root"),
                     StandardMuonCuts = cms.untracked.vstring("GlobalMuonPromptTight","TMLastStationLoose","TMLastStationTight","TMLastStationAngLoose","TMLastStationAngTight"),
                     maxChamberDist = cms.untracked.double(-3.),
                     maxChamberDistPull = cms.untracked.double(-3.),
                     tracksTag = cms.untracked.InputTag("globalMuons")
)

process.MuonFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("muons"),
    minNumber = cms.uint32(1))

process.MuonMC = cms.Path(process.MuonFilter+process.mix*process.trackingParticlesNoSimHits*process.MuonSelector)


from FWCore.MessageLogger.MessageLogger_cfi import *
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000   #quench the message logger (optional)
