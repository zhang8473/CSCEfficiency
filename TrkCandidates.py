import FWCore.ParameterSet.Config as cms
process = cms.Process("TrackerMu")

#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
#process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
#process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cff")

"""
process.TrackAssociatorByChi2ESProducer = cms.ESProducer("TrackAssociatorByChi2ESProducer",
     chi2cut = cms.double(50.0),
     beamSpot = cms.InputTag("offlineBeamSpot"),
     onlyDiagonal = cms.bool(False)
)
"""

process.TrackAssociatorByHits = cms.ESProducer("TrackAssociatorByHitsESProducer",
    Quality_SimToReco = cms.double(0.5),
    associateRecoTracks = cms.bool(True),
    UseGrouped = cms.bool(True),
    associatePixel = cms.bool(True),
    ROUList = cms.vstring('TrackerHitsTIBLowTof', 
        'TrackerHitsTIBHighTof', 
        'TrackerHitsTIDLowTof', 
        'TrackerHitsTIDHighTof', 
        'TrackerHitsTOBLowTof', 
        'TrackerHitsTOBHighTof', 
        'TrackerHitsTECLowTof', 
        'TrackerHitsTECHighTof', 
        'TrackerHitsPixelBarrelLowTof', 
        'TrackerHitsPixelBarrelHighTof', 
        'TrackerHitsPixelEndcapLowTof', 
        'TrackerHitsPixelEndcapHighTof'),
    UseSplitting = cms.bool(True),
    ComponentName = cms.string('TrackAssociatorByHits'),                                                
    UsePixels = cms.bool(True),
    ThreeHitTracksAreSpecial = cms.bool(True),
    AbsoluteNumberOfHits = cms.bool(False),
    associateStrip = cms.bool(True),
    Purity_SimToReco = cms.double(0.75),
    Cut_RecoToSim = cms.double(0.5),           
    SimToRecoDenominator = cms.string('sim')) ##"reco"

#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source('PoolSource',
#fileNames = cms.untracked.vstring("file:/home/zhangjin/CSCEffStudy/testCSCEff.root")
fileNames = cms.untracked.vstring("rfio:/castor/cern.ch/user/z/zhangjin/testCSCEff.root")
#,eventsToProcess = cms.untracked.VEventRange("1:10925272-1:10925272")
#,skipEvents = cms.untracked.uint32(1700)
#fileNames = cms.untracked.vstring("file:/afs/cern.ch/user/z/zhangjin/scratch0/ProblemEvt.root")
)

process.TracksSelector = cms.EDFilter('CSCPriEff',FileName = cms.untracked.string("TrkCands.root")
,RecoToSimAlgorithm=cms.untracked.string("TrackAssociatorByHits")
)


#process.p2 = cms.Path(process.trackingParticleRecoTrackAsssociation)
process.p = cms.Path(process.TracksSelector)

from FWCore.MessageLogger.MessageLogger_cfi import *
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000   #quench the message logger (optional)
