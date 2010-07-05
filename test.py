import FWCore.ParameterSet.Config as cms
process = cms.Process("PROC")

#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
#process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
#process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cff")


#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source('PoolSource',
#fileNames = cms.untracked.vstring("file:/home/zhangjin/CSCEffStudy/testCSCEff.root")
fileNames = cms.untracked.vstring("file:/afs/cern.ch/user/z/zhangjin/scratch0/testCSCEff.root")
)

process.TrackAssociatorByChi2ESProducer.chi2cut=cms.double(25.0)

process.TracksSelector = cms.EDFilter('CSCPriEff',
FileName = cms.untracked.string("TrkCands.root"),
LocalRun = cms.untracked.bool(True),
MaxDR = cms.untracked.double(0.15),
MaxRelpT = cms.untracked.double(0.2)
)


#process.p2 = cms.Path(process.trackingParticleRecoTrackAsssociation)
process.p = cms.Path(process.TracksSelector)

from FWCore.MessageLogger.MessageLogger_cfi import *
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000   #quench the message logger (optional)
