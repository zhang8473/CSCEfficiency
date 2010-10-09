import FWCore.ParameterSet.Config as cms
process = cms.Process("TestParticle")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

process.source = cms.Source('PoolSource',
fileNames = cms.untracked.vstring("file:/home/zhangjin/CSCEffStudy/testCSCEff.root")
)
#process.source = cms.Source('PoolSource',
#fileNames = cms.untracked.vstring("file:/home/zhangjin/Test_CMSSW333_DYmumu_Mcut500_7TeV_n10_MC31XV9_RECO.root")
#)

process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
    src = cms.InputTag("genParticles"),
    maxEventsToPrint  = cms.untracked.int32(1),
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(False),
    printVertex = cms.untracked.bool(False),
    printStatus = cms.untracked.bool(False),
    printIndex  = cms.untracked.bool(False)
)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printList = cms.EDAnalyzer("ParticleListDrawer",
  maxEventsToPrint  = cms.untracked.int32(10),
  printVertex = cms.untracked.bool(False),
  src = cms.InputTag("genParticles")
)


process.printEventNumber = cms.OutputModule("AsciiOutputModule")

#process.p1 = cms.Path(process.printTree)
process.p2 = cms.Path(process.printList)

process.outpath = cms.EndPath(process.printEventNumber)
process.MessageLogger.destinations = cms.untracked.vstring('cout','cerr')

