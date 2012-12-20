#################################################
#             Tree Fitter                       #
#################################################
#Author: Jinzhong Zhang(zhangjin@cern.ch), Northeastern University, Boston, U.S.A
#Based on: Chaouki's old tree fitter settings
#This script should be run by SelectAndFit.py
#Usage: cmsRun tagandProbe.py TreeName(string, e.g. ME-1_1_3 ) ClassifiedTreeFile(optional) IsZPeak(ture/false)
import sys

if sys.argv[0] == "python" or sys.argv[0] == "cmsRun":
    args=sys.argv[1:]
else:
    args=sys.argv[0:]

TreeName=args[1]

if len(args)>2:
    ClassifiedTreeFile=args[2]
else:
    from Config import ClassifiedTreeFile as ClassifiedTreeFile

if len(args)>3:
    IsZPeak=args[3]
else:
    from Config import AnalyzeZPeak as IsZPeak

from Config import CalculateSystematic as CalculateSystematic
from Config import RunOnMC as RunOnMC

if TreeName[2] in ("1","2","3","4"):
    station=int(TreeName[2])
elif TreeName[3] in ("1","2","3","4"):
    station=int(TreeName[3])
else:
    print "Incorrect TreeName: cannot parse the station number from TreeName. Is that MEx or ME+x, or ME-x?"
    sys.exit()

print "Running TagProbeFitTreeAnalyzer on "+TreeName+", station%d...."%(station),"on","Z" if IsZPeak else "JPsi"

import FWCore.ParameterSet.Config as cms
process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

from Config import TagProbeFitResult as TagProbeFitResult
process.TagProbeFitTreeAnalyzer = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    # IO parameters:
    InputFileNames     = cms.vstring(ClassifiedTreeFile),
    InputDirectoryName = cms.string("Probes"),
    InputTreeName      = cms.string(TreeName),
    OutputFileName     = cms.string("%s%s.root" % (TagProbeFitResult,TreeName) ), 

    #numbrer of CPUs to use for fitting
    NumCPU = cms.uint32(3),
    # specifies wether to save the RooWorkspace containing the data for each bin and
    # the pdf object with the initial and final state snapshots
    SaveWorkspace = cms.bool(True),
    binsForMassPlots = cms.uint32(40),
    # defines all the real variables of the probes available in the input tree and intended for use in the efficiencies
    Variables = cms.PSet(
        invMass    = cms.vstring("Tag-Probe Mass", "75", "110", "GeV/c^{2}") if IsZPeak else cms.vstring("Tag-Probe Mass", "2.5", "3.6", "GeV/c^{2}"),
        tracks_pt  = cms.vstring("Probe p_{T}", "0", "1000", "GeV/c"),
        tracks_eta = cms.vstring("Probe #eta", "-2.5", "2.5", ""),
#        abseta     = cms.vstring("Probe #eta", "0.", "4.8", ""),
        tracks_phi = cms.vstring("Probe #phi", "0.", "6.32", ""),
        tracks_e   = cms.vstring("Probe p", "1.", "1000.", "GeV/c"),
        weight = cms.vstring("weight","0","1000.","")
    ),
    WeightVariable = cms.string("weight"),
    # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculations
    Categories = cms.PSet(), #will be defined later
    # defines all the PDFs that will be available for the efficiency calculations; uses RooFit's "factory" syntax;
    # each pdf needs to define "signal", "backgroundPass", "backgroundFail" pdfs, "efficiency[0.9,0,1]" and "signalFractionInPassing[0.9]" are used for initial values  
    PDFs = cms.PSet(
        GaussianPlusExpo = cms.vstring(
            "Gaussian::signal(invMass, mean[91.2, 88.70, 93.70], sigma[2.3, 0.5, 10.0])" if IsZPeak else "Gaussian::signal(invMass, mean[3.10, 3.0, 3.2], sigma[0.03,0.01,0.05])",
            "Exponential::backgroundPass(invMass, lp[0,-5,5])",
            "Exponential::backgroundFail(invMass, lf[0,-5,5])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
            ),
        VoigtianPlusQuadratic = cms.vstring(
            "Voigtian::signal(invMass, mean[91.2, 88.70, 93.70], width[2.495], sigma[2.3, 0.5, 10.0])" if IsZPeak else "Voigtian::signal(invMass, mean[3.10, 3.0, 3.2], width[0.0929,0.09,0.96], sigma[0.03,0.01,0.05])",
            "Chebychev::backgroundPass(invMass, {cPass1[0,-1,1], cPass2[0,-1,1]})",
            "Chebychev::backgroundFail(invMass, {cFail1[0,-1,1], cFail2[0,-1,1]})",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
            ),
        VoigtianPlusExpo = cms.vstring(
            "Voigtian::signal(invMass, mean[91.2, 88.70, 93.70], width[2.495], sigma[2.3, 0.5, 10.0])" if IsZPeak else "Voigtian::signal(invMass, mean[3.10, 3.0, 3.2], width[0.0929,0.09,0.96], sigma[0.03,0.01,0.05])",
            "Exponential::backgroundPass(invMass, lp[0,-5,5])",
            "Exponential::backgroundFail(invMass, lf[0,-5,5])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
            ),
        GaussianPlusQuadratic = cms.vstring(
            "Gaussian::signal(invMass, mean[91.2, 88.70, 93.70], sigma[2.3, 0.5, 10.0])" if IsZPeak else "Gaussian::signal(invMass, mean[3.10, 3.0, 3.2], sigma[0.03,0.01,0.05])",
            "Chebychev::backgroundPass(invMass, {cPass1[0,-1,1], cPass2[0,-1,1]})",
            "Chebychev::backgroundFail(invMass, {cFail1[0,-1,1], cFail2[0,-1,1]})",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
            ),
        GaussianPlusLinear = cms.vstring(
            "Gaussian::signal(invMass, mean[3.1,3.0,3.2], sigma[0.03,0.01,0.05])",
            "Chebychev::backgroundPass(invMass, cPass[0,-1,1])",
            "Chebychev::backgroundFail(invMass, cFail[0,-1,1])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
            ),
        VoigtianPlusLinear = cms.vstring(
            "Voigtian::signal(invMass, mean[91.2, 88.70, 93.70], width[2.495], sigma[2.3, 0.5, 10.0])" if IsZPeak else "Voigtian::signal(invMass, mean[3.10, 3.0, 3.2], width[0.0929,0.09,0.96], sigma[0.03,0.01,0.05])",
            "Chebychev::backgroundPass(invMass, cPass[0,-1,1])",
            "Chebychev::backgroundFail(invMass, cFail[0,-1,1])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
            )
        ),
    Efficiencies = cms.PSet() # will be filled later
    # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
    # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting. 
)

exec( "process.TagProbeFitTreeAnalyzer.Categories = cms.PSet( foundLCTSt%d  = cms.vstring(\"lctStation%d\",  \"dummy[true=1,false=0]\"), foundSEGSt%d  = cms.vstring(\"segStation%d\",  \"dummy[true=1,false=0]\") )"%(station,station,station,station) )

#the name of the parameter set becomes the name of the directory
process.TagProbeFitTreeAnalyzer.Efficiencies = cms.PSet(
    lct_effV = cms.PSet(
        EfficiencyCategoryAndState = cms.vstring("foundLCTSt%d"%(station),"true"),
        UnbinnedVariables = cms.vstring("invMass","weight"),
        BinnedVariables = cms.PSet(
            tracks_pt        = cms.vdouble(0., 1000.),
            tracks_e        = cms.vdouble(15., 1000.),
            tracks_phi       = cms.vdouble(-7.0, 7.0),
            ),
        BinToPDFmap = cms.vstring("VoigtianPlusExpo") if IsZPeak else cms.vstring("VoigtianPlusLinear")
        ),
    seg_effV = cms.PSet(
        EfficiencyCategoryAndState = cms.vstring("foundSEGSt%d"%(station),"true"),
        UnbinnedVariables = cms.vstring("invMass","weight"),
        BinnedVariables = cms.PSet(
            tracks_pt        = cms.vdouble(0., 1000.),
            tracks_e        = cms.vdouble(15., 1000.),
            tracks_phi       = cms.vdouble(-7.0, 7.0),
            ),
        BinToPDFmap = cms.vstring("VoigtianPlusExpo") if IsZPeak else cms.vstring("VoigtianPlusLinear")
        )
)

print process.TagProbeFitTreeAnalyzer.PDFs.VoigtianPlusExpo

if CalculateSystematic and (not RunOnMC):
    print "I am going to calculate the effficency with different signal and background modeling..."
    process.TagProbeFitTreeAnalyzer.Efficiencies = cms.PSet(
        process.TagProbeFitTreeAnalyzer.Efficiencies,
        lct_effV_SigModeling = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("foundLCTSt%d"%(station),"true"),
            UnbinnedVariables = cms.vstring("invMass"),
            BinnedVariables = cms.PSet(
                tracks_pt        = cms.vdouble(0., 1000.),
                tracks_e        = cms.vdouble(15., 1000.),
                tracks_phi       = cms.vdouble(-7.0, 7.0),
                ),
            BinToPDFmap = cms.vstring("GaussianPlusExpo") if IsZPeak else cms.vstring("GaussianPlusLinear")
            ),
        seg_effV_SigModeling = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("foundSEGSt%d"%(station),"true"),
            UnbinnedVariables = cms.vstring("invMass"),
            BinnedVariables = cms.PSet(
                tracks_pt        = cms.vdouble(0., 1000.),
                tracks_e        = cms.vdouble(15., 1000.),
                tracks_phi       = cms.vdouble(-7.0, 7.0),
                ),
            BinToPDFmap = cms.vstring("GaussianPlusExpo") if IsZPeak else cms.vstring("GaussianPlusLinear")
            ),
        lct_effV_BkgModeling = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("foundLCTSt%d"%(station),"true"),
            UnbinnedVariables = cms.vstring("invMass"),
            BinnedVariables = cms.PSet(
                tracks_pt        = cms.vdouble(0., 1000.),
                tracks_e        = cms.vdouble(15., 1000.),
                tracks_phi       = cms.vdouble(-7.0, 7.0),
                ),
            BinToPDFmap = cms.vstring("VoigtianPlusQuadratic")
            ),
        seg_effV_BkgModeling = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("foundSEGSt%d"%(station),"true"),
            UnbinnedVariables = cms.vstring("invMass"),
            BinnedVariables = cms.PSet(
                tracks_pt        = cms.vdouble(0., 1000.),
                tracks_e        = cms.vdouble(15., 1000.),
                tracks_phi       = cms.vdouble(-7.0, 7.0),
                ),
            BinToPDFmap = cms.vstring("VoigtianPlusQuadratic")
            )
        )
"""
if RunOnMC: # it will crash. Problems are in TagProbeFitTreeAnalyzer.
    process.TagProbeFitTreeAnalyzer.Categories = cms.PSet(
        process.TagProbeFitTreeAnalyzer.Categories,
        isTrueMuMuPair = cms.vstring("MCtrue", "dummy[true=1,false=0]")
        )
    process.TagProbeFitTreeAnalyzer.Efficiencies = cms.PSet(
        process.TagProbeFitTreeAnalyzer.Efficiencies,
        lct_effV_MCTruth = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("foundLCTSt%d"%(station),"true"),
            UnbinnedVariables = cms.vstring("invMass","weight"),
            BinnedVariables = cms.PSet(
                isTrueMuMuPair = cms.vstring("true"),
                tracks_pt        = cms.vdouble(0., 1000.),
                tracks_e        = cms.vdouble(15., 1000.),
                tracks_phi       = cms.vdouble(-7.0, 7.0),
                )
            ),
        seg_effV_MCTruth = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("foundSEGSt%d"%(station),"true"),
            UnbinnedVariables = cms.vstring("invMass","weight"),
            BinnedVariables = cms.PSet(
                isTrueMuMuPair = cms.vstring("true"),
                tracks_pt        = cms.vdouble(0., 1000.),
                tracks_e        = cms.vdouble(15., 1000.),
                tracks_phi       = cms.vdouble(-7.0, 7.0),
                )
            )
        )
"""

process.fitness = cms.Path(
    process.TagProbeFitTreeAnalyzer
)
