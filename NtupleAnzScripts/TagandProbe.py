#################################################
#             Tree Fitter                       #
#################################################
#Author: Jinzhong Zhang(zhangjin@cern.ch), Northeastern University, Boston, U.S.A
#Based on: Chaouki's old tree fitter settings
#This script should be run by SelectAndFit.py
#Usage: cmsRun tagandProbe.py File Station IsZPeak(ture/false)
#Station=[0-4] station 0 means "or" combination of several stations
import sys

if sys.argv[0] == "python" or sys.argv[0] == "cmsRun":
    args=sys.argv[1:]
else:
    args=sys.argv[0:]

from Config import TreeName

if len(args)>1:
    ClassifiedTreeFile=args[1]

if len(args)>2:
    from Config import ptbin,etabin1,etabin2,etabin3,etabin4,phibin  #import the pt,eta,phi binning
    station=int(args[2])
    if station<=1:
        etabin=etabin1
    elif station==2:
        etabin=etabin2
    elif station==3:
        etabin=etabin3
    elif station==4:
        etabin=etabin4
    else:
        print "No input file. Exit."
        sys.exit()
else:
    print "No input station number. Exit."
    sys.exit()

if len(args)>3:
    IsZPeak=args[3]
else:
    from Config import Resonance
    IsZPeak=Resonance is "Z"

from Config import CalculateSystematic as CalculateSystematic
from Config import RunOnMC as RunOnMC
from Config import TemporaryOutputFile
print "Running TagProbeFitTreeAnalyzer on ","station%d on"%(station), "Z" if IsZPeak else "JPsi","pole."


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
    OutputFileName     = cms.string("%s%s" % (TagProbeFitResult,ClassifiedTreeFile.replace(TemporaryOutputFile[:-5],"")) ), 

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
#       tracks_eta = cms.vstring("Probe #eta", "-2.5", "2.5", ""),
        abseta     = cms.vstring("Probe |#eta|", "0.", "2.5", ""),
        shiftedphi = cms.vstring("Probe #phi", "-0.0872664626", "6.28318531", ""),
#       tracks_phi = cms.vstring("Probe #phi", "0.", "6.32", ""),
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

EfficiencyBins = cms.PSet(
            tracks_pt        = cms.vdouble(ptbin),
            abseta        = cms.vdouble(etabin),
            #tracks_phi       = cms.vdouble(phibin)
            shiftedphi = cms.vdouble(phibin)
)
#the name of the parameter set becomes the name of the directory
process.TagProbeFitTreeAnalyzer.Efficiencies = cms.PSet(
    lct_effV = cms.PSet(
        EfficiencyCategoryAndState = cms.vstring("foundLCTSt%d"%(station),"true"),
        UnbinnedVariables = cms.vstring("invMass","weight"),
        BinnedVariables = cms.PSet(EfficiencyBins),
        BinToPDFmap = cms.vstring("VoigtianPlusExpo") if IsZPeak else cms.vstring("VoigtianPlusLinear")
        ),
    seg_effV = cms.PSet(
        EfficiencyCategoryAndState = cms.vstring("foundSEGSt%d"%(station),"true"),
        UnbinnedVariables = cms.vstring("invMass","weight"),
        BinnedVariables = cms.PSet(EfficiencyBins),
        BinToPDFmap = cms.vstring("VoigtianPlusExpo") if IsZPeak else cms.vstring("VoigtianPlusLinear")
        )
)

if CalculateSystematic and (not RunOnMC):
    print "I am going to calculate the effficency with different signal and background modeling..."
    process.TagProbeFitTreeAnalyzer.Efficiencies = cms.PSet(
        process.TagProbeFitTreeAnalyzer.Efficiencies,
        lct_effV_SigModeling = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("foundLCTSt%d"%(station),"true"),
            UnbinnedVariables = cms.vstring("invMass"),
            BinnedVariables = cms.PSet(EfficiencyBins),
            BinToPDFmap = cms.vstring("GaussianPlusExpo") if IsZPeak else cms.vstring("GaussianPlusLinear")
            ),
        seg_effV_SigModeling = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("foundSEGSt%d"%(station),"true"),
            UnbinnedVariables = cms.vstring("invMass"),
            BinnedVariables = cms.PSet(EfficiencyBins),
            BinToPDFmap = cms.vstring("GaussianPlusExpo") if IsZPeak else cms.vstring("GaussianPlusLinear")
            ),
        lct_effV_BkgModeling = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("foundLCTSt%d"%(station),"true"),
            UnbinnedVariables = cms.vstring("invMass"),
            BinnedVariables = cms.PSet(EfficiencyBins),
            BinToPDFmap = cms.vstring("VoigtianPlusQuadratic")
            ),
        seg_effV_BkgModeling = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("foundSEGSt%d"%(station),"true"),
            UnbinnedVariables = cms.vstring("invMass"),
            BinnedVariables = cms.PSet(EfficiencyBins),
            BinToPDFmap = cms.vstring("VoigtianPlusQuadratic")
            )
        )

if RunOnMC: # BinToPDFmap has to be defined, otherwise it will crash. Problems are in TagProbeFitTreeAnalyzer. Only the CntEfficiency (counting) will be used later.
    process.TagProbeFitTreeAnalyzer.Categories = cms.PSet(
        process.TagProbeFitTreeAnalyzer.Categories,
        mcTrue = cms.vstring("mcTrue", "dummy[true=1,false=0]")
        )
    process.TagProbeFitTreeAnalyzer.Efficiencies = cms.PSet(
        process.TagProbeFitTreeAnalyzer.Efficiencies,
        lct_effV_MCTruth = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("foundLCTSt%d"%(station),"true"),
            UnbinnedVariables = cms.vstring("invMass","weight"),
            BinnedVariables = cms.PSet(EfficiencyBins,mcTrue = cms.vstring("true")),
            BinToPDFmap = cms.vstring("VoigtianPlusExpo") if IsZPeak else cms.vstring("VoigtianPlusLinear")
            ),
        seg_effV_MCTruth = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("foundSEGSt%d"%(station),"true"),
            UnbinnedVariables = cms.vstring("invMass","weight"),
            BinnedVariables = cms.PSet(EfficiencyBins,mcTrue = cms.vstring("true")),
            BinToPDFmap = cms.vstring("VoigtianPlusExpo") if IsZPeak else cms.vstring("VoigtianPlusLinear")
            )
        )

print process.TagProbeFitTreeAnalyzer.PDFs.VoigtianPlusExpo
print process.TagProbeFitTreeAnalyzer.Efficiencies

process.fitness = cms.Path(
    process.TagProbeFitTreeAnalyzer
)
