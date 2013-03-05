########################
# User Specifications  #
########################
RunOnMC=False
AnalyzeZPeak=True#otherwise JPsi
AnalyzeStations=True#otherwise chambers
CalculateSystematic=True#whether calculate the systematic for data; for MC, it will be turned off automatically
DataPileupRootFileName="CMSRun2012C_Cert_201554-202305_8TeV_PromptReco_Collisions12_JSON_MuonPhys.root"
pileup_mc=[2.344E-05,2.344E-05,2.344E-05,2.344E-05,4.687E-04,4.687E-04,7.032E-04,9.414E-04,1.234E-03,1.603E-03,2.464E-03,3.250E-03,5.021E-03,6.644E-03,8.502E-03,1.121E-02,1.518E-02,2.033E-02,2.608E-02,3.171E-02,3.667E-02,4.060E-02,4.338E-02,4.520E-02,4.641E-02,4.735E-02,4.816E-02,4.881E-02,4.917E-02,4.909E-02,4.842E-02,4.707E-02,4.501E-02,4.228E-02,3.896E-02,3.521E-02,3.118E-02,2.702E-02,2.287E-02,1.885E-02,1.508E-02,1.166E-02,8.673E-03,6.190E-03,4.222E-03,2.746E-03,1.698E-03,9.971E-04,5.549E-04,2.924E-04,1.457E-04,6.864E-05,3.054E-05,1.282E-05,5.081E-06,1.898E-06,6.688E-07,2.221E-07,6.947E-08,2.047E-08]
#from SimGeneral.MixingModule.mix_2012_Startup_50ns_PoissonOOTPU_cfi import mix
#pileup_mc=mix.input.nbPileupEvents.probValue

#################################################################################################
# Temporary output files --- find enough space to have them (about the input file size for each)#
#################################################################################################
ClassifiedTreeFile="/uscms/home/zhangjin/nobackup/Select_Chambers.root"# This file includes trees for every station/chamber. The size is about the same with input file size.
TagProbeFitResult="/uscms/home/zhangjin/nobackup/TnP_"#Those files are the TagProbeFitTreeAnalyzer outputs.
ResultPlotsFileName="resultplots.root"#This is the file name of the final result. The path will be Prefix+ResultPlotsFileName. Prefix is the input argument of the Step2_PlotAll.py.
TemporaryOutputFile="/tmp/Tmp_AddedBranchTree.root"#New branches will be added to the orginal tree and saved in this file. The size is about the same with input file size.

###########################################################
# A normal user is not supposed to change anything below  #
###########################################################
#---note: The logical expressions are in C style
DenominatorRequire="( !CSCCBad# && CSCProjDistEdge#<-5 &&  CSCProjDistEdge#> -100 && dRTkMu# < 10. && dRTkMu# > 0.2 && CSCDyProjHVGap#>1. && CSCDyProjHVGap#>3*CSCDyErrProjHVGap# && CSCProjDistEdge#<-3*CSCProjDistErrEdge# )"#default
ProbeSegment="( CSCDxyTTSeg#<40. && CSCDxyTTSeg#>0.)"#default
ProbeLCT="( CSCDxyTTLCT#<40. && CSCDxyTTLCT#>0. )"#default

#DenominatorRequire="( !CSCCBad# && CSCProjDistEdge#<-5 &&  CSCProjDistEdge#> -100 && dRTkMu# < 10. && dRTkMu# > 0.2 ) "#old one
#DenominatorRequire="( !CSCCBad# && CSCProjDistEdge#<-5 &&  CSCProjDistEdge#> -100 && dRTkMu# < 10. && CSCDyProjHVGap#>1. && CSCDyProjHVGap#>3*CSCyErrProjLc# && CSCProjDistEdge#<-3*CSCProjDistErrEdge# ) "# for calculate systematic: close-by muons
#ProbeSegment="( CSCDrTTSeg#<200. && CSCDrTTSeg#>0.)"# for calculate systematic: Segment to TT distance cut
#ProbeLCT="( CSCDrTTLCT#<200. && CSCDrTTLCT#>0. )"# for calculate systematic: LCT to TT distance cut

#ProbeSegment="( CSCDrTTSeg#<40. && CSCDrTTSeg#<5*CSCDrErrTTSeg# )"# too strict cut, not using
#ProbeLCT="( CSCDrTTLCT#<40. && CSCDrTTLCT#<5*CSCDrErrTTLCT# )"# too strict cut, not using

if AnalyzeZPeak:
  MCTruth="( MuTagtracktruth_type==11 && tracktruth_type==11 && ! MuTagtracktruth_isPileup && ! tracktruth_isPileup )"#both from Z
  InvariantMass="( invMass>75. )" #Z peak
else:
  MCTruth="( MuTagtracktruth_type==12 && tracktruth_type==12 && ! MuTagtracktruth_isPileup && ! tracktruth_isPileup )"#both from JPsi
  InvariantMass="(invMass > 2.5 &&  invMass < 3.6)" # JPsi

MuTrackPairCut=InvariantMass+"&&"+"( iSameVtx && minDRHLTAllSingleMu < 0.01 && dRTkMu1 < 10. && dRTkMu1 > 0.2 ) "#default

if RunOnMC:
  CalculateSystematic=False

from  ROOT import *
#define the station loop here, content: (cut,legendname,color,station)
stations={
  1:("(CSCRg1==1 || CSCRg1==4)","ME11",kMagenta,1),
  2:("(CSCRg1==2 || CSCRg1==3)","ME12+13",kBlack,1),
  3:("CSCRg2>0","ME2",kRed,2),
  4:("CSCRg3>0","ME3",kGreen,3),
  5:("CSCRg4>0","ME4",kBlue,4)
  }
n_stations=len(stations)

chambers=[]
for ec_ in (True,False):
  for st_ in (1,2,3,4):
    if st_==1:
      maxrg_=4
    else:
      maxrg_=3

    for rg_ in range(1,maxrg_):
      if st_!=1 and rg_==1:
        chambers_=range(1,19)
      elif st_==4 and rg_==2:
        if ec_:
          chambers_=(9,10,11,12,13)
        else:
          chambers_=()
      else:
        chambers_=range(1,37)
          
      for ch_ in chambers_:
        chambers.append( "ME%s%d_%d_%d"%( '+' if (ec_) else '-', st_, rg_, ch_ ) )

n_chambers=len(chambers)

def Getpuweight(DataPileupRootFileName_=DataPileupRootFileName,pileup_mc_=pileup_mc):
  pileup_file=TFile.Open(DataPileupRootFileName_,"read")
  pileup_data=pileup_file.Get("pileup")
  maxPU=len(pileup_mc_)
  puweight=[0]*maxPU
  renormalize=100/pileup_data.Integral()
  checkplot=TH1F("pileupweight","pileupweight",maxPU,0.,float(maxPU))
  for nPU in range(0,maxPU):
    if pileup_mc_[nPU]>0:
      puweight[nPU]=pileup_data.GetBinContent(pileup_data.FindBin(nPU))/pileup_mc_[nPU]*renormalize
      checkplot.SetBinContent(nPU,puweight[nPU])
  pileup_file.Close()
  return puweight

def ConvertCLogicalExp(Expression_):
  return Expression_.replace("&&"," and ").replace("||"," or ").replace("!="," is not ").replace("!"," not ")
