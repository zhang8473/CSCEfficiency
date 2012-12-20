#!/usr/bin/python
#Author: Jinzhong Zhang(zhangjin@cern.ch), Northeastern University, Boston, U.S.A
#Usage: python Step1.py haddedrootfile.root (no more option)

from  ROOT import *
from  numpy import *
from Config import *

import sys,os
if (sys.argv[0] == "python"): args=sys.argv[2:]
else: args=sys.argv[1:]

if len(args)<2:
  Fit_="fit"
else:
  Fit_=args[1]

if RunOnMC:
  puweight=Getpuweight(DataPileupRootFileName,pileup_mc)

def CalculateOverallScaleNum(tree_):
  overall=0.
  for n in range(tree_.GetEntries()):
    tree_.GetEntry(n)
    overall+=tree_.mcweight*puweight[int(tree_.numberOfPUVerticesMixingTruth)]
  return tree_.GetEntries()/overall

def AddNewBranchesToTree(tree_,sampleweight=1.,segSelect_=[],lctSelect_=[]):
# create 1 dimensional float arrays (python's float datatype corresponds to c++ doubles) as fill variables
  tree_cloned_=tree_.CloneTree(0)
  print "Now adding:\n abseta,\n weight",
  abseta = zeros(1, dtype=float)
  absetaBranch = tree_cloned_.Branch('abseta',abseta, 'abseta/D')
  weight = zeros(1, dtype=float)
  weightBranch = tree_cloned_.Branch('weight',weight, 'weight/D')
  if RunOnMC:
    print ",\n isTrueMuMuPair:"
    MCTruth_=ConvertCLogicalExp(MCTruth).replace("MuTag"," tree_.MuTag").replace(" track"," tree_.track")
    print "\t MC truth criteria: \n",MCTruth_,","
    isTrueMuMuPair = zeros(1, dtype=int)
    isTrueMuMuPairBranch = tree_cloned_.Branch('isTrueMuMuPair',isTrueMuMuPair, 'isTrueMuMuPair/I')
  else:
    print "(=1),"
  if len(segSelect_)>0:
    for st_ in (1,2,3,4):
      exec( "foundSEGSt%d = zeros(1, dtype=int)"%(st_) )
      exec( "foundLCTSt%d = zeros(1, dtype=int)"%(st_) )
      exec( "segBranch%d = tree_cloned_.Branch('foundSEGSt%d',foundSEGSt%d, 'foundSEGSt%d/I' )"%(st_,st_,st_,st_) )
      exec( "lctBranch%d = tree_cloned_.Branch('foundLCTSt%d',foundLCTSt%d, 'foundLCTSt%d/I' )"%(st_,st_,st_,st_) )
    print " foundSEGSt1-4 (seg probe),\n foundLCTSt1-4 (lct probe)"
  print "to the tree.\nPlease wait about half to more than one hours......"

  for n in range(tree_.GetEntries()):
    tree_.GetEntry(n)
    abseta[0] = abs(tree_.tracks_eta)
    if RunOnMC:
      weight[0] = tree_.mcweight*puweight[int(tree_.numberOfPUVerticesMixingTruth)]*sampleweight
      exec ( "isTrueMuMuPair[0] = 1 if %s else 0"%(MCTruth_) )
    else:
      weight[0] = 1.
    if len(segSelect_)>0:
      for st_ in (1,2,3,4):
        exec( "foundSEGSt%d[0] = 1 if %s else 0"%(st_,segSelect_[st_-1]) )
        exec( "foundLCTSt%d[0] = 1 if %s else 0"%(st_,lctSelect_[st_-1]) )
    tree_cloned_.Fill()
  return tree_cloned_

segSelect=[""]*4
lctSelect=[""]*4
for st_ in (1,2,3,4):
  segSelect[st_-1] = ConvertCLogicalExp(ProbeSegment).replace("CSC","tree_.CSC").replace("#",str(st_))
  lctSelect[st_-1] = ConvertCLogicalExp(ProbeLCT).replace("CSC","tree_.CSC").replace("#",str(st_))
print "segment selection criteria: \n",segSelect
print "local charged trigger selection criteria: \n",lctSelect

file_in = TFile.Open(args[0],"read")
tmpoutputfile=TFile.Open(TemporaryOutputFile,'RECREATE')
tree_in = file_in.Get("Fraction")
if RunOnMC:
  overallweight=CalculateOverallScaleNum(tree_in)
else:
  overallweight=1.
tree_tmpout=AddNewBranchesToTree( tree_in, overallweight, segSelect, lctSelect )
tree_tmpout.AutoSave()
print "Output to ", tmpoutputfile.GetName()     
tmpoutputfile.Close()
file_in.Close()

if AnalyzeStations:
  os.system( "python SelectAndFit.py "+Fit_+" &" )
  print "\"nohup python SelectAndFit.py "+Fit_+" &\" is running..."
else:
  job_scheduler=[(0,99),(100,199),(200,299),(300,399),(400,n_chambers-1)]
  for ajob in job_scheduler:
    command="nohup python SelectAndFit.py "+Fit_+" "+TemporaryOutputFile+" "+str(ajob[0])+","+str(ajob[1])+" &"
    os.system( command )
    print "\"%s\" is running..." % ( command )
    
