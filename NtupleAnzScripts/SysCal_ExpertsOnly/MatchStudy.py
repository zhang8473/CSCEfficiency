#!/usr/bin/python
#######################################################
# Tracker track to LCT/Segments distance study script #
# Author: Zhang Jinzhong, zhangjin@cern.ch            #
#######################################################
#Usage: python MatchStudy.py TheNtupleRootFile.root
from  ROOT import *
from Config import *
import sys

gROOT.SetStyle("Plain")
gStyle.SetOptStat(0)

if sys.argv[0] == "python":
    args=sys.argv[2:]
else:
    args=sys.argv[1:]

def AddWeightToTree(tree_,sampleweight=1.):
  pileup_file=TFile.Open(DataPileupRootFileName,"read")
  pileup_data=pileup_file.Get("pileup")
  maxPU=len(pileup_mc)
  puweight=[0]*maxPU
  renormalize=100/pileup_data.Integral()
  for nPU in range(0,maxPU):
    if pileup_mc[nPU]>0:
      puweight[nPU]=pileup_data.GetBinContent(pileup_data.FindBin(nPU))/pileup_mc[nPU]*renormalize
  pileup_file.Close()
  tree_cloned_=tree_.CloneTree(0)
# create 1 dimensional float arrays (python's float datatype corresponds to c++ doubles) as fill variables
  weight = zeros(1, dtype=float)
  weightBranch = tree_cloned_.Branch('weight',weight, 'weight/D')
  for n in range(tree_.GetEntries()):
    tree_.GetEntry(n)
    weight[0] = tree_.mcweight*puweight[int(tree_.numberOfPUVerticesMixingTruth)]*sampleweight
    tree_cloned_.Fill()
  return tree_cloned_

file_in = TFile.Open(args[0],"read")
tree_in = file_in.Get("Fraction")
if RunOnMC:
    tmpoutputfile=TFile.Open(TemporaryOutputFile,'RECREATE')
    tree_in = AddWeightToTree(tree_in)
    tree_in.AutoSave()

#content order: (formula,x_min,x_max,y_max)
plots={
       "TTtoSegX":('CSCDxTTSeg#',-100.,100.),
       "TTtoSegY":('CSCDyTTSeg#',-100.,100.),
#       "TTtoSegXY":('CSCDxyTTSeg#',-1.,150.),
       "TTtoLCTX":('CSCDxTTLCT#',-100.,100.),
       "TTtoLCTY":('CSCDyTTLCT#',-100.,100.),
#       "TTtoLCTXY":('CSCDxyTTLCT#',0.,150.),
       "TTtoSegXY":('sqrt(pow(CSCDxTTSeg#,2)+pow(CSCDyTTSeg#,2))',-1.,150.),# should be the same with TTtoSegR
       "TTtoLCTXY":('sqrt(pow(CSCDxTTLCT#,2)+pow(CSCDyTTLCT#,2))',-1.,150.)# should be the same with TTtoLCTR
#       "SegtoLCTX":('CSCSegxLc#-CSCLCTxLc#',-20.,20.),
#       "SegtoLCTY":('CSCSegyLc#-CSCLCTyLc#',-60.,60.),
#       "SegtoLCTR":('sqrt(pow(CSCSegxLc#-CSCLCTxLc#,2)+pow(CSCSegyLc#-CSCLCTyLc#,2))',-1.,40.,30.)
       }

def MakeAPlot(Name_,Forumla_,Cut_,color_):
    if RunOnMC:
        tree_in.Project(Name_,Forumla_,'weight*(%s)'%(Cut_))
    else:
        tree_in.Project(Name_,Forumla_,Cut_)
    exec( Name_+'.Scale( 100/'+Name_+'.Integral() )' )
    exec( Name_+'.SetLineColor(color_)' )
    exec( Name_+'.SetMarkerColor(color_)' )
    exec( 'result=(%s.GetMean(),%s.GetRMS())'%(Name_,Name_) )
    print Name_,result
    return result
    
for plot in plots:
    exec( plot+'=TCanvas("%s","%s",500,500)'%(plot,plot) )
    exec( plot+"_Max=0" )
    exec( 'legend_%s=TLegend(0.6,0.6,0.8,0.8)'%(plot) )

for idx in range(1,n_stations+1):
    Cut=stations[idx][0]+" && "+DenominatorRequire.replace("#",str(stations[idx][3]))
    for plot in plots:
        Name='%s_%d'%(plot,idx)
        exec( Name+'=TH1D(Name,plot+";cm;scaled number of entries",100,%f,%f)'%(plots[plot][1],plots[plot][2]) )
        MeanAndRMS=MakeAPlot(Name,plots[plot][0].replace('#',str(stations[idx][3])),Cut,stations[idx][2])
        exec( "max="+Name+".GetBinContent("+Name+".GetMaximumBin())" )
        exec( "if "+plot+"_Max<max: "+plot+"_Max=max" )
        exec( 'legend_%s.AddEntry('%(plot)+Name+',\"'+stations[idx][1]+': mean:%.2fcm;RMS:%.2fcm\",\"l\")'%(MeanAndRMS[0],MeanAndRMS[1]) )

for plot in plots:
    exec( plot+'.cd()')
    DrawOption="H"
    for idx in range(1,n_stations+1):
        Name='%s_%d'%(plot,idx)
        exec( Name+'.Draw(\"'+DrawOption+'\")' )
        if DrawOption is "H":
            exec( Name+'.SetMaximum(%s_Max*1.05)'%(plot) )
            DrawOption="H,same"
        exec( 'legend_%s.Draw()'%(plot) )
        exec( plot+'.SetLogy()')

raw_input("Please adjust plots. Press ENTER when finish.")

for plot in plots:
  exec( plot+'.cd()')
  exec( plot+'.SaveAs("%s.pdf")'%(plot) )

file_in.Close()
if RunOnMC:
    tmpoutputfile.Close()
    import os
    os.system( "rm "+TemporaryOutputFile )
#print "It will probably crash while running sys.exit(). Just ignore that."
