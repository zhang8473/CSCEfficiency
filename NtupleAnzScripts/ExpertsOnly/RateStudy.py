#!/usr/bin/python
#######################################################
# Tracker track to LCT/Segments distance study script #
# Author: Zhang Jinzhong, zhangjin@cern.ch            #
#######################################################
#Usage: python MatchStudy.py TheNtupleRootFile.root
from  ROOT import *
from Config import *
from  numpy import *
import sys

ReweightMC=False
RunOnMC=True
gROOT.SetStyle("Plain")
gStyle.SetOptStat(0)
plotmin=0
plotmax=1.0

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
  tmpoutputfile=TFile.Open(TemporaryOutputFile,'RECREATE')
  tmpoutputfile.cd()
  tree_cloned_=tree_.CloneTree(0)
# create 1 dimensional float arrays (python's float datatype corresponds to c++ doubles) as fill variables
  weight = zeros(1, dtype=float)
  weightBranch = tree_cloned_.Branch('weight',weight, 'weight/D')
  for n in range(tree_.GetEntries()):
    tree_.GetEntry(n)
    weight[0] = tree_.cweight*puweight[int(tree_.numberOfPUVerticesMixingTruth)]*sampleweight
    tree_cloned_.Fill()
  tree_cloned_.AutoSave()
  tmpoutputfile.Close()
  return tree_cloned_

file_in = TFile.Open(args[0],"read")
tree_in = file_in.Get("Fraction")
if ReweightMC:
    tree_in = AddWeightToTree(tree_in)
    print "Please set ReweightMC=False and run \"python RateStudy.py\"",TemporaryOutputFile
    sys.exit()
#content order: name:(formula,x_min,x_max,y_max)

plots={
#    "tracks_IsoR03Ratio":('tracks_IsoR03Ratio',10,0,1.,DenominatorRequire+"&&"+MCTruth,DenominatorRequire+"&&"+MCTruth+"&& ("+ProbeLCT+")"),
    "tracks_ptRevUncertainty":('tracks_ptError/tracks_pt',10,0,0.1,DenominatorRequire+"&&"+MCTruth,DenominatorRequire+"&&"+MCTruth+"&& ("+ProbeLCT+")"),
#    "tracks_dxy":('tracks_dxy',10,0,1.,DenominatorRequire+"&&"+MCTruth,DenominatorRequire+"&&"+MCTruth+"&& ("+ProbeLCT+")"),
#    "TracktoLCTX_in_cm":('CSCDxTTLCT#',12,0,40.,DenominatorRequire+"&&"+MCTruth,DenominatorRequire+"&&"+MCTruth+"&& ("+ProbeLCT+")"),
#    "TracktoLCTY_in_cm":('CSCDyTTLCT#',12,0,40.,DenominatorRequire+"&&"+MCTruth,DenominatorRequire+"&&"+MCTruth+"&& ("+ProbeLCT+")"),
#    "SegmenttoLCTX_in_cm":('CSCDxTTSeg#',12,0,40.,DenominatorRequire+"&&"+MCTruth,DenominatorRequire+"&&"+MCTruth+"&& ("+ProbeSegment+")"),
#    "SegmenttoLCTY_in_cm":('CSCDyTTSeg#',12,0,40.,DenominatorRequire+"&&"+MCTruth,DenominatorRequire+"&&"+MCTruth+"&& ("+ProbeSegment+")"),
#    "DisttoEdge_in_cm":('CSCProjDistEdge#',10,-25,0.,DenominatorRequire+"&&"+MCTruth,DenominatorRequire+"&&"+MCTruth+"&& ("+ProbeLCT+")"),
#    "abseta":('abs(tracks_eta)',10,1.5,2.5,DenominatorRequire+"&&"+MCTruth,DenominatorRequire+"&&"+MCTruth+"&& ("+ProbeLCT+")")
    }

def MakeARatioPlot(Name_,plot_,st_,tree_in):
    Name_Num=Name_+"_Num"
    Name_Den=Name_+"_Den"
    Forumla_=plot_[0].replace('#',str(st_[3]))
    Den_=st_[0]+" && "+plot_[4].replace("#",str(st_[3]))
    Num_=st_[0]+" && "+plot_[5].replace("#",str(st_[3]))
    exec( Name_Num+'=TH1D(Name+"_Num",plot+";cm;efficiency",%d,%f,%f)'%(plot_[1],plot_[2],plot_[3]) )
    exec( Name_Den+'=TH1D(Name+"_Den",plot+";cm;number of entires",%d,%f,%f)'%(plot_[1],plot_[2],plot_[3]) )
    tree_in.Project(Name_Den,Forumla_,'weight*(%s)'%(Den_) if RunOnMC else Den_ )
    tree_in.Project(Name_Num,Forumla_,'weight*(%s)'%(Num_) if RunOnMC else Num_ )
    exec( 'NEntries='+Name_Den+'.GetEntries()')
    if NEntries>0:
        exec( "eff=TEfficiency("+Name_Num+","+Name_Den+")" )
        exec( "eff.SetName(\""+Name_+"\")" )
        eff.SetLineColor(st_[2])
        eff.SetMarkerColor(st_[2])
        print Name_,"eff=",
        for i in range(plot_[1]):
            print "%.3f"%(eff.GetEfficiency(i)),
        print
        return eff
    else:
        print "Skip "+Name_+", denominator is 0."
        return NULL
    
for plot in plots:
    exec( plot+'=TCanvas("%s","%s",500,500)'%(plot,plot) )
    exec( plot+"_Max=0" )
    exec( 'legend_%s=TLegend(0.6,0.6,0.8,0.8)'%(plot) )
    exec( 'legend_%s.SetFillColor(0)'%(plot) )

for idx in range(1,n_stations+1):
    for plot in plots:
        Name='%s_%d'%(plot,idx)
        exec( Name+"=MakeARatioPlot(Name,plots[plot],stations[idx],tree_in)" )
        exec( 'legend_%s.AddEntry('%(plot)+Name+',\"'+stations[idx][1]+'","l")' )

for plot in plots:
    exec( plot+'.cd()')
    exec( plot+'_dummy=TH1D(plot+"_dummy",plot+";"+plot.replace("_"," ")+";efficiency",plots[plot][1],plots[plot][2],plots[plot][3])' )
    exec( plot+'_dummy.Draw("")' )
    exec( plot+'_dummy.SetMinimum(plotmin)' )
    exec( plot+'_dummy.SetMaximum(plotmax)' )
    for idx in range(n_stations,0,-1):
        Name='%s_%d'%(plot,idx)
        exec( Name+'.Draw("P,same")' )
        exec( 'legend_%s.Draw()'%(plot) )

raw_input("Please adjust plots. Press ENTER when finish.")

for plot in plots:
  exec( plot+'.cd()')
  exec( plot+'.SaveAs("%s.pdf")'%(plot) )
  exec( plot+'.SaveAs("%s.C")'%(plot) )

file_in.Close()
