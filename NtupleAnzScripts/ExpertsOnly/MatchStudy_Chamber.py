#!/usr/bin/python
#######################################################
# Tracker track to LCT/Segments distance study script #
# Author: Zhang Jinzhong, zhangjin@cern.ch            #
#######################################################
#Usage: python MatchStudy.py TheNtupleRootFile.root
from  ROOT import *
from Config import Group,RunOnMC,chambers,DenominatorRequire
import sys,os

gROOT.SetStyle("Plain")
gStyle.SetOptStat(0)

if sys.argv[0] == "python":
    args=sys.argv[2:]
else:
    args=sys.argv[1:]

def Getch():
  import tty, termios
  fd = sys.stdin.fileno()
  old_settings = termios.tcgetattr(fd)
  try:
    tty.setraw(sys.stdin.fileno())
    ch = sys.stdin.read(1)
  finally:
    termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
  return ch

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

#content order: name:(formula,x_min,x_max,y_max)
plots={
       "TracktoSegX":('CSCDxTTSeg#',-40.,40.),
       "TracktoSegY":('CSCDyTTSeg#',-40.,40.),
#       "TracktoSegXY":('CSCDxyTTSeg#',-1.,40.),
       "TracktoLCTX":('CSCDxTTLCT#',-40,40.),
       "TracktoLCTY":('CSCDyTTLCT#',-40.,40.),
#       "TracktoLCTXY":('CSCDxyTTLCT#',0.,40.),
#       "TracktoSegXY_Check":('sqrt(pow(CSCDxTTSeg#,2)+pow(CSCDyTTSeg#,2))',-1.,150.),# should be the same with TracktoSegR
#       "TracktoLCTXY_Check":('sqrt(CSCDxTTLCT#*CSCDxTTLCT#+CSCDyTTLCT#*CSCDyTTLCT#)',-1.,150.)# should be the same with TracktoLCTR
#       "SegtoLCTX":('CSCSegxLc#-CSCLCTxLc#',-40.,40.),
#       "SegtoLCTY":('CSCSegyLc#-CSCLCTyLc#',-40.,40.),
#       "SegtoLCTXY":('sqrt(pow(CSCSegxLc#-CSCLCTxLc#,2)+pow(CSCSegyLc#-CSCLCTyLc#,2))',-1.,40.)
#        "tracks_phi":('tracks_phi',-7.,7.)
       }

chambers_to_draw=["ME-1_1_29","ME-1_1_30","ME-1_1_31","ME-1_1_32"]
colors_to_draw=[kMagenta,kBlack,kRed-9,kRed,kGreen,kBlue]
ncolors_to_draw=len(colors_to_draw)

def MakeAPlot(Name_,Forumla_,Cut_,color_):
    if RunOnMC:
        tree_in.Project(Name_,Forumla_,'weight*(%s)'%(Cut_))
    else:
        tree_in.Project(Name_,Forumla_,Cut_)
    exec( 'NEntries='+Name_+'.GetEntries()')
    if NEntries>0:
        exec( Name_+'.Scale( 100/'+Name_+'.Integral() )' )
        exec( Name_+'.SetLineColor(color_)' )
        exec( Name_+'.SetMarkerColor(color_)' )
        exec( 'result=(%s.GetMean(),%s.GetRMS())'%(Name_,Name_) )
        print Name_,result
        return result
    else:
        print "Skip "+Name_+", no entry."
        return (0,0)
    
for plot in plots:
    exec( plot+'=TCanvas("%s","%s",500,500)'%(plot,plot) )
    exec( plot+"_Max=0" )
    exec( 'legend_%s=TLegend(0.6,0.6,0.8,0.8)'%(plot) )
    exec( 'legend_%s.SetFillColor(0)'%(plot) )

icurve=0
for idx in chambers:
    Skip=True
    for ichamber in chambers_to_draw:
        if ichamber==idx:
            Skip=False
            break
    if Skip:
        continue
    ec=idx[2] is '+'
    st=int(idx[3])
    rg=int(idx[5])
    ch=int(idx[7:])
    Cut="%sCSCEndCapPlus && CSCRg%s==%s && CSCCh%s==%s &&"%("" if (ec) else '!',st,rg,st,ch)+DenominatorRequire.replace("#",idx[3])#+" && "+InvariantMass
    for plot in plots:
        Name='%s_%s'%( plot,idx.replace("+","plus").replace("-","minus") )
        exec( Name+'=TH1D(Name,plot+";cm;scaled number of entries",100,%f,%f)'%(plots[plot][1],plots[plot][2]) )
        MeanAndRMS=MakeAPlot(Name,plots[plot][0].replace('#',idx[3]),Cut,colors_to_draw[icurve%ncolors_to_draw])
        exec( "max="+Name+".GetBinContent("+Name+".GetMaximumBin())" )
        exec( "if "+plot+"_Max<max: "+plot+"_Max=max" )
        exec( 'legend_%s.AddEntry('%(plot)+Name+',\"'+idx+': mean:%.2fcm;RMS:%.2fcm\",\"l\")'%(MeanAndRMS[0],MeanAndRMS[1]) )
    icurve += 1

#raw_input("Please adjust plots. Press ENTER when finish.") # exit here if you only want to print out mean and RMS

for plot in plots:
    exec( plot+'.cd()')
    DrawOption="H"
    for idx in chambers:
        Name='%s_%s'%( plot,idx.replace("+","plus").replace("-","minus") )
        if Name not in globals():
            continue
        exec( Name+'.Draw(\"'+DrawOption+'\")' )
        if DrawOption is "H":
            exec( Name+'.SetMaximum(%s_Max*1.05)'%(plot) )
            DrawOption="H,same"
        exec( 'legend_%s.Draw()'%(plot) )
        exec( plot+'.SetLogy()')

print "Save plots(y/n)?",
ch=Getch()
if ch=="y":
    for plot in plots:
        exec( plot+'.cd()')
        exec( plot+'.SaveAs("%s.pdf")'%(plot) )
        exec( plot+'.SaveAs("%s.C")'%(plot) )
file_in.Close()
if RunOnMC:
    tmpoutputfile.Close()
    import os
    os.system( "rm "+TemporaryOutputFile )
