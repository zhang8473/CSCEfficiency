#!/usr/bin/python
##################################################
# Calculate the chamber efficiency uncertainies  #
# Author: Jinzhong Zhang(zhangjin@cern.ch),      #
# Northeastern University, Boston, U.S.A         #
##################################################

#Usage: python Systematic2D.py
#All result files must be ready by Step2_PlotAll.py

from  ROOT import *
import sys
if (sys.argv[0] == "python"): args=sys.argv[2:]
else: args=sys.argv[1:]

Data_default="DATA_Chamber_Def/resultplots.root"
Data_bkgModel="DATA_Chamber_Def/resultplots_BkgModeling.root"
Data_sigModel="DATA_Chamber_Def/resultplots_SigModeling.root"
Data_CloseByMu="DATA_Chamber_CloseMu/resultplots.root"
Data_LooseProbe="DATA_Chamber_LooseProbe/resultplots.root"
MC_default="MC_Chamber_Def/resultplots.root"
MC_Counting="MC_Chamber_Def/resultplots_Counting.root"

gROOT.SetStyle("Plain")
gStyle.SetPaintTextFormat("0.3g")
gStyle.SetOptStat(0)

def DiffBetweenPlots(file1,file2,plotname,sysname="sys"):
    p1=TFile.Open(file1).Get(plotname)
    p2=TFile.Open(file2).Get(plotname)
    p=p1-p2
    p.SetName(sysname)
    for x in range( p.GetNbinsX() ):
        for y in range( p.GetNbinsY() ):
            if p1.GetBinContent(x,y)>0 or p2.GetBinContent(x,y)>0:
                val=abs(p.GetBinContent(x,y))
                if val<1E-3:
                    val=1E-5
                p.SetBinContent( x,y,val )
    p.SetMaximum(100)
    p.SetMinimum(0)
    exec( sysname+"=TCanvas(\"%s\",\"%s\",500,500)"%(sysname,sysname) )
    exec( sysname+".cd()")
    p.Draw("COLZ,TEXT")
    raw_input(sysname+": Press ENTER to continue")
    return p

def Systematic(ListofPlots,sysname="sys_tot"):
    for iplot in ListofPlots:
        nbins=iplot.GetNbinsX()*iplot.GetNbinsY()
    p=ListofPlots[0].Clone(sysname)
    p.SetTitle(sysname)
    for x in range( 1,nbins+1 ):
        sumover=0
        for iplot in ListofPlots:
            sumover+=iplot.GetBinContent(x)*iplot.GetBinContent(x)
        p.SetBinContent( x,sqrt(sumover) )
    return p

#DiffBetweenPlots(MC_default,Data_default,"SEGEff","SEGDATAMC")
#DiffBetweenPlots(MC_default,Data_default,"LCTEff","LCTDATAMC")

SEGSysList=[ DiffBetweenPlots(MC_default,MC_Counting,"SEGEff","SEGMCTruth"),
#             DiffBetweenPlots(Data_default,MC_default,"SEGEff","SEGDATAMC"),
             DiffBetweenPlots(Data_default,Data_bkgModel,"SEGEff","SEG_BKGModeling"),
             DiffBetweenPlots(Data_default,Data_sigModel,"SEGEff","SEG_SIGModeling"),
             DiffBetweenPlots(Data_default,Data_CloseByMu,"SEGEff","SEG_CloseByMu"),
             DiffBetweenPlots(Data_default,Data_LooseProbe,"SEGEff","SEG_LooseProbe")
             ]

LCTSysList=[ DiffBetweenPlots(MC_default,MC_Counting,"LCTEff","LCTMCTruth"),
             # DiffBetweenPlots(Data_default,MC_default,"LCTEff_St","LCTDATAMC"),
             DiffBetweenPlots(Data_default,Data_bkgModel,"LCTEff","LCT_BKGModeling"),
             DiffBetweenPlots(Data_default,Data_sigModel,"LCTEff","LCT_SIGModeling"),
             DiffBetweenPlots(Data_default,Data_CloseByMu,"LCTEff","LCT_CloseByMu"),
             DiffBetweenPlots(Data_default,Data_LooseProbe,"LCTEff","LCT_LooseProbe")
             ]

SegCanvas=TCanvas("segment efficiency systematic","segment efficiency systematic",500,500)
LCTCanvas=TCanvas("LCT efficiency systematic","LCT efficiency systematic",500,500)

SegSysTotal=Systematic(SEGSysList,"SEG Total (%)")
SegCanvas.cd()
#SEGSysList[0].Draw("COLZ,TEXT")
SegSysTotal.Draw("COLZ,TEXT")

LCTSysTotal=Systematic(LCTSysList,"LCT Total (%)")
LCTCanvas.cd()
#LCTSysList[0].Draw("COLZ,TEXT")
LCTSysTotal.Draw("COLZ,TEXT")
raw_input("Press ENTER to exit")
