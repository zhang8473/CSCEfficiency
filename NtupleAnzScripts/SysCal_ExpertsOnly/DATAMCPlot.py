#!/usr/bin/python
##################################################
#       Make the data MC comparison plots        #
# Author: Jinzhong Zhang(zhangjin@cern.ch),      #
# Northeastern University, Boston, U.S.A         #
##################################################

#Usage: python DATAMCPlot.py

from  ROOT import *
from Config import *
import sys
if (sys.argv[0] == "python"): args=sys.argv[2:]
else: args=sys.argv[1:]

gROOT.SetStyle("Plain")
gStyle.SetOptStat(0)
gStyle.SetFillColor(10)
gStyle.SetOptTitle(0)

if AnalyzeZPeak:
    Data_default="DATA_Sta_Def/resultplots.root"
    MC_default="MC_Sta_Def/resultplots.root"
else:
    Data_default="DATA_Sta_Def/resultplots_JPsi.root"
    MC_default="MC_Sta_Def/resultplots_JPsi.root"

def SetPlotStyle(p,color=kBlack):
    p.SetLineColor(color)
    p.SetMarkerColor(color)
    p.SetMarkerSize(.5)
    p.SetMarkerStyle(8)

def DiffBetweenPlots(file1,file2,plotname,sysname="sys",name1="DATA",name2="MC"):
    print sysname,":"
    p1=TFile.Open(file1).Get(plotname)
    p2=TFile.Open(file2).Get(plotname)
    nbins=p1.GetN()
    val1=p1.GetY()
    errhi1=p1.GetEYhigh()
    errlo1=p1.GetEYlow()
    val2=p2.GetY()
    errhi2=p2.GetEYhigh()
    errlo2=p2.GetEYlow()
    print "p1","&",
    for x in range(nbins):
        print "%.2f" % val1[x],"&",
    print "\np2","&",
    for x in range(nbins):
        print "%.2f" % val2[x],"&",
    print "\nsys","&",
    for x in range(nbins):
        val=abs(val1[x]-val2[x])
        print "%.3f" % val,"&",
    #Here the way to calculate errors is not correct, but just approximation
    print "\nsys_err_high","&",
    for x in range(nbins):
        errhi=max(errhi1[x],errhi2[x])
        print "%.2f" % errhi,"&",
    print "\nsys_err_low","&",
    for x in range(nbins):
        errlo=max(errlo1[x],errlo2[x])
        print "%.2f" % errlo,"&",
    print
    canvas=TCanvas(sysname,sysname,500,500)
    legend=TLegend(0.6,0.2,0.8,0.4)
    p1.Draw("APE")
    p2.Draw("PE,same")
    legend.AddEntry(p1,name1,"LPE")
    legend.AddEntry(p2,name2,"LPE")
    legend.Draw()
    SetPlotStyle(p1,kBlack)
    SetPlotStyle(p2,kGreen)
    raw_input("Press ENTER to continue")
    return canvas

#DiffBetweenPlots("ZMuMuFit_OldSelect/resultplots.root","ZMuMuFit_OldSelect/resultplots_Counting.root","SEGEff_St","SEGMCTruth")
if AnalyzeZPeak:
    DiffBetweenPlots(Data_default,MC_default,"SEGEff_St","DATAMC_SEG","DATA","Z->#mu#mu MC")
    DiffBetweenPlots(Data_default,MC_default,"LCTEff_St","DATAMC_LCT","DATA","Z->#mu#mu MC")
else:
    DiffBetweenPlots(Data_default,MC_default,"SEGEff_St","DATAMC_SEG","DATA","J/#Psi->#mu#mu MC")
    DiffBetweenPlots(Data_default,MC_default,"LCTEff_St","DATAMC_LCT","DATA","J/#Psi->#mu#mu MC")    
