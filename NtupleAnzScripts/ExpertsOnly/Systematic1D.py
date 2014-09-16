#!/usr/bin/python
##################################################
#  Calculate the station efficiency uncertainies #
# Author: Jinzhong Zhang(zhangjin@cern.ch),      #
# Northeastern University, Boston, U.S.A         #
##################################################

#Usage: python Systematic1D.py
#All result files must be ready by Step2_PlotAll.py

from  ROOT import *
from Config import *
import sys
if (sys.argv[0] == "python"): args=sys.argv[2:]
else: args=sys.argv[1:]

if Resonance is "Z":
    Data_default="data_Z_Stations/resultplots_data_Z_Stations.root"
    Data_bkgModel="data_Z_Stations/resultplots_data_Z_Stations_BkgModeling.root"
    Data_sigModel="data_Z_Stations/resultplots_data_Z_Stations_SigModeling.root"
    Data_lct_default="data_lct_Z_Stations/resultplots_data_lct_Z_Stations.root"
    Data_lct_bkgModel="data_lct_Z_Stations/resultplots_data_lct_Z_Stations_BkgModeling.root"
    Data_lct_sigModel="data_lct_Z_Stations/resultplots_data_lct_Z_Stations_SigModeling.root"
#    Data_CloseByMu="DATA_Sta_CloseByMu/resultplots.root"
#    Data_LooseProbe="DATA_Sta_LooseProbe/resultplots.root"
    MC_default="mc_Z_Stations/resultplots_mc_Z_Stations.root"
    MC_Counting="mc_Z_Stations/resultplots_mc_Z_Stations_MCTruth.root"
    MC_lct_default="mc_lct_Z_Stations/resultplots_mc_lct_Z_Stations.root"
    MC_lct_Counting="mc_lct_Z_Stations/resultplots_mc_lct_Z_Stations_MCTruth.root"
else:
    Data_default="DATA_Sta_Def/resultplots_JPsi.root"
    Data_bkgModel="DATA_Sta_Def/resultplots_JPsi_BkgModeling.root"
    Data_sigModel="DATA_Sta_Def/resultplots_JPsi_SigModeling.root"
#    Data_CloseByMu="DATA_Sta_CloseByMu/resultplots_JPsi.root"
#    Data_LooseProbe="DATA_Sta_LooseProbe/resultplots_JPsi.root"
    MC_default="MC_Sta_Def/resultplots_JPsi.root"
    MC_Counting="MC_Sta_Def/resultplots_JPsi_Counting.root"    

def DiffBetweenPlots(file1,file2,plotname,sysname="sys"):
    print sysname,":"
    p1=TFile.Open(file1).Get(plotname)
    p2=TFile.Open(file2).Get(plotname)
    nbins=p1.GetN()
    p=TH1D(sysname,sysname,nbins,0,nbins)
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
        p.SetBinContent(x,val)
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
    return p

def Systematic(ListofPlots,Stat,OtherSys=0.,sysname="sys_tot"):
    print "Summary :",
    for iplot in ListofPlots:
        print "\n",iplot.GetTitle(),"&",
        nbins=iplot.GetNbinsX()*iplot.GetNbinsY()
        if iplot.GetNbinsY()==1:
            for x in range(nbins):
                print "%.3f" % iplot.GetBinContent(x),"&",
    print "\n",sysname,"&",
    p=ListofPlots[0].Clone(sysname)
    for x in range(nbins):
        sumover=0
        for iplot in ListofPlots:
            sumover+=iplot.GetBinContent(x)*iplot.GetBinContent(x)
        sumover+=OtherSys*OtherSys
        if p.GetNbinsY()==1:
            print "%.3f" % sqrt(sumover),"&",
        p.SetBinContent( x,sqrt(sumover) )
    print "\n=============Final Result","(adding other systematic uncertainty,",str(OtherSys)+"%)","=============="
    val=Stat.GetY()
    errhi=Stat.GetEYhigh()
    errlo=Stat.GetEYlow()
    print "stat high error&",
    for x in range(nbins):
        print "%.3f" % errhi[x],"&",
    print "\nstat low error&",
    for x in range(nbins):
        print "%.3f" % errlo[x],"&",
    print "\nStat=max(hi,lo) &",
    staterr=[0]*nbins
    for x in range(nbins):
        staterr[x]=max(errhi[x],errlo[x])
        print "%.3f" % staterr[x],"&",
    print "\nValue&",
    for x in range(nbins):
        print "%.3f" % val[x],"&",
    print "\nStat+Sys &",
    errtot=[0]*nbins
    for x in range(nbins):
        syserr=p.GetBinContent(x)
        errtot[x]=sqrt(staterr[x]*staterr[x]+syserr*syserr)
        print "%.3f" % errtot[x],"&",
    print "\nFinal Result (val$\pm$Sys$\pm$Stat) &",
    for x in range(nbins):
        print "$%.2f\pm%.2f\pm%.2f$" % (val[x],p.GetBinContent(x),staterr[x]), "&",
    print "\nFinal Result &",
    for x in range(nbins):
        print "$%.2f\pm%.2f$" % (val[x],errtot[x]), "&",
    print "\n\n"
    return p

#DiffBetweenPlots("ZMuMuFit_OldSelect/resultplots.root","ZMuMuFit_OldSelect/resultplots_Counting.root","SEGEff_St","SEGMCTruth")

SEGSysList=[ DiffBetweenPlots(MC_default,MC_Counting,"SEGEff","SEGMCTruth"),
#             DiffBetweenPlots(Data_default,MC_default,"SEGEff_St","SEGDATAMC"),
             DiffBetweenPlots(Data_default,Data_bkgModel,"SEGEff","SEG_BKGModeling"),
             DiffBetweenPlots(Data_default,Data_sigModel,"SEGEff","SEG_SIGModeling"),
#             DiffBetweenPlots(Data_default,Data_CloseByMu,"SEGEff_St","SEG_CloseByMu"),
#             DiffBetweenPlots(Data_default,Data_LooseProbe,"SEGEff_St","SEG_LooseProbe")
             ]

LCTSysList=[ DiffBetweenPlots(MC_lct_default,MC_lct_Counting,"LCTEff","LCTMCTruth"),
#             DiffBetweenPlots(Data_default,MC_default,"LCTEff_St","LCTDATAMC"),
             DiffBetweenPlots(Data_lct_default,Data_lct_bkgModel,"LCTEff","LCT_BKGModeling"),
             DiffBetweenPlots(Data_lct_default,Data_lct_sigModel,"LCTEff","LCT_SIGModeling"),
#             DiffBetweenPlots(Data_default,Data_CloseByMu,"LCTEff_St","LCT_CloseByMu"),
#             DiffBetweenPlots(Data_default,Data_LooseProbe,"LCTEff_St","LCT_LooseProbe")
             ]

print
Systematic(SEGSysList, TFile.Open(Data_default).Get("SEGEff"),2, "SEG Total" )
Systematic(LCTSysList, TFile.Open(Data_lct_default).Get("LCTEff"),2, "LCT Total")
