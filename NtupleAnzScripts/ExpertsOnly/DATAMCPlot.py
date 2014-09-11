#!/usr/bin/python
##################################################
#       Make the data MC comparison plots        #
# Author: Jinzhong Zhang(zhangjin@cern.ch),      #
# Northeastern University, Boston, U.S.A         #
##################################################

#Usage: python DATAMCPlot.py datafile mcfile plotname

#AUXPtDATAFile="" #no JPsi file for low pt
#AUXPtDATAFile="data_lct_JPsi_Stationspt/resultplots_data_lct_JPsi_Stationspt.root" # for lct
#AUXPtDATAFile="data_JPsi_Stationspt/resultplots_data_JPsi_Stationspt.root" # for seg
#AUXPtMCFile=""  #no JPsi file for low pt
#AUXPtMCFile="mc_lct_JPsi_Stationspt/resultplots_mc_lct_JPsi_Stationspt.root" # for lct
#AUXPtMCFile="mc_JPsi_Stationspt/resultplots_mc_JPsi_Stationspt.root" # for seg
plotmin=0.88
plotmax=1.0


from  ROOT import *
import sys,os

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

if (sys.argv[0] == "python"): args=sys.argv[2:]
else: args=sys.argv[1:]

gROOT.SetStyle("Plain")
gStyle.SetOptStat(0)
gStyle.SetFillColor(10)
gStyle.SetTitleAlign(33)
gStyle.SetTitleX(.99)
gStyle.SetPadLeftMargin(0.15)
#gStyle.SetOptTitle(0)

if len(args)<3:
  print "#Usage: python DATAMCPlot.py datafile mcfile plotname"
  sys.exit(0)

#-----Get the data file name and path, import the Config file from the data path if possible
if len(args)>0:
  DATA=args[0]
  path="/".join(DATA.split("/")[0:-1])

if os.path.isfile(path+"/"+"__init__.py"):
    exec( "from %s.Config import *"%(path) )
else:
    from Config import *

#-----Get the MC file name and path, draw MCTruth if the file name includes the keyword "MCTruth"
if len(args)>1:
    MC=args[1]
if "MCTruth" in MC:
    CompareMCTruth=True
    if "AUXPtMCFile" in globals():
      print "change", AUXPtMCFile, "to",
      AUXPtMCFile=AUXPtMCFile[0:-5]+"_MCTruth.root"
      print AUXPtMCFile
else:
    CompareMCTruth=False

#-----Get the plotname, it should be the same in both data and mc files
if len(args)>2:
    objname=args[2]
else:
    objname="ME4lct_effV"

def SetPlotStyle(p,color=kBlack):
    p.SetLineColor(color)
    p.SetLineWidth(2)
    p.SetMarkerColor(color)
    p.SetFillColor(color)
    p.SetMarkerSize(.5)
    p.SetMarkerStyle(8)
    p.SetFillStyle(3001)
    if "pt" in DATA:
        p.GetXaxis().SetLimits(6., 100.);
    if "eta" in DATA:
        p.GetXaxis().SetLimits(0.9, 2.4);
    if "phi" in DATA:
        p.GetXaxis().SetLimits(-0.0872664626,6.1959188464);
    p.SetMinimum(plotmin)
    p.SetMaximum(plotmax)

#Method 1: for only the first few low pt bins until the Z plot is better
def CombinedPtPlots(file1,file2,plotname,name1="DATA",name2="MC"):
    pdata=TFile.Open(file1).Get(plotname)
    pmc=TFile.Open(file2).Get(plotname+"_MCTruth" if CompareMCTruth else plotname)
    pdata_JPsi=TFile.Open(AUXPtDATAFile).Get(plotname)
    pmc_JPsi=TFile.Open(AUXPtMCFile).Get(plotname+"_MCTruth" if CompareMCTruth else plotname)
    nbins=min(pdata.GetN(),pdata_JPsi.GetN())
    x=pdata.GetX()
    xtitle=pdata.GetXaxis().GetTitle()
    ytitle=pdata.GetYaxis().GetTitle()
    errhi1=pdata.GetEYhigh()
    errlo1=pdata.GetEYlow()
    errhi1_JPsi=pdata_JPsi.GetEYhigh()
    errlo1_JPsi=pdata_JPsi.GetEYlow()
    for ibin in range(nbins):
      if pdata.GetX()[ibin]!=pdata_JPsi.GetX()[ibin] or pdata.GetErrorXlow(ibin)!=pdata_JPsi.GetErrorXlow(ibin):
        print "Warning: bin",ibin,":Z pt:",pdata.GetX()[ibin],"; JPsi pt:",pdata_JPsi.GetX()[ibin]
        continue
      print "pt=",pdata.GetX()[ibin],"eff(Z)=",pdata.GetY()[ibin],"+/-",max(errhi1[ibin],errlo1[ibin]),"(MC)=",pmc.GetY()[ibin],"; eff(JPsi)=",pdata_JPsi.GetY()[ibin],"+/-",max(errhi1[ibin],errlo1_JPsi[ibin]),"(MC)=",pmc_JPsi.GetY()[ibin]
      #return canvas
      #if pdata.GetX()[ibin]<20:
      if pdata.GetX()[ibin]<10 or max(errhi1[ibin],errlo1[ibin])>max(errhi1_JPsi[ibin],errlo1_JPsi[ibin]):
        pdata.SetPoint(ibin,x[ibin],pdata_JPsi.GetY()[ibin])
        pdata.SetPointError(ibin,pdata_JPsi.GetErrorXlow(ibin),pdata_JPsi.GetErrorXhigh(ibin),errlo1_JPsi[ibin],errhi1_JPsi[ibin])
        pmc.SetPoint(ibin,x[ibin],pmc_JPsi.GetY()[ibin])
      else:
        ibin-=1
        break
        
    splitbin=ibin
    canvas=TCanvas(plotname,plotname,500,500)
    legend=TLegend(0.6,0.2,0.8,0.4)
    SetPlotStyle(pdata,kBlack)
    SetPlotStyle(pmc,kRed)
    #root BUG, after SetPoint, title is missing
    pdata.GetXaxis().SetTitle(xtitle)
    pdata.GetYaxis().SetTitle(ytitle)
    pdata.GetYaxis().SetTitleOffset(1.53)
    pdata.Draw("A2")
    pmc.Draw("E,same" if CompareMCTruth else "E,same")
    legend.AddEntry(pdata,name1,"F")
    legend.AddEntry(pmc,name2,"L")
    legend.Draw()
    if splitbin>=0:
      line = TLine(x[splitbin]+pdata.GetErrorXhigh(splitbin),plotmin,x[splitbin]+pdata.GetErrorXhigh(splitbin),plotmax);
      line.SetLineColor(kBlue)
      line.SetLineWidth(2)
      line.SetLineStyle(2)
      line.Draw()
      text=TLatex()
      text.SetTextColor(kBlue)
      text.DrawLatex(x[splitbin]+pdata.GetErrorXhigh(splitbin)-15,plotmax+0.002,"J/#Psi#leftarrow #rightarrow Z")
    print "Save it(y/n)?",
    ch=Getch()
    if ch=="y":
        canvas.SaveAs(objname+".pdf")
    return canvas

#Method2: For pt plots, if difference between data and mc is larger than the one of JPsi or errorbar is larger then the one of JPsi then take the JPsi bin, For all bins
def CombinedPtPlots2(file1,file2,plotname,name1="DATA",name2="MC"):
    pdata=TFile.Open(file1).Get(plotname)
    pmc=TFile.Open(file2).Get(plotname+"_MCTruth" if CompareMCTruth else plotname)
    pdata_JPsi=TFile.Open(AUXPtDATAFile).Get(plotname)
    pmc_JPsi=TFile.Open(AUXPtMCFile).Get(plotname+"_MCTruth" if CompareMCTruth else plotname)
    Nbins=min(pdata.GetN(),pdata_JPsi.GetN())
    x=pdata.GetX()
    xtitle=pdata.GetXaxis().GetTitle()
    ytitle=pdata.GetYaxis().GetTitle()
    errhi1=pdata.GetEYhigh()
    errlo1=pdata.GetEYlow()
    errhi1_JPsi=pdata_JPsi.GetEYhigh()
    errlo1_JPsi=pdata_JPsi.GetEYlow()
    
    N_removedpoints_Z=0
    N_removedpoints_JPsi=0
    for ibin in range(Nbins):
      ibin_Z=ibin-N_removedpoints_Z
      ibin_JPsi=ibin-N_removedpoints_JPsi
      print "eff(Z)=",pdata.GetY()[ibin_Z],"+/-",max(errhi1[ibin_Z],errlo1[ibin_Z]),"(MC)=",pmc.GetY()[ibin_Z],"; eff(JPsi)=",pdata_JPsi.GetY()[ibin_JPsi],"+/-",max(errhi1_JPsi[ibin_JPsi],errlo1_JPsi[ibin_JPsi]),"(MC)=",pmc_JPsi.GetY()[ibin_Z]
      if pdata.GetX()[ibin_Z]!=pdata_JPsi.GetX()[ibin_JPsi] or pdata.GetErrorXlow(ibin_Z)!=pdata_JPsi.GetErrorXlow(ibin_JPsi):
        print "Warning: bin",ibin,":Z pt:",pdata.GetX()[ibin_Z],"; JPsi pt:",pdata_JPsi.GetX()[ibin_JPsi]
        continue
      if abs(pdata.GetY()[ibin_Z]-pmc.GetY()[ibin_Z])>abs(pdata_JPsi.GetY()[ibin_JPsi]-pmc_JPsi.GetY()[ibin_JPsi]) or max(errhi1[ibin_Z],errlo1[ibin_Z])>max(errhi1_JPsi[ibin_JPsi],errlo1_JPsi[ibin_JPsi]):
        pdata.RemovePoint(ibin_Z)
        pmc.RemovePoint(ibin_Z)
        N_removedpoints_Z+=1
      else:
        pdata_JPsi.RemovePoint(ibin_JPsi)
        pmc_JPsi.RemovePoint(ibin_JPsi)
        N_removedpoints_JPsi+=1
        
    canvas=TCanvas(plotname,plotname,500,500)
    legend=TLegend(0.6,0.2,0.8,0.4)
    SetPlotStyle(pdata,kBlack)
    SetPlotStyle(pdata_JPsi,kBlue)
    SetPlotStyle(pmc,kRed)
    SetPlotStyle(pmc_JPsi,kRed)
    #root BUG, after SetPoint, title is missing
    pdata.GetXaxis().SetTitle(xtitle)
    pdata.GetYaxis().SetTitle(ytitle)
    pdata.GetYaxis().SetTitleOffset(1.53)
    pdata.Draw("A2")
    pdata_JPsi.Draw("2,same")
    pmc.Draw("E,same" if CompareMCTruth else "E,same")
    pmc_JPsi.Draw("E,same" if CompareMCTruth else "E,same")
    legend.AddEntry(pdata,name1+"(Z)","F")
    legend.AddEntry(pdata_JPsi,name1+"(J/#Psi)","F")
    legend.AddEntry(pmc,name2,"L")
    legend.Draw()
    print "Save it(y/n)?",
    ch=Getch()
    if ch=="y":
        canvas.SaveAs(objname+".pdf")
    return canvas

def DiffBetweenPlots(file1,file2,plotname,name1="DATA",name2="MC"):
    pdata=TFile.Open(file1).Get(plotname)
    pmc=TFile.Open(file2).Get(plotname+"_MCTruth" if CompareMCTruth else plotname)
    if pdata.__class__.__name__!="TGraphAsymmErrors":
      print "\033[91m Error:",plotname,"does not exist in",file1,"\033[0m"
      return
    if pmc.__class__.__name__!="TGraphAsymmErrors":
      print "\033[91m Error:",plotname+"_MCTruth" if CompareMCTruth else plotname,"does not exist in",file1,"\033[0m"
      return
    nbins=pdata.GetN()
    x=pdata.GetX()
    val1=pdata.GetY()
    errhi1=pdata.GetEYhigh()
    errlo1=pdata.GetEYlow()
    val2=pmc.GetY()
    errhi2=pmc.GetEYhigh()
    errlo2=pmc.GetEYlow()
    print " x","&",
    for ibin in range(nbins):
        print "%.2f" % x[ibin],"&",
    print "\n pdata","&",
    for ibin in range(nbins):
        print "%.2f" % val1[ibin],"&",
    print "\n pmc","&",
    for ibin in range(nbins):
        print "%.2f" % val2[ibin],"&",
    print "\ndif","&",
    for ibin in range(nbins):
        val=abs(val1[ibin]-val2[ibin])
        print "%.2f" % val,"&",
    #Here the way to calculate errors is not correct, but just approximation
    print "\nmax_err_high","&",
    for ibin in range(nbins):
        errhi=max(errhi1[ibin],errhi2[ibin])
        print "%.2f" % errhi,"&",
    print "\nmax_err_low ","&",
    for ibin in range(nbins):
        errlo=max(errlo1[ibin],errlo2[ibin])
        print "%.2f" % errlo,"&",
    print
    canvas=TCanvas(plotname,plotname,500,500)
    legend=TLegend(0.6,0.2,0.8,0.4)
    pdata.GetYaxis().SetTitleOffset(1.53)
    pdata.Draw("A2")
    if CompareMCTruth:
        pmc.Draw("E,same")
    else:
        pmc.Draw("E,same")
    legend.AddEntry(pdata,name1,"F")
    legend.AddEntry(pmc,name2,"L")
    legend.Draw()
    SetPlotStyle(pdata,kBlack)
    SetPlotStyle(pmc,kRed)
    print "Save it(y/n)?",
    ch=Getch()
    if ch=="y":
        canvas.SaveAs(objname+".pdf")
    return canvas

#DiffBetweenPlots("ZMuMuFit_OldSelect/resultplots.root","ZMuMuFit_OldSelect/resultplots_Counting.root","SEGEff_St","SEGMCTruth")
if "pt" in Group and os.path.isfile(AUXPtDATAFile) and os.path.isfile(AUXPtMCFile):
  print "combine the results from \n DATA: ",DATA,AUXPtDATAFile,"\n MC: ",MC,AUXPtMCFile
  CombinedPtPlots(DATA,MC,objname,"DATA","MC"+("Truth" if CompareMCTruth else "") )
elif Resonance=="Z":
  DiffBetweenPlots(DATA,MC,objname,"DATA","Z->#mu#mu MC"+("Truth" if CompareMCTruth else "") )
else:
  DiffBetweenPlots(DATA,MC,objname,"DATA","J/#Psi->#mu#mu MC"+("Truth" if CompareMCTruth else "") )  
