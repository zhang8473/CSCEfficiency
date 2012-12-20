#!/usr/bin/python
#Author: Jinzhong Zhang(zhangjin@cern.ch), Northeastern University, Boston, U.S.A
#################################################################################
#Advanced Usage:
#python Step2_PlotAll.py arg1 arg2
#   arg1 is the "TnP_" files directory path;
#   arg2 is the postfix of root path name "lct_effV" and "seg_effV";
#   arg2 can be specified as "bkg" or "sig" for background and signal modeling
#Example1(plot default efficiencies): python Step2_PlotAll.py
#Example2(systematic -- bkg modeling): python Step2_PlotAll.py /uscms/home/zhangjin/nobackup/ bkg
#Example3(systematic -- sig modeling): python Step2_PlotAll.py /uscms/home/zhangjin/nobackup/ sig
#Example2 and 3 are used for systematic calculation.
##################################################################################

from  ROOT import *
from  numpy import *
from Config import *
gROOT.SetStyle("Plain")
gStyle.SetPaintTextFormat("0.3g")
gStyle.SetOptStat(0)

from array import array as arr
Red = arr('d',[0.00, 0.00, 0.00, 1.00, 1.00])
Green = arr('d',[0.00, 0.10, 1.00, 1.00, 0.00])
Blue = arr('d',[1.00, 0.90, 0.00, 0.00, 0.00])
Length = arr('d',[0.00, 0.85, 0.90, 0.95, 1.00])
TColor.CreateGradientColorTable(5,Length,Red,Green,Blue,500)
gStyle.SetNumberContours(500)

import sys,os
if (sys.argv[0] == "python"): args=sys.argv[2:]
else: args=sys.argv[1:]

Prefix=""
Postfix=""
if len(args)>0:
    Prefix=args[0]
    if Prefix[-1] != "/":
        Prefix+="/"
    TagProbeFitResult=TagProbeFitResult.split("/")[-1]
    ResultPlotsFileName=ResultPlotsFileName.split("/")[-1]
    if len(args)>1:
        if args[1] == "bkg":
            Postfix="_BkgModeling"
        elif args[1] == "sig":
            Postfix="_SigModeling"
        else:
            Postfix=args[1]
    ResultPlotsFileName=Prefix+ResultPlotsFileName.replace(".root",Postfix+".root")

file_out=TFile.Open(ResultPlotsFileName,'RECREATE')

def GetEff(filename_,lctpath="lct_effV",segpath="seg_effV"):
    import os
    try:
        if os.path.isfile(filename_):
            f_in =TFile(filename_,"READ");
            LCTEff=f_in.Get("Probes/"+lctpath+"/fit_eff").get().find("efficiency")
            Segeff=f_in.Get("Probes/"+segpath+"/fit_eff").get().find("efficiency")
            return [Segeff.getVal(),abs(Segeff.getErrorLo()),Segeff.getErrorHi(),LCTEff.getVal(),abs(LCTEff.getErrorLo()),LCTEff.getErrorHi()]
        else:
            print filename_+" is not found, skip.. "
            return [0.]*6
    except:
        print "\033[97mOther problems, skip",filename_,"\033[0m"
        #os.remove(filename_)
        return [0.]*6 
    
if AnalyzeStations:
    fitEffs=[]
    for idx in range(1,n_stations+1):
      fitEffs.append( GetEff( "%s%s%s.root" % (Prefix,TagProbeFitResult,stations[idx][1]),"lct_effV"+Postfix,"seg_effV"+Postfix ) )
    fitEffs=array(fitEffs).transpose()*100.
    xval=array(range(1,n_stations+1))*1.0
    xerr=zeros(n_stations, dtype=float)
    SegEff=TGraphAsymmErrors(n_stations, xval, array(fitEffs[0]), xerr, xerr, array(fitEffs[1]), array(fitEffs[2]))
    LCTEff=TGraphAsymmErrors(n_stations, xval, array(fitEffs[3]), xerr, xerr, array(fitEffs[4]), array(fitEffs[5]))
    SegCanvas=TCanvas("segment efficiency","segment efficiency",500,500)
    SegCanvas.cd()
    SegEff.SetMaximum(100)
    SegEff.SetMinimum(90)
    LCTEff.SetMaximum(100)
    LCTEff.SetMinimum(90)
    SegEff.Draw("AP")
    LCTCanvas=TCanvas("lct efficiency","lct efficiency",500,500)
    LCTCanvas.cd()
    LCTEff.Draw("AP")
    for st in range(1,n_stations+1):
       binnum=SegEff.GetXaxis().FindBin(st)
       SegEff.GetXaxis().SetBinLabel( binnum,stations[st][1] )
       LCTEff.GetXaxis().SetBinLabel( binnum,stations[st][1] )
    file_out.cd()
    SegEff.Write("SEGEff")
    LCTEff.Write("LCTEff")
else:
    SEGEff=TH2F("SEGEff","segment efficiency",36,1,37,18,-9,9)
    SEGEff.SetMarkerSize(0.7)
    SEGEff.SetContour(500)
    SEGEff_upErr=TH2F("SEGEff_upErr","segment efficiency uperror",36,1,37,18,-8.7,9.3)
    SEGEff_upErr.SetMarkerSize(0.45)
    SEGEff_downErr=TH2F("SEGEff_downErr","segment efficiency loerror",36,1,37,18,-9.3,8.7)
    SEGEff_downErr.SetMarkerSize(0.45)
    Chambers_  = ["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36"]
    Rings_ = ["ME-42","ME-41","ME-32","ME-31","ME-22","ME-21","ME-13","ME-12","ME-11","ME+11","ME+12","ME+13","ME+21","ME+22","ME+31","ME+32","ME+41","ME+42"]
    for ich in range(36):
        SEGEff.GetXaxis().SetBinLabel(ich+1,Chambers_[ich])
        SEGEff_upErr.GetXaxis().SetBinLabel(ich+1,Chambers_[ich])
        SEGEff.GetXaxis().SetBinLabel(ich+1,Chambers_[ich])
    for irg in range(18):
        SEGEff.GetYaxis().SetBinLabel(irg+1,Rings_[irg])
        SEGEff_upErr.GetYaxis().SetBinLabel(irg+1,Rings_[irg])
        SEGEff_downErr.GetYaxis().SetBinLabel(irg+1,Rings_[irg])
    LCTEff=SEGEff.Clone("LCTEff")
    LCTEff.SetTitle("LCT efficiency")
    LCTEff_upErr=SEGEff_upErr.Clone("LCTEff_upErr")
    LCTEff_downErr=SEGEff_downErr.Clone("LCTEff_downErr")
    RingToYMap={(1,1):0,(1,2):1,(1,3):2,(2,1):3,(2,2):4,(3,1):5,(3,2):6,(4,1):7,(4,2):8}
  #split tree to chamber
    for idx in range(n_chambers):
        ec=chambers[idx][2] == '+'
        st=int(chambers[idx][3])
        rg=int(chambers[idx][5])
        ch=int(chambers[idx][7:])
        fitEffs=GetEff( "%s%s%s.root" % (Prefix,TagProbeFitResult,chambers[idx]),"lct_effV"+Postfix,"seg_effV"+Postfix )
        iBin_y=RingToYMap[(st,rg)]
        iBin_y=10+iBin_y if ec else 9-iBin_y
        SEGEff.SetBinContent(ch,iBin_y,fitEffs[0]*100.);
        SEGEff_downErr.SetBinContent(ch,iBin_y,fitEffs[1]*100.);
        SEGEff_upErr.SetBinContent(ch,iBin_y,fitEffs[2]*100.);
        LCTEff.SetBinContent(ch,iBin_y,fitEffs[3]*100.);
        LCTEff_downErr.SetBinContent(ch,iBin_y,fitEffs[4]*100.);
        LCTEff_upErr.SetBinContent(ch,iBin_y,fitEffs[5]*100.);
    SegCanvas=TCanvas("segment efficiency","segment efficiency",1500,1000)
    SegCanvas.cd()
    SEGEff.Draw("COLZ,TEXT")
    SEGEff_upErr.Draw("TEXT,SAME")
    SEGEff_downErr.Draw("TEXT,SAME")
    SEGEff.SetMaximum(100)
    SEGEff.SetMinimum(0)
    LCTCanvas=TCanvas("lct efficiency","lct efficiency",1500,1000)
    LCTCanvas.cd()
    LCTEff.Draw("COLZ,TEXT")
    LCTEff_upErr.Draw("TEXT,SAME")
    LCTEff_downErr.Draw("TEXT,SAME")
    LCTEff.SetMaximum(100)
    LCTEff.SetMinimum(0)
    file_out.cd()
    SegCanvas.Write()
    LCTCanvas.Write()
    SEGEff.Write()
    SEGEff_upErr.Write()
    SEGEff_downErr.Write()
    LCTEff.Write()
    LCTEff_upErr.Write()
    LCTEff_downErr.Write()
raw_input("Plots are saved in "+ResultPlotsFileName+". Press ENTER to exit")
file_out.Close()
