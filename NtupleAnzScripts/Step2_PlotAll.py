#!/usr/bin/python
#Author: Jinzhong Zhang(zhangjin@cern.ch), Northeastern University, Boston, U.S.A
#################################################################################
#Advanced Usage:
#python Step2_PlotAll.py arg1 arg2
#   arg1 is directory of the files given by the TagandProbe cmssw package. The file names have to match what is defined in Config.py;
#   arg2 "lct_effV"+arg2 and "seg_effV"+arg2 are the directory name in the TagandProbe result file;
#   arg2 can be specified as "bkg" or "sig" for background and signal modeling
#Example1(plot default efficiencies): python Step2_PlotAll.py
#Example2(systematic -- bkg modeling): python Step2_PlotAll.py /uscms/home/zhangjin/nobackup/ bkg
#Example3(systematic -- sig modeling): python Step2_PlotAll.py /uscms/home/zhangjin/nobackup/ sig
#Example4(systematic -- MCTruth     ): python Step2_PlotAll.py /uscms/home/zhangjin/nobackup/ mc
#Example2-4 are used for systematic calculation.
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
Length = arr('d',[0.00, 0.40, 0.60, 0.80, 1.00])
TColor.CreateGradientColorTable(5,Length,Red,Green,Blue,500)
gStyle.SetNumberContours(500)

import sys,os,re
if (sys.argv[0] == "python"): args=sys.argv[2:]
else: args=sys.argv[1:]

Prefix=""
Postfix=""
TagProbeFitResult=TagProbeFitResult.split("/")[-1]
ResultPlotsFileName=ResultPlotsFileName.split("/")[-1]
if len(args)>0:
    Prefix=args[0]
    if Prefix[-1] != "/":
        Prefix+="/"
    if len(args)>1:
        if args[1] == "bkg":
            Postfix="_BkgModeling"
        elif args[1] == "sig":
            Postfix="_SigModeling"
        elif args[1] == "mc":
            Postfix="_MCTruth"
        else:
            Postfix=args[1]
    ResultPlotsFileName=Prefix+ResultPlotsFileName.replace(".root",Postfix+".root")

file_out=TFile.Open(ResultPlotsFileName,'RECREATE')

etascheme="abseta"
#etascheme="tracks_eta"
phischeme="shiftedphi"
#phischeme="tracks_phi"
if "pt" in Group:
    binning="pt"
    plotname="tracks_pt_PLOT_"+etascheme+"_bin0_&_"+phischeme+"_bin0"
elif "eta" in Group:
    binning="eta"
    plotname=etascheme+"_PLOT_"+phischeme+"_bin0_&_tracks_pt_bin0"
elif "phi" in Group:
    binning="phi"
    plotname=phischeme+"_PLOT_"+etascheme+"_bin0_&_tracks_pt_bin0"
else:
    plotname=phischeme+"_bin0__"+etascheme+"_bin0__tracks_pt_bin0__VoigtianPlusExpo"

if Postfix=="_MCTruth":
    plotname+="_&_mcTrue_true"

def GetEff(f_in,path="lct_effV",effcat="fit_eff"):
    try:
        eff=f_in.Get("Probes/"+path+"/"+effcat).get().find("efficiency")
        return [eff.getVal(),abs(eff.getErrorLo()),eff.getErrorHi()]
    except:
        print "\033[97mOther problems, skip",f_in.GetName(),"\033[0m"
        return [0.]*3

def GetBinnedEffPlot(f_in,path="lct_effV",effcat="fit_eff",st_=0,name_=plotname):
    canvas_=f_in.Get("Probes/"+path+"/"+effcat+"_plots/"+name_)
    if not canvas_:
        print "\033[91m Warning: Probes/"+path+"/"+effcat+"_plots/"+name_," does not exist in",f_in.GetName(),"\033[0m"
        return NULL
    dummyplot_=canvas_.GetListOfPrimitives().At(0)
    plot_=canvas_.FindObject("hxy_"+effcat).Clone();
    #we are going to fix the bugs in tagandprobe package in the following code
    #1 - recreate the arrays
    nbins=plot_.GetN()
    xval=zeros(nbins, dtype=float)
    xerr=zeros(nbins, dtype=float)
    yerrhi=zeros(nbins, dtype=float)
    yerrlo=zeros(nbins, dtype=float)
    #2 - the y values are correct
    Y=plot_.GetY()
    #3 - find the corresponding first bin in the predefined bins for the plot_ first bin0
    exec( "bins=%sbin%s"%(binning,str(st_) if binning=="eta" else "") )
    X=plot_.GetX()
    for abin in bins:
        if X[0]<abin:
            firstbin=bins.index(abin)-1
            break
    #4 - fill the yerror bars from the correct input (only for fit efficiency)
    if effcat=="fit_eff":
        list_=f_in.Get("Probes/"+path).GetListOfKeys()
        ikey=list_.First()
        while (ikey!=list_.After(list_.Last())):
            dirname_=ikey.GetName()
            binnumber=re.match(".*"+binning+"_bin(\d*)_.*",dirname_)
            if binnumber:
                ibin=int(binnumber.group(1))-firstbin
                if ibin<nbins and ibin>=0:
                    result_=f_in.Get("Probes/"+path+"/"+dirname_+"/fitresults")
                    if result_:
                        value_=result_.floatParsFinal().find("efficiency")
                        yerrlo[ibin]=abs(value_.getErrorLo())
                        yerrhi[ibin]=value_.getErrorHi()
                        """
                        if Y[ibin]!=value_.getVal(): #show differences (should less than 1E-6)
                            print Y[ibin],
                            value_.Print()
                        """
                        if Y[ibin]<0.999 and yerrhi[ibin]<1E-7:
                            yerrhi[ibin]=yerrlo[ibin] #sometime the result misses the high error, we make an approximation: ErrorHi=ErrorLo in this case
                        if Y[ibin]+yerrhi[ibin]>1.:
                            yerrhi[ibin]=1.-Y[ibin] # it happens sometime when ErrorHi=ErrorLo
            ikey = list_.After(ikey);
    #5 - fill the correct x values from the binning 
    for ibin in range(nbins):
        xval[ibin]=(bins[ibin+firstbin]+bins[ibin+firstbin+1])/2.
        xerr[ibin]=abs(bins[ibin+firstbin+1]-bins[ibin+firstbin])/2.
    #6 - remake the TGraph
    plotname_=f_in.GetName().replace(Prefix+TagProbeFitResult,"")[:-5]+path
    outputplot=TGraphAsymmErrors(nbins, xval, Y, xerr, xerr, yerrlo, yerrhi)
    outputplot.SetName(plotname_)
    outputplot.SetTitle(outputplot.GetName())
    outputplot.GetXaxis().SetTitle(dummyplot_.GetXaxis().GetTitle())
    outputplot.GetYaxis().SetTitle(dummyplot_.GetYaxis().GetTitle())
    outputplot.GetYaxis().SetTitleOffset(1.2)
    #outputplot.SetMarkerStyle(8)
    #outputplot.SetMarkerSize(.5)
    """
    outputplot.SetMinimum(0.9)
    outputplot.SetMaximum(1.0)
    EffCanvas=TCanvas(plotname_,plotname_,500,500)
    EffCanvas.cd()
    if binning=="pt":
        outputplot.GetXaxis().SetLimits(10., 100.);
    #gPad.SetLogx()
    outputplot.Draw("AP")
    raw_input("pause")
    """
    return outputplot

if "Stations" in Group:
    Effs=[]
    for idx in range(1,n_stations+1):
        filename_=Prefix+TagProbeFitResult+stations[idx][1]+".root"
        if not os.path.isfile(filename_):
            print filename_+" is not found, skip.. "
            Effs.append([0.]*6)
            continue
        f_in=TFile(filename_,"READ");
        categoryname="cnt_eff" if Postfix=="_MCTruth" else "fit_eff"
        if "pt" in Group or "eta" in Group or "phi" in Group:
            LCTEff=GetBinnedEffPlot(f_in, "lct_effV"+Postfix,categoryname,stations[idx][3])
            SEGEff=GetBinnedEffPlot(f_in, "seg_effV"+Postfix,categoryname,stations[idx][3])
            file_out.cd()
            if LCTEff:
                LCTEff.Write()
            if SEGEff:
                SEGEff.Write()
        else:
            Effs.append( GetEff(f_in, "lct_effV"+Postfix,categoryname)+GetEff(f_in,"seg_effV"+Postfix,categoryname) )
            f_in.Close()
    if not ("pt" in Group or "eta" in Group or "phi" in Group):
        Effs=array(Effs).transpose()*100.
        xval=array(range(1,n_stations+1))*1.0
        xerr=zeros(n_stations, dtype=float)
        SEGEff=TGraphAsymmErrors(n_stations, xval, array(Effs[0]), xerr, xerr, array(Effs[1]), array(Effs[2]))
        LCTEff=TGraphAsymmErrors(n_stations, xval, array(Effs[3]), xerr, xerr, array(Effs[4]), array(Effs[5]))
        SegCanvas=TCanvas("segment efficiency","segment efficiency",500,500)
        SegCanvas.cd()
        SEGEff.SetMaximum(100)
        SEGEff.SetMinimum(90)
        LCTEff.SetMaximum(100)
        LCTEff.SetMinimum(90)
        SEGEff.SetMarkerStyle(8)
        SEGEff.SetMarkerSize(.5)
        SEGEff.Draw("AP")
        LCTCanvas=TCanvas("lct efficiency","lct efficiency",500,500)
        LCTCanvas.cd()
        LCTEff.SetMarkerStyle(8)
        LCTEff.SetMarkerSize(.5)
        LCTEff.Draw("AP")
        for st in range(1,n_stations+1):
           binnum=SEGEff.GetXaxis().FindBin(st)
           SEGEff.GetXaxis().SetBinLabel( binnum,stations[st][1] )
           LCTEff.GetXaxis().SetBinLabel( binnum,stations[st][1] )
        file_out.cd()
        if LCTEff:
           LCTEff.Write("LCTEff")
        if SEGEff:
           SEGEff.Write("SEGEff")
elif "Chambers" in Group:
    SEGEff=TH2F("SEGEff","segment efficiency",36,1,37,20,-9,9)
    SEGEff.SetMarkerSize(0.7)
    SEGEff.SetContour(500)
    SEGEff_upErr=TH2F("SEGEff_upErr","segment efficiency uperror",36,1,37,20,-8.7,9.3)
    SEGEff_upErr.SetMarkerSize(0.45)
    SEGEff_downErr=TH2F("SEGEff_downErr","segment efficiency loerror",36,1,37,20,-9.3,8.7)
    SEGEff_downErr.SetMarkerSize(0.45)
    SEGEff.GetYaxis().SetTickLength(0)
    Chambers_  = ["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36"]
    Rings_ = ["ME-42","ME-41","ME-32","ME-31","ME-22","ME-21","ME-13","ME-12","ME-11B","ME-11A","ME+11A","ME+11B","ME+12","ME+13","ME+21","ME+22","ME+31","ME+32","ME+41","ME+42"]
    for ich in range(36):
        SEGEff.GetXaxis().SetBinLabel(ich+1,Chambers_[ich])
        SEGEff_upErr.GetXaxis().SetBinLabel(ich+1,Chambers_[ich])
        SEGEff.GetXaxis().SetBinLabel(ich+1,Chambers_[ich])
    for irg in range(20):
        SEGEff.GetYaxis().SetBinLabel(irg+1,Rings_[irg])
        SEGEff_upErr.GetYaxis().SetBinLabel(irg+1,Rings_[irg])
        SEGEff_downErr.GetYaxis().SetBinLabel(irg+1,Rings_[irg])
    LCTEff=SEGEff.Clone("LCTEff")
    LCTEff.SetTitle("LCT efficiency")
    LCTEff_upErr=SEGEff_upErr.Clone("LCTEff_upErr")
    LCTEff_downErr=SEGEff_downErr.Clone("LCTEff_downErr")
    LCTEff.GetYaxis().SetTickLength(0)
    RingToYMap={(1,4):0,(1,1):1,(1,2):2,(1,3):3,(2,1):4,(2,2):5,(3,1):6,(3,2):7,(4,1):8,(4,2):9}
  #split tree to chamber
    for idx in range(n_chambers):
        ec=chambers[idx][2] == '+'
        st=int(chambers[idx][3])
        rg=int(chambers[idx][5])
        ch=int(chambers[idx][7:])
        filename_="%s%s%s.root" % (Prefix,TagProbeFitResult,chambers[idx])
        if not os.path.isfile(filename_):
            print filename_+" is not found, skip.. "
            Effs.append([0.]*6)
            continue
        f_in=TFile(filename_,"READ");
        if Postfix=="_MCTruth":
            Effs=GetEff(f_in, "lct_effV"+Postfix,"cnt_eff")+GetEff(f_in,"seg_effV"+Postfix,"cnt_eff")
        else:
            Effs=GetEff(f_in, "lct_effV"+Postfix,"fit_eff")+GetEff(f_in,"seg_effV"+Postfix,"fit_eff" )
        f_in.Close()
        iBin_y=RingToYMap[(st,rg)]
        iBin_y=11+iBin_y if ec else 10-iBin_y
        LCTEff.SetBinContent(ch,iBin_y,Effs[0]*100.);
        LCTEff_downErr.SetBinContent(ch,iBin_y,Effs[1]*100.);
        LCTEff_upErr.SetBinContent(ch,iBin_y,Effs[2]*100.);
        SEGEff.SetBinContent(ch,iBin_y,Effs[3]*100.);
        SEGEff_downErr.SetBinContent(ch,iBin_y,Effs[4]*100.);
        SEGEff_upErr.SetBinContent(ch,iBin_y,Effs[5]*100.);
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
elif "pt" in Group or "eta" in Group or "phi" in Group:
    filename_=Prefix+TagProbeFitResult+"AllStations.root"
    if not os.path.isfile(filename_):
        print filename_+" is not found, skip.. "
    else:
        f_in=TFile(filename_,"READ");
        categoryname="cnt_eff" if Postfix=="_MCTruth" else "fit_eff"
        LCTEff=GetBinnedEffPlot(f_in, "lct_effV"+Postfix,categoryname,stations[idx][3])
        SEGEff=GetBinnedEffPlot(f_in, "seg_effV"+Postfix,categoryname,stations[idx][3])
        f_in.Close()
        file_out.cd()
        if LCTEff:
           LCTEff.Write("LCTEff")
        if SEGEff:
           SEGEff.Write("SEGEff")
#raw_input("Plots are saved in "+ResultPlotsFileName+". Press ENTER to exit")
print "Plots are saved in",ResultPlotsFileName+"."
file_out.Close()
