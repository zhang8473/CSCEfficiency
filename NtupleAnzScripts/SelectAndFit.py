#!/usr/bin/python
#Author: Jinzhong Zhang(zhangjin@cern.ch), Northeastern University, Boston, U.S.A
#This script should be run by Step1.py
#Usage: python SelectAndFit.py ["debugfit","fit","nofit"] [ TemporaryOutputFile (optional) ] [ chambermin,chambermax (optinal) ]
#Example1 (select and fit all): python SelectAndFit.py fit
#Example2 (select and fit all): python SelectAndFit.py fit SelectedTrees.root
#Example3 (caculate counting efficiency only): python SelectAndFit.py nofit
#Example4 (refit existing selected mu pairs for chambers): python SelectAndFit.py debugfit AnyName(won't use).root 1,100
#debugfit: The ClassifiedTreeFile exists. Only implement TagProbeFitTreeAnalyzer.
#fit: The TemporaryOutputFile exists. After classified the tree, implement TagProbeFitTreeAnalyzer.
#nofit: The TemporaryOutputFile exists. Only classify the tree.

from  ROOT import *
from  numpy import *
from Config import *
gROOT.SetStyle("Plain")

import sys,os,gc
if (sys.argv[0] == "python"): args=sys.argv[2:]
else: args=sys.argv[1:]

BadChambers_=[]

def TreeFitter():
  #I have to do it in this way to avoid python memory allocation problems
  TmpFileName="tmp_Fitters_Ch%d_Ch%d.sh" % (chambermin,chambermax)
  f = open(TmpFileName, 'w')
  f.write("#!/bin/bash\n")
  if AnalyzeStations:
    for st in range(1,n_stations+1):
      f.write( "nohup cmsRun TagandProbe.py %s %s &\n"%(stations[st][1],ClassifiedTreeFile) )
      f.write( "sleep 3\n" )
      print "Fitting station ",stations[st][1],"..."
  else:
    for idx in range(chambermin,chambermax+1):
      if idx not in BadChambers_ :
        if os.path.isfile("%s%s.root" % (TagProbeFitResult,chambers[idx]) ):
          print "%s%s.root" % (TagProbeFitResult,chambers[idx]),"exists. Delete the old file first"
        else:
	  command="cmsRun TagandProbe.py %s %s &"%(chambers[idx],ClassifiedTreeFile)
          f.write( command+"\n" )
          f.write( "sleep 1\n" )
          print "Added \"", command, "\" to tmp_Fitters.sh"
      else:
        print "Skip bad chamber ",chambers[idx]
  f.write("rm "+TmpFileName)
  f.close()
  os.chmod( TmpFileName,0o777)# set file executable by all
  if os.system( "./"+TmpFileName ) != 0:#typically -1 for not enough memory
    print "\033[93mPlease run the following commands manually:\033[0m"
    print "./"+TmpFileName
    print "\033[93mThat is because of a memory problem.\033[0m"
  else:
    print "\033[92mRunning Processes:\033[0m"
    os.system( "ps -f")

def GetCountingEff(tree_,probe):#MCTruth Counting Efficiency
  passed = TH1D('passed','0',1,0,2)
  total = TH1D('total','2',1,0,2)
  total_eq = TH1D('total_eq','3',1,0,2)#total equivalent number of events
  tree_.Project('passed','1.0','weight*('+probe+'&&'+MCTruth+')')
  tree_.Project('total','1.0','weight*'+MCTruth)
  tree_.Project('total_eq','1.0','weight*weight*'+MCTruth)
  total_int=total.Integral()
  if total_int>0:
    eff=passed.Integral()/total_int
  else:
    return [0,1.,0.]
  ntotal_eq = sqrt(total_eq.Integral())
  npassed_eq = ntotal_eq*eff
  effErr_up=TEfficiency.ClopperPearson(int(ntotal_eq),int(npassed_eq),0.683,True)-eff
  effErr_down=eff-TEfficiency.ClopperPearson(int(ntotal_eq),int(npassed_eq),0.683,False)
#  effErr_up=effErr_down=sqrt( eff*(1-eff)/ntotal_eq )
  return [eff,effErr_up,effErr_down]

def AnalyzeTree( tree_,Title,Denominator,st_ ):
  print "Treename:",Title,MuTrackPairCut+"&&"+Denominator.replace("#",str(st_))
  try:
    tree_out=tree_.CopyTree( MuTrackPairCut+"&&"+Denominator.replace("#",str(st_)) )
  except:
    print "\033[93mWarning:I got memory problem. I am going to skip this tree.\033[0m"
    return [-1.]*6
  if tree_out.GetEntries() > 0:
    print ", keep."
    tree_out.SetName(Title)
    tree_out.AutoSave()
    if RunOnMC:
      result=GetCountingEff( tree_out,ProbeSegment.replace("#",str(st_)) )
      result.extend( GetCountingEff( tree_out,ProbeLCT.replace("#",str(st_)) ) )
      del tree_out#release memory
      return result
    else:
      del tree_out#release memory
      return [0.]*6
  print ", deleted."
  del tree_out
  return [-1.]*6

chambermin=0
chambermax=n_chambers-1
if len(args)>1:
  TemporaryOutputFile=args[1]
  if len(args)>2:
    chambermin=int(args[2].split(',')[0])
    chambermax=int(args[2].split(',')[1])
    ClassifiedTreeFile=ClassifiedTreeFile.replace( ".root","_Ch%dToCh%d.root"%(chambermin,chambermax) )

if ( args[0] == "debugfit" ):
  if not AnalyzeStations:
    file_in = TFile.Open(ClassifiedTreeFile,"read")
    for idx in range(chambermin,chambermax+1):
      try:
        test=file_in.Get("Probes/"+chambers[idx])
        if test.GetEntries() == 0:
          BadChambers_.append(idx)
        del test
      except:
        BadChambers_.append(idx)
    file_in.Close()
    del file_in
  TreeFitter()
  sys.exit()

file_in = TFile.Open(TemporaryOutputFile,"read")
tree_in=file_in.Get("Fraction")
outputfile=TFile.Open(ClassifiedTreeFile,'RECREATE')
outputdir=outputfile.mkdir("Probes")
outputdir.cd()

if AnalyzeStations:
  countingEffs=[]
  for st in range(1,n_stations+1):
    countingEffs.append( AnalyzeTree( tree_in,stations[st][1],stations[st][0]+"&&"+DenominatorRequire,stations[st][3]) )
  if RunOnMC:
    countingEffs=array(countingEffs).transpose()*100.
    xval=array(range(1,n_stations+1))*1.0
    xerr=zeros(n_stations, dtype=float)
    SEGEff=TGraphAsymmErrors(n_stations, xval, array(countingEffs[0]), xerr, xerr, array(countingEffs[1]), array(countingEffs[2]))
    LCTEff=TGraphAsymmErrors(n_stations, xval, array(countingEffs[3]), xerr, xerr, array(countingEffs[4]), array(countingEffs[5]))
    for st in range(1,n_stations+1):
       binnum=SEGEff.GetXaxis().FindBin(st)
       SEGEff.GetXaxis().SetBinLabel( binnum,stations[st][1] )
       LCTEff.GetXaxis().SetBinLabel( binnum,stations[st][1] )
    resultplotfile=TFile.Open(ResultPlotsFileName.replace(".root","_Counting.root"),'RECREATE')
    resultplotfile.cd()
    SEGEff.Write("SEGEff_St")
    LCTEff.Write("LCTEff_St")
    resultplotfile.Close()
    print "Counting Efficiency Plots are saved in ", resultplotfile.GetName()     

else:
  if RunOnMC:
    SEGEff=TH2F("SEGEff","segment efficiency",36,1,37,18,-9,9);
    Chambers_  = ["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36"]
    Rings_ = ["ME-42","ME-41","ME-32","ME-31","ME-22","ME-21","ME-13","ME-12","ME-11","ME+11","ME+12","ME+13","ME+21","ME+22","ME+31","ME+32","ME+41","ME+42"]
    for ich in range(36):
      SEGEff.GetXaxis().SetBinLabel(ich+1,Chambers_[ich])
    for irg in range(18):
      SEGEff.GetYaxis().SetBinLabel(irg+1,Rings_[irg])
    SEGEff.SetTitle("Efficiency")
    SEGEff_upErr=SEGEff.Clone("SEGEff_upErr")
    SEGEff_downErr=SEGEff.Clone("SEGEff_downErr")
    LCTEff=SEGEff.Clone("LCTEff")
    LCTEff_upErr=LCTEff.Clone("LCTEff_upErr")
    LCTEff_downErr=LCTEff.Clone("LCTEff_downErr")
    RingToYMap={(1,1):0,(1,2):1,(1,3):2,(2,1):3,(2,2):4,(3,1):5,(3,2):6,(4,1):7,(4,2):8}
  #split tree to chamber
  for idx in range(chambermin,chambermax+1):
    ec_=chambers[idx][2] == '+'
    st_=int(chambers[idx][3])
    rg_=int(chambers[idx].split('_')[1])
    ch_=int(chambers[idx].split('_')[2])
    countingEffs=AnalyzeTree( tree_in, chambers[idx],
                              "%sCSCEndCapPlus && CSCRg%s==%s && CSCCh%s==%s &&"%("" if (ec_) else '!',st_,rg_,st_,ch_)+DenominatorRequire,
                              st_ )
    gc.collect()#release memory
    if countingEffs[0] < -0.5 :
      BadChambers_.append(idx)
    else:
      if RunOnMC:
        iBin_y=RingToYMap[(st_,rg_)]
        iBin_y=10+iBin_y if ec_ else 9-iBin_y
        SEGEff.SetBinContent(ch_,iBin_y,countingEffs[0]*100.);
        SEGEff_upErr.SetBinContent(ch_,iBin_y,countingEffs[1]*100.);
        SEGEff_downErr.SetBinContent(ch_,iBin_y,countingEffs[2]*100.);
        LCTEff.SetBinContent(ch_,iBin_y,countingEffs[3]*100.);
        LCTEff_upErr.SetBinContent(ch_,iBin_y,countingEffs[4]*100.);
        LCTEff_downErr.SetBinContent(ch_,iBin_y,countingEffs[5]*100.);
  if RunOnMC:
    resultplotfile=TFile.Open(ResultPlotsFileName.replace(".root","_Counting.root"),'RECREATE')
    resultplotfile.cd()
    SEGEff.Write()
    SEGEff_upErr.Write()
    SEGEff_downErr.Write()
    LCTEff.Write()
    LCTEff_upErr.Write()
    LCTEff_downErr.Write()
    resultplotfile.Close()
    print "Counting Efficiency Plots are saved in ", resultplotfile.GetName()     

print "Output to ", outputfile.GetName()     
outputfile.Close()
file_in.Close()

if args[0]!="nofit":
  TreeFitter()
