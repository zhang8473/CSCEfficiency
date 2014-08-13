#!/usr/bin/python
#Author: Jinzhong Zhang(zhangjin@cern.ch), Northeastern University, Boston, U.S.A
#Usage: python Step1.5.py
#It only does the tagandprobe fit, in case that the classified Ntuple files have been made)"
from Config import *

def Fit(filename_,st_):
  if os.path.isfile(filename_):
      cmd="nohup cmsRun TagandProbe.py %s %d &"%(filename_,st_)
      subprocess.Popen([cmd], shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
      os.system("sleep 1")
      status = subprocess.Popen(["ps -f|grep -v 'grep'|grep '%s'"%(cmd[6:-2])], shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE).stdout.read()
      if status:
          print "\033[92m",status,"\033[0m"
      else:
          print "\033[91mThe job for",filename_,"is NOT running. Memory or disk full?\033[0m"
  else:
      print "\033[91m",filename_," does not exist.\033[0m"

if "Stations" in Group:
    for st in range(1,n_stations+1):
        Fit( TemporaryOutputFile[:-5]+stations[st][1]+'.root',stations[st][3] )
elif "Chambers" in Group:
    for idx in range(n_chambers):
        Fit(TemporaryOutputFile[:-5]+chambers[idx]+'.root',st)
else:
    Fit( TemporaryOutputFile[:-5]+"AllStations.root",0 )
