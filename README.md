# CSC Primitive Efficiency Measurement Package
====================

## About
--- tested in CMSSW_5_x_x, not tested in CMSSW_6_x_x

## Table of Contents
- [Table of Contents](#table-of-contents)
    - [About](#about)
    - [Installation](#installation)
    - [Make the Ntuple](#make-the-ntuple)
    - [Example](#example)
    - [Limitations](#limitations)
   
## Installation:
<pre>
cd CMSSW_5_x_x
cmsenv
mkdir CSCEfficiency
cd CSCEfficiency
git clone git@github.com:zhang8473/CSCEfficiency.git
scramv1 b
</pre>

## Make the Ntuple
1. Config the variable `datatype` in [NtupleMaker.py](https://github.com/zhang8473/CSCEfficiency/blob/master/NtupleMaker.py): 
<pre>datatype="RAW"#
#Candidates are for data: "RAW" "RAW-RECO" "FEVT"
#for mc: in order of suggestions: "GEN-RAWDEBUG"(mc) "GEN-SIM-RAW"(mc) "GEN-RAW"(mc) "GEN-SIM"
</pre>
The The default output file name is 'CSCPFG_Ineff_DATA.root' ---
<pre>
process.aoddump.rootFileName=cms.untracked.string('CSCPFG_Ineff_DATA.root')
</pre>

2. Run [NtupleMaker.py](https://github.com/zhang8473/CSCEfficiency/blob/master/NtupleMaker.py) using Crab or locally. After all jobs finished, combine the output root files into one:
<pre>
hadd CSCPFG_Ineff_DATA.root CSCPFG_Ineff_DATA*.root
</pre>

## Make the efficiency plots
1. Setup  [Config.py](https://github.com/zhang8473/CSCEfficiency/blob/master/NtupleAnzScripts/Config.py).
   1. Data or MC:
   <pre>
    RunOnMC=False # or True
   </pre>
   2. Use Z resonance or J/Psi resonance:
   <pre>
   Resonance="Z" # or "JPsi"
   </pre>
   3. How to categorize the data:
   <pre>
   Group="Stations" #x axis: stations; y axis: efficiency
   #options are
   </pre>
   <table style="width:100%"  align="center">
    <tr>
      <th> Group </th><th> x axis </th><th>y axis</th><th>z axis</th><th> Comments</th>
    </tr>
    <tr>
      <td>"Stations"</td><td>stations</td><td>efficiency</td><td></td>
    </tr>
    <tr>
      <td>"Chambers"</td><td>chamber number</td><td>stations</td><td>efficiency</td>
    </tr>
    <tr>
      <td>"Stationspt"</td><td>pt in GeV</td><td>efficiency</td><td></td><td>make a plot for each station</td>
    </tr>
    <tr>
      <td>"Stationseta"</td><td>|η|</td><td>efficiency</td><td></td><td>make a plot for each station</td>
    </tr>
    <tr>
      <td>"Stationsphi"</td><td>ϕ</td><td>efficiency</td><td></td><td>make a plot for each station</td>
    </tr>
    <tr>
      <td>"pt","eta", or "phi"</td><td>pt,|η|,ϕ</td><td>efficiency</td><td></td><td>make one plot for all stations</td>
    </tr>
   </table> 
   4. Temporary output file (you may change the path but not the file name): 
      It may take two times the size of the Ntuple file space. The variable is `TemporaryOutputFile`. By default,        it will use the linux temporary path: /tmp/.
   5. Tag and probe file (do not need to change): The variable is `TagProbeFitResult`.
   6. result file (do not need to change): The variable is `ResultPlotsFileName`.
   7. For simulation sample, setup the pileup weight file `DataPileupRootFileName`. This is the file made by [estimatePileup2.py](https://cmssdt.cern.ch/SDT/lxr/source/RecoLuminosity/LumiDB/scripts/estimatePileup2.py) in CMSSW. The `pileup_mc` is the pileup weight used to generate the simulation sample, one can get it directly from the CMSSW package, for example:
   <pre>
   from SimGeneral.MixingModule.mix_2012_Startup_50ns_PoissonOOTPU_cfi import mix
   pileup_mc=mix.input.nbPileupEvents.probValue
   </pre>
  8. Station categorizing method: variable `station` is a python dictionary. The format of each component in the dictionary is "key(index):(logic expression in C style,name,color,station number)". e.g.,
   <pre>
   stations={
    ......
    2:("( CSCRg1==1 )","ME11B",kRed-9,1),
    ......}
   </pre>
2. <pre> python Step1.py Ntuple.root </pre>
3. Wait until all jobs finished. Use `ps -f` to check.
4. <pre> python Step2_PlotAll.py </pre>
5. Plots are in the result root file

## Organizing Result Plots:
* [DATAMCPlot.py](https://github.com/zhang8473/CSCEfficiency/blob/master/NtupleAnzScripts/SysCal_ExpertsOnly/DATAMCPlot.py)

## Distance Matching Study
