# CSC Primitive Efficiency Measurement Package
====================

## About
--- tested in CMSSW_5_x_x, not tested in CMSSW_6_x_x
* It is based on the tag-and-probe method using the Z pole or the J/ψ pole;
* The efficiency obtained is the CSC detector efficiency times the efficiency that the muon is not scattered.

## Table of Contents
- [Table of Contents](#table-of-contents)
    - [About](#about)
    - [Installation](#installation)
    - [Make the Ntuple](#make-the-ntuple)
    - [Make the Efficiency Plots](#make-the-efficiency-plots)
    - [Organize the Result Plots](#organize-the-result-plots)
    - [Study the Variables in the Ntuple](#study-the-variables-in-the-ntuple)
   
## Installation
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

## Make the Efficiency plots
1. Setup  [Config.py](https://github.com/zhang8473/CSCEfficiency/blob/master/NtupleAnzScripts/Config.py).
   1. Setup for Data or for MC:
   <pre>
    RunOnMC=False # or True
   </pre>
   2. Setup using Z resonance or using J/ψ resonance:
   <pre>
   Resonance="Z" # or "JPsi"
   </pre>
   3. Setup how to categorize the data:
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
   4. Arrange space for the temporary file (you may change the path but not the file name): 
      It may take two times the size of the Ntuple file space. The variable is `TemporaryOutputFile`. By default,        it will use the linux temporary path: /tmp/.
   5. Tag and probe file (do not need to change): The variable is `TagProbeFitResult`.
   6. result file (do not need to change): The variable is `ResultPlotsFileName`.
   7. **Setup the pileup reweighting scheme for simulation sample**: The variable of the data pileup weight file is  `DataPileupRootFileName`. This is the file made by [estimatePileup2.py](https://cmssdt.cern.ch/SDT/lxr/source/RecoLuminosity/LumiDB/scripts/estimatePileup2.py) in CMSSW. One can see [the twiki page](https://twiki.cern.ch/twiki/bin/view/CMS/PileupMCReweightingUtilities) to learn how to make such file corresponding to data. The `pileup_mc` is the pileup weight used to generate the simulation sample, one can get it directly from the CMSSW package, for example:
   <pre>
   from SimGeneral.MixingModule.mix_2012_Startup_50ns_PoissonOOTPU_cfi import mix
   pileup_mc=mix.input.nbPileupEvents.probValue
   </pre>
   The mixing module used for each MC sample is described in [the pileup information twiki page](https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupInformation).
  8. Station categorizing method: variable `station` is a python dictionary. The format of each component in the dictionary is "key(index):(logic expression in C style,name,color,station number)". e.g.,
   <pre>
   stations={
    ......
    2:("( CSCRg1==1 )","ME11B",kRed-9,1),
    ......}
   </pre>
2. Categorize the data and run the tag-and-probe package in CMSSW:
   <pre> python Step1.py Ntuple.root </pre>
3. Wait until all jobs finished. Use `ps -f` to check.
4. Make the plot:
   <pre> python Step2_PlotAll.py </pre>
   Advanced Usage of [Step2_PlotAll.py](https://github.com/zhang8473/CSCEfficiency/blob/master/NtupleAnzScripts/Step2_PlotAll.py):
   <pre> python Step2_PlotAll.py arg1 arg2 </pre>
   * arg1 is the name of the directory that stores the TagandProbe result files;
   * arg2: "lct_effV"+arg2 and "seg_effV"+arg2 are the root TDirectory name in the TagandProbe result root file. arg2 can be specified as "bkg" or "sig" for background and signal modeling;
   * Example1(plot default efficiencies): python Step2_PlotAll.py
   * Example2(for systematic -- bkg modeling): python Step2_PlotAll.py . bkg
   * Example3(for systematic -- sig modeling): python Step2_PlotAll.py . sig
   * Example4(MCTruth): python Step2_PlotAll.py ~/home/xxxxx/ mc
5. Plots are in the result root file

## Organize the Result Plots
To combine the data and MC results into one plot, one can use [DATAMCPlot.py](https://github.com/zhang8473/CSCEfficiency/blob/master/NtupleAnzScripts/SysCal_ExpertsOnly/DATAMCPlot.py).

## Study the Variables in the Ntuple
With the following python scipts, one can study the variables and their correlations in the Ntuple, e.g., the distance between the track and the LCT/segment.
* [MatchStudy.py](https://github.com/zhang8473/CSCEfficiency/blob/master/NtupleAnzScripts/SysCal_ExpertsOnly/MatchStudy.py) is for stations
* [MatchStudy_Chamber.py](https://github.com/zhang8473/CSCEfficiency/blob/master/NtupleAnzScripts/SysCal_ExpertsOnly/MatchStudy_Chamber.py) is for one or multiple chambers.

