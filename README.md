# CSC Primitive Efficiency Measurement Package
====================

## About
--- tested in CMSSW_5_x_x, not tested in CMSSW_6_x_x
* It is based on the tag-and-probe method using the Z pole or the J/ψ pole;
* The efficiency obtained is the CSC detector efficiency times the efficiency that the muon is not scattered.
* Need RAW information to get the LCT efficiency. RECO or AOD sample is not good.

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
1. Config the variable `datatype` in [NtupleMaker.py](NtupleMaker.py): 
<pre>datatype="RAW"#
#Candidates are for data: "RAW" "RAW-RECO" "FEVT"
#for mc: in order of suggestions: "GEN-RAWDEBUG"(mc) "GEN-SIM-RAW"(mc) "GEN-RAW"(mc) "GEN-SIM"
</pre>
The default output file name is 'CSCPFG_Ineff_DATA.root' ---
<pre>
process.aoddump.rootFileName=cms.untracked.string('CSCPFG_Ineff_DATA.root')
</pre>

2. Run [NtupleMaker.py](NtupleMaker.py) using Crab or locally. After all jobs finished, combine the output root files into one:
<pre>
hadd Ntuple.root CSCPFG_Ineff_DATA*.root
</pre>

## Make the Efficiency plots
1. Setup  [Config.py](NtupleAnzScripts/Config.py).
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
      <td>"Chambers"</td><td>chamber number</td><td>rings</td><td>efficiency</td>
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
      <td>"pt","eta", or "phi"</td><td>pt,|η|,ϕ</td><td>efficiency</td><td></td><td>make one plot for all stations, not tested yet</td>
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
   <h1 id="AdvUsage">Advanced Usage of</h1> [Step2_PlotAll.py](NtupleAnzScripts/Step2_PlotAll.py):
   <pre> python Step2_PlotAll.py arg1 arg2 </pre>
   * arg1 is the name of the directory that stores the TagandProbe result files;
   * arg2 is the postfix of the root TDirectory name in the TagandProbe result root file, for lct, the TDirectory name is "lct_effV"+arg2 and for segment, the TDirectory name is "seg_effV"+arg2. Moreover, arg2 can also be specified as "bkg" or "sig" for background and signal modeling;
   * Example1(plot default efficiencies): python Step2_PlotAll.py
   * Example2(for systematic -- bkg modeling): python Step2_PlotAll.py . bkg
   * Example3(for systematic -- sig modeling): python Step2_PlotAll.py . sig
   * Example4(MCTruth): python Step2_PlotAll.py ~/home/xxxxx/ mc
5. Plots are in the result root file

## Organize the Result Plots
To combine the data and MC results into one plot, one can use [DATAMCPlot.py](NtupleAnzScripts/DATAMCPlot.py). It oragnizes the plots made by [Step2_PlotAll.py](NtupleAnzScripts/Step2_PlotAll.py). The usage is
<pre>
python DATAMCPlot.py datafile mcfile plotname
</pre>
* datafile is the result root file from data;
* mcfile is the result root file from simulation. If the keyword `MCTruth` appears in the file name, the simulation couting efficiency for real muons will be plotted. In that case, the plotname will be changed to plotname+"_MCTruth" automatically in the script. So one should still use the same plotname while calculating MCTruth.
* plotname is the name of the plot saved in the result root file, e.g. "ME12+13seg_effV" for segment efficiency or "ME12+13lct_effV" for lct efficiency.

I suggest to put the datafile and the mcfile in different directories. This script will use the  [Config.py](NtupleAnzScripts/Config.py) in the datafile directory. If no Config.py or no \__init\__.py is found in the datafile directory, it will use the Config.py in the current directory.

## Study the Variables in the Ntuple
This part is only for **experts** who want to find out a problem or know more. Here only list a breif discription for each script because **experts** are able to read the python script themselves. With the following python scipts, one can study the variables and their correlations in the Ntuple, e.g., the distance between the track and the LCT/segment.
* [MatchStudy.py](NtupleAnzScripts/ExpertsOnly/MatchStudy.py) can be used to study the variables in category of stations. While using this, the `Group` should be set to "Stations" in [Config.py](NtupleAnzScripts/Config.py).
* [MatchStudy_Chamber.py](NtupleAnzScripts/ExpertsOnly/MatchStudy_Chamber.py) can be used to study the variables in category of chambers. While using this, the `Group` should be set to "Chambers" in [Config.py](NtupleAnzScripts/Config.py).
* [RateStudy.py](NtupleAnzScripts/ExpertsOnly/RateStudy.py) can be used to study the simulation truth counting efficiency versus any variables in the Ntuple. The purpose is to find out the correlations between the efficiency and variables in the Ntuple.
* [Systematic1D.py](NtupleAnzScripts/ExpertsOnly/ Systematic1D.py) is to calculate the systematic uncertainties from different sources and get the final results with both systematic and statistic uncertainties. The input files are the plots made by [Step2_PlotAll.py](NtupleAnzScripts/Step2_PlotAll.py). See [the advanced usage](#AdvUsage) of [Step2_PlotAll.py](NtupleAnzScripts/Step2_PlotAll.py):
