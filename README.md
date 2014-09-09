# CSCEfficiency
=============

## Make the Ntuple
1. Installation:
<pre>
cd CMSSW_5_x_x
cmsenv
mkdir CSCEfficiency
cd CSCEfficiency
git clone git@github.com:zhang8473/CSCEfficiency.git
scramv1 b
</pre>

2. Config the variable `datatype` in [NtupleMaker.py](https://github.com/zhang8473/CSCEfficiency/blob/master/NtupleMaker.py): 
<pre>datatype="RAW"#
#Candidates are data: "RAW" "RAW-RECO" "FEVT"
#mc: in order of suggestions: "GEN-RAWDEBUG"(mc) "GEN-SIM-RAW"(mc) "GEN-RAW"(mc) "GEN-SIM"
</pre>
The The default output file name is 'CSCPFG_Ineff_DATA.root' ---
<pre>
process.aoddump.rootFileName=cms.untracked.string('CSCPFG_Ineff_DATA.root')
</pre>

3. Run [NtupleMaker.py](https://github.com/zhang8473/CSCEfficiency/blob/master/NtupleMaker.py) using Crab or locally. After all jobs finished, combine the output root files into one:
<pre>
hadd CSCPFG_Ineff_DATA.root CSCPFG_Ineff_DATA*.root
</pre>

## Make the efficiency plots
1. Setup  [Config.py](https://github.com/zhang8473/CSCEfficiency/blob/master/NtupleAnzScripts/Config.py).
   1. Data or MC:
   <pre>
    RunOnMC=False
   </pre>
   2. Use Z resonance or J/Psi resonance:
   <pre>
   Resonance="Z"#options are "Z","JPsi"
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
      <td>"Stationseta"</td><td>|&eta|</td><td>efficiency</td><td></td><td>make a plot for each station</td>
    </tr>
   </table> 
   4. Temporary output file: 
      It may take two times the size of the Ntuple file space. The variable is `TemporaryOutputFile`. By default,        it will use the linux temporary path: /tmp/.
2. <pre> python Step1.py Ntuple.root </pre>
3. 

