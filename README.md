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
hadd Ntuple.root AllRetrievedJobs*.root
</pre>

## Make station efficiency plots using Z peak


