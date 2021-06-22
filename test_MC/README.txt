This branch contains the necessary files to run over monte carlo and data samples. It includes the important features for the Displaced Tops analysis like :
- An ntuple creator with the information about tracks (gen, simu, and reco + associations) 
- Information about jets (AK4 and AK8 for calo and PF jets)
- Informations about muons, electrons, MET ... 
- Information about the bs, primary vertices + refit 
- Trigger bits (single muon)
- ZMu candidates (because of interest in ZMu skim data sample)

Recipe to run the ntuple :
mkdir NewEnv 
cd NewEnv 
cmsrel CMSSW_10_6_0
cd CMSSW_10_6_0/src
cmsenv 
git cms-init 
git pull https://github.com/doujadarej/cmssw.git NTuple_DisplacedTop 
git cms-addpkg TrackingPerf 
git cms-addpkg RecoJets/JetProducers 
git cms-addpkg DPGAnalysis/Skims 
#other packages may be needed if displaced tracking also wanted
git cms-addpkg test_MC 
git cms-addpkg test_data 

scramv1 b -j16 

to run on MC : 
cd test_MC 
#replace name of input file by the right one (e.g: using step2.root from benchmarks)
cmsRun reco_Ntuple_withZmuSkim.py

to run on data :
cd test_data 
cmsRun reReco_Data_Ntuple.py
