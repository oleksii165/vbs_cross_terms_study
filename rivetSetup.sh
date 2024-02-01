export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh 

asetup main,latest,AthGeneration
source setupRivet

lsetup "panda 1.5.68"
lsetup "rucio -w"
voms-proxy-init -voms atlas --valid 48:0 # 48h grid certif and not 12h default

source /exp/atlas/kurdysh/vbs_cross_terms_study/python_packages/setup.sh # this is for sympy installed with pipInstall
source /exp/atlas/kurdysh/vbs_cross_terms_study/python_package2/setup.sh # for plotting modules of sympy
