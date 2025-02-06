# for json in c++
cd ../vcpkg
./bootstrap-vcpkg.sh
./vcpkg install nlohmann-json
cd ../vbs_cross_terms_study

setupATLAS # defined in startup bash script
lsetup git

asetup 23.6.26,AthGeneration # el9
source setupRivet

#lsetup "panda"
#lsetup "rucio -w"
#voms-proxy-init -voms atlas --valid 48:0 # 48h grid certif and not 12h defaulti
