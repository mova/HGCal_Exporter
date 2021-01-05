```
cmsrel CMSSW_11_1_6
cd CMSSW_11_1_6/src/
cmsenv
git clone git@github.com:SohamBhattacharya/HGCal_GAN.git .
scram b -j8 USER_CXXFLAGS="-fopenmp"
./run.sh
```
