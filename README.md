# EleTrkAna

export SCRAM_ARCH=slc7_amd64_gcc700

cmsrel CMSSW_10_6_26

cd CMSSW_10_6_26/src/

cmsenv

git cms-init

git clone git@github.com:cms-egamma/EgammaPostRecoTools.git  EgammaUser/EgammaPostRecoTools

cd  EgammaUser/EgammaPostRecoTools

git checkout master

cd -

git clone git@github.com:zhangcg123/EleTrkAna.git

scram b -j 8
