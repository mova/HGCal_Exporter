# If on DESY NFS, use inputFiles:/nfs/dust/cms/user/korcariw/generation/Guns/eGun10k/step3.root

cmsRun EDAnalyzers/TreeMaker/python/ConfFile_cfg.py \
    inputFiles=/afs/cern.ch/work/w/wkorcari/public/eGun10k/step3.root \
    genEleFilter=0 \
    genPartonFilter=0 \
    isGunSample=1 \
    onRaw=0 \
    storeSimHit=1 \
    storeRecHit=1 \
    debugFile=0 \
    maxEvents=100 \
