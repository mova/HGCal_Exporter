import os

import FWCore.ParameterSet.Config as cms

processName = "Demo"

from Configuration.StandardSequences.Eras import eras
process = cms.Process(processName, eras.Phase2C9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "auto:phase2_realistic_T15", "")


process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49_cff')


############################## Parse arguments ##############################

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing("analysis")


options.register("sourceFile",
    "", # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.string, # string, int, or float
    "File containing list of input files" # Description
)

options.register("outputDir",
    "", # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.string, # string, int, or float
    "Output directory" # Description
)

options.register("outFileNumber",
    -1, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "File number (will be added to the filename if >= 0)" # Description
)

options.register("eventRange",
    "", # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.string, # string, int, or float
    "Syntax: Run1:Event1-Run2:Event2 (includes both)" # Description
)

options.register("debugFile",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Create debug file" # Description
)

options.register("onRaw",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Running on RAW" # Description
)

options.register("storeSimHit",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Store sim-hits" # Description
)

options.register("storeRecHit",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Store rec-hits" # Description
)

options.register("isGunSample",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Is it a particle gun sample" # Description
)

options.register("genEleFilter",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Apply gen-electron filter" # Description
)

options.register("genPhoFilter",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Apply gen-photon filter" # Description
)

options.register("genPartonFilter",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Apply gen-parton filter" # Description
)

options.register("trace",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Trace modules" # Description
)

options.register("memoryCheck",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Check memory usage" # Description
)

options.register("printTime",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Print timing information" # Description
)

options.register("depGraph",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Produce dependency graph only" # Description
)

options.parseArguments()


#maxEvents = -1
#options.maxEvents = 15
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))


sourceFile = "sourceFiles/SingleElectronFlatPtGun_fpantale_pT-0-200_eta-1p5-3p0_GEN-SIM-RECO/SingleElectronFlatPtGun_fpantale_pT-0-200_eta-1p5-3p0_GEN-SIM-RECO_mod.txt"


if (len(options.sourceFile)) :
    
    sourceFile = options.sourceFile


fNames = []

if (len(options.inputFiles)) :
    
    fNames = options.inputFiles

else :
    
    with open(sourceFile) as f:
        
        fNames = f.readlines()


for iFile, fName in enumerate(fNames) :
    
    if (
        "file:" not in fName and
        "root:" not in fName
    ) :
        
        fNames[iFile] = "file:%s" %(fName)


outFileSuffix = ""


if (options.onRaw) :
    
    outFileSuffix = "%s_onRaw" %(outFileSuffix)


if (options.outFileNumber >= 0) :
    
    outFileSuffix = "%s_%d" %(outFileSuffix, options.outFileNumber)


outFile = "ntupleTree%s.root" %(outFileSuffix)

if (len(options.outputDir)) :
    
    os.system("mkdir -p %s" %(options.outputDir))
    
    outFile = "%s/%s" %(options.outputDir, outFile)


sourceFileNames = cms.untracked.vstring(fNames)
#print sourceFileNames

process.source = cms.Source("PoolSource",
    fileNames = sourceFileNames,
    
    # Run1:Event1 to Run2:Event2
    #eventsToProcess = cms.untracked.VEventRange("1:78722-1:78722"),
    
    #duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
)


if (len(options.eventRange)) :
    
    process.source.eventsToProcess = cms.untracked.VEventRange(options.eventRange)


#process.options = cms.untracked.PSet(
#    #SkipEvent = cms.untracked.vstring("ProductNotFound"),
#    
#    #printDependencies = cms.untracked.bool(True),
#)


if (options.depGraph) :
    
    process.DependencyGraph = cms.Service("DependencyGraph")
    process.source = cms.Source("EmptySource")
    process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(0))



process.treeMaker = cms.EDAnalyzer(
    "TreeMaker",
    
    ############################## My stuff ##############################
    debug = cms.bool(False),
    
    isGunSample = cms.bool(bool(options.isGunSample)),
    
    storeSimHit = cms.bool(bool(options.storeSimHit)),
    storeRecHit = cms.bool(bool(options.storeRecHit)),
    
    
    ############################## GEN ##############################
    
    label_generator = cms.InputTag("generator"),
    label_genParticle = cms.InputTag("genParticles"),
    
    
    ############################## RECO ##############################
    
    label_pileup = cms.InputTag("addPileupInfo"),
    label_rho = cms.InputTag("fixedGridRhoFastjetAll"),
    
    label_HGCEESimHit = cms.InputTag("g4SimHits", "HGCHitsEE"),
    #label_HGCEESimHit = cms.InputTag("g4SimHits", "HGCHitsEE", "SIM"),
    
    label_HGCEERecHit = cms.InputTag("HGCalRecHit", "HGCEERecHits"),
    label_HGCHEFRecHit = cms.InputTag("HGCalRecHit", "HGCHEFRecHits"),
    label_HGCHEBRecHit = cms.InputTag("HGCalRecHit", "HGCHEBRecHits"),
    
    label_PFRecHit = cms.InputTag("particleFlowRecHitHGC", "Cleaned"),
    
    label_simCluster = cms.InputTag("mix", "MergedCaloTruth"),
    
    label_caloParticle = cms.InputTag("mix", "MergedCaloTruth"),
)


########## Filters ##########

from EDFilters.MyFilters.GenParticleFilter_cfi import *

# Gen-ele filter
process.GenParticleFilter_ele = GenParticleFilter.clone()
process.GenParticleFilter_ele.atLeastN = cms.int32(1)
process.GenParticleFilter_ele.pdgIds = cms.vint32(11)
process.GenParticleFilter_ele.minPt = cms.double(10)
process.GenParticleFilter_ele.minEta = cms.double(1.479)
process.GenParticleFilter_ele.maxEta = cms.double(3.1)
process.GenParticleFilter_ele.isGunSample = cms.bool(bool(options.isGunSample))
#process.GenParticleFilter_ele.debug = cms.bool(True)

process.filter_seq_genEle = cms.Sequence()

if (options.genEleFilter) :

    process.filter_seq_genEle = cms.Sequence(process.GenParticleFilter_ele)


# Gen-pho filter
process.GenParticleFilter_pho = GenParticleFilter.clone()
process.GenParticleFilter_pho.atLeastN = cms.int32(1)
process.GenParticleFilter_pho.pdgIds = cms.vint32(22)
process.GenParticleFilter_pho.minPt = cms.double(10)
process.GenParticleFilter_pho.minEta = cms.double(1.479)
process.GenParticleFilter_pho.maxEta = cms.double(3.1)
process.GenParticleFilter_pho.isGunSample = cms.bool(bool(options.isGunSample))
#process.GenParticleFilter_pho.debug = cms.bool(True)

process.filter_seq_genPho = cms.Sequence()

if (options.genPhoFilter) :

    process.filter_seq_genPho = cms.Sequence(process.GenParticleFilter_pho)


# Gen-parton filter
process.GenParticleFilter_part = GenParticleFilter.clone()
process.GenParticleFilter_part.atLeastN = cms.int32(1)
process.GenParticleFilter_part.pdgIds = cms.vint32(1, 2, 3, 4, 5, 21)
process.GenParticleFilter_part.minPt = cms.double(10)
process.GenParticleFilter_part.minEta = cms.double(1.479)
process.GenParticleFilter_part.maxEta = cms.double(3.1)
process.GenParticleFilter_part.isGunSample = cms.bool(bool(options.isGunSample))
#process.GenParticleFilter_part.debug = cms.bool(True)

process.filter_seq_genPart = cms.Sequence()

if (options.genPartonFilter) :

    process.filter_seq_genPart = cms.Sequence(process.GenParticleFilter_part)



#process.filter_path = cms.Path(
#    process.filter_seq_genEle *
#    process.filter_seq_genPart
#)


print "Deleting existing output file."
os.system("rm %s" %(outFile))


# Output file name modification
if (outFile.find("/eos/cms") ==  0) :
    
    outFile = outFile.replace("/eos/cms", "root://eoscms.cern.ch//eos/cms")


# Output
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(outFile)
)


process.schedule = cms.Schedule()


# Aging
from SLHCUpgradeSimulations.Configuration.aging import customise_aging_1000
customise_aging_1000(process)


process.reco_seq = cms.Sequence()

if (options.onRaw) :
    
    process.reco_seq = cms.Sequence(
        process.RawToDigi *
        process.L1Reco *
        process.reconstruction_mod
    )


###### PixelCPE issue
process.TrackProducer.TTRHBuilder = "WithTrackAngle"
process.PixelCPEGenericESProducer.UseErrorsFromTemplates = False
process.PixelCPEGenericESProducer.LoadTemplatesFromDB = False
process.PixelCPEGenericESProducer.TruncatePixelCharge = False
process.PixelCPEGenericESProducer.IrradiationBiasCorrection = False
process.PixelCPEGenericESProducer.DoCosmics = False
process.PixelCPEGenericESProducer.Upgrade = cms.bool(True) 
######


process.p = cms.Path(
    #process.filter_seq *
    
    process.filter_seq_genEle *
    process.filter_seq_genPho *
    process.filter_seq_genPart *
    
    process.reco_seq *
    
    process.treeMaker
)

process.schedule.insert(0, process.p)

print "\n"
print "*"*50
print "process.schedule:", process.schedule
print "*"*50
#print "process.schedule.__dict__:", process.schedule.__dict__
#print "*"*50
print "\n"


# Tracer
if (options.trace) :
    
    process.Tracer = cms.Service("Tracer")


if (options.memoryCheck) :
    
    process.SimpleMemoryCheck = cms.Service(
        "SimpleMemoryCheck",
        moduleMemorySummary = cms.untracked.bool(True),
    )


#Timing
if (options.printTime) :

    process.Timing = cms.Service("Timing",
        summaryOnly = cms.untracked.bool(False),
        useJobReport = cms.untracked.bool(True)
    )


# Debug
if (options.debugFile) :
    
    process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string("debug.root")
    )
    
    process.output_step = cms.EndPath(process.out)
    process.schedule.extend([process.output_step])


process.MessageLogger = cms.Service(
    "MessageLogger",
    destinations = cms.untracked.vstring(
        "cerr",
    ),
    cerr = cms.untracked.PSet(
        #threshold  = cms.untracked.string("ERROR"),
        DEBUG = cms.untracked.PSet(limit = cms.untracked.int32(0)),
        WARNING = cms.untracked.PSet(limit = cms.untracked.int32(0)),
        ERROR = cms.untracked.PSet(limit = cms.untracked.int32(0)),
    )
)


#from FWCore.ParameterSet.Utilities import convertToUnscheduled
#process = convertToUnscheduled(process)


# Add early deletion of temporary data products to reduce peak memory need
#from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
#process = customiseEarlyDelete(process)
# End adding early deletion
