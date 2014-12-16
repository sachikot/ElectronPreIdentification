import FWCore.ParameterSet.Config as cms

process = cms.Process("RERECO")

#import of standard configuration
# process.load('Configuration.StandardSequences.Services_cff')
# process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
# process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
# process.load('Configuration.StandardSequences.EndOfProcess_cff')
# process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# process.load('TrackingTools.GsfTracking.FwdAnalyticalPropagator_cfi')
# process.load('RecoParticleFlow.PFTracking.trackerDrivenElectronSeeds_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(50)
    )

#Input source 
process.source = cms.Source("PoolSource",
#   secondaryFileNames = cms.untracked.vstring(),
   fileNames = cms.untracked.vstring(
#    INPUTFILE

#        "root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/QCD_Pt_30to80_bcToE_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0A802DC9-296F-E411-9107-0025905964A6.root"
        "root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0432E62A-7A6C-E411-87BB-002590DB92A8.root"
    )
)

#process.options = cms.untracked.PSet()

# Production Info
# process.configurationMetadata = cms.untracked.PSet(
#     version = cms.untracked.string('$Revision: 1.5 $'),
#     annotation = cms.untracked.string('reco nevts:1'),
#     name = cms.untracked.string('Applications')
# )


# Other statements
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')

# Path and EndPath definitions
#process.reconstruction_step = cms.Path(process.reconstruction_fromRECO)
#process.endjob_step = cms.EndPath(process.endOfProcess)

process.demo = cms.EDAnalyzer(
 'ElePFAnalyzer',
 genParticles = cms.InputTag("prunedGenParticles"),
 gsfTracks = cms.InputTag("gsfElectrons")
)

process.TFileService = cms.Service(
 "TFileService",
 fileName = cms.string('test_DY.root') 
)


process.run_analyzer = cms.Path(process.demo)

# Schedule definition
#process.schedule = cms.Schedule(
#    process.reconstruction_step,
#    process.make_pftracks,
#    process.run_analyzer
#    )



