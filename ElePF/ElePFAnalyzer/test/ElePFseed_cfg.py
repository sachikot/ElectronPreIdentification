import FWCore.ParameterSet.Config as cms

process = cms.Process("RERECO")

#import of standard configuration
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('TrackingTools.GsfTracking.FwdAnalyticalPropagator_cfi')
process.load('RecoParticleFlow.PFTracking.trackerDrivenElectronSeeds_cfi')
process.load('SimTracker.TrackAssociation.TrackAssociatorByHits_cfi')

#from SimTracker.TrackAssociation.TrackAssociatorByHits_cfi import *
process.TrackAssociatorByHits.Quality_SimToReco = cms.double(0.01)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5)
    )

#Input source 
process.source = cms.Source("PoolSource",
   secondaryFileNames = cms.untracked.vstring(),
   fileNames = cms.untracked.vstring(
    "root://cmssrv32.fnal.gov//store/relval/CMSSW_7_0_0_pre11/RelValQCD_Pt_30_80_BCtoE_8TeV/GEN-SIM-DIGI-RAW-HLTDEBUG/PU50ns_START70_V4-v6/00000//BC6832CE-6283-E311-851C-002618943911.root"
    )
)

process.options = cms.untracked.PSet()

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.5 $'),
    annotation = cms.untracked.string('reco nevts:1'),
    name = cms.untracked.string('Applications')
)


# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')

#added by Sachiko
process.MessageLogger = cms.Service("MessageLogger",
     default = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
)

process.load("SimGeneral.TrackingAnalysis.simHitTPAssociation_cfi");

# Path and EndPath definitions
process.reconstruction_step = cms.Path(process.reconstruction_fromRECO)
process.endjob_step = cms.EndPath(process.endOfProcess)

process.load("RecoParticleFlow.PFTracking.particleFlowTrack_cff")
process.make_pftracks = cms.Path(process.pfTrackingGlobalReco)

process.demo = cms.EDAnalyzer(
 'ElePFAnalyzer',
 simtracksTag = cms.InputTag("g4SimHits"),
 tracksTag = cms.InputTag("generalTracks"),
 label_tp=cms.InputTag("mix","MergedTrackTruth"),
 EtaMap = cms.string('RecoParticleFlow/PFBlockProducer/data/resmap_ECAL_eta.dat'),
 PhiMap = cms.string('RecoParticleFlow/PFBlockProducer/data/resmap_ECAL_phi.dat'),
 HcalWindow=cms.double(0.184),
 GsfTrackModuleLabel = cms.InputTag("gsfElectrons"),
 PFRecTrackLabel = cms.InputTag("particleFlowRecHitECAL"),
 PFNuclear = cms.InputTag("particleFlowBlock"),
 useNuclear = cms.bool(False),
 TkColList = cms.VInputTag(cms.InputTag("generalTracks")),
 UseQuality = cms.bool(True),
 TrackQuality = cms.string('highPurity'),
 Smoother = cms.string('GsfTrajectorySmoother_forPreId'),
 Fitter = cms.string('GsfTrajectoryFitter_forPreId'),
 PFEcalClusterLabel = cms.InputTag("particleFlowClusterECAL"),
 PFHcalClusterLabel = cms.InputTag("particleFlowClusterHCAL"),
 PSThresholdFile = cms.string('RecoParticleFlow/PFTracking/data/PSThreshold.dat'),
 ClusterThreshold = cms.double(0.5),
 MinEOverP = cms.double(0.3),
 MaxEOverP = cms.double(3.0),
 MinPt = cms.double(2.0),
 MaxPt = cms.double(50.0),
 MaxEta = cms.double(2.4),
 UsePreShower =cms.bool(False),
 ApplyIsolation = cms.bool(False),
 EcalStripSumE_deltaPhiOverQ_minValue = cms.double(-0.1),
 EcalStripSumE_minClusEnergy = cms.double(0.1),
 EcalStripSumE_deltaEta = cms.double(0.03),
 EcalStripSumE_deltaPhiOverQ_maxValue = cms.double(0.5),
 PFPSClusterLabel = cms.InputTag("particleFlowClusterPS"),
 EOverPLead_minValue = cms.double(0.95),
 HOverPLead_maxValue = cms.double(0.05),
 ThresholdFile = cms.string('RecoParticleFlow/PFTracking/data/Threshold.dat'),
 UseTMVA = cms.untracked.bool(True),
 filename = cms.string("test.root"),
 TMVAMethod = cms.string('BDT'),
 PreCkfLabel = cms.string('SeedsForCkf'),
 PreGsfLabel = cms.string('SeedsForGsf'),
 PreIdLabel = cms.string('preid'),
 NHitsInSeed = cms.int32(3),
 ProducePreId = cms.untracked.bool(True)
)


process.eana = cms.EDAnalyzer("EventContentAnalyzer")

process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)


process.run_analyzer = cms.Path(#process.eana*
                                #process.printTree*
                                process.simHitTPAssocProducer*
                                process.demo
                                )


# Schedule definition
process.schedule = cms.Schedule(
    process.raw2digi_step,
    process.L1Reco_step,
    process.reconstruction_step,
    process.make_pftracks,
    process.run_analyzer
    )



