# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --conditions auto:run2_mc --pileup_input das:/RelValMinBias_13/CMSSW_10_5_0_pre1-103X_mcRun2_asymptotic_v3-v1/GEN-SIM -n 10 --era Run2_2016 --eventcontent RECOSIM,MINIAODSIM,DQM --runUnscheduled -s RAW2DIGI,L1Reco,RECO,RECOSIM,EI,PAT,VALIDATION:@standardValidation+@miniAODValidation,DQM:@standardDQM+@ExtraHLT+@miniAODDQM --datatier GEN-SIM-RECO,MINIAODSIM,DQMIO --pileup AVE_35_BX_25ns --number=1 --filein file:step2.root --fileout file:step3.root
import FWCore.ParameterSet.Config as cms



from Configuration.Eras.Era_Run2_2018_cff import Run2_2018

process = cms.Process('RECO',Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.PATMC_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('DQMServices.Core.DQMStoreNonLegacy_cff')
process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(8)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:step2.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:step3.root'),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

process.MINIAODSIMoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('MINIAODSIM'),
        filterName = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    eventAutoFlushCompressedSize = cms.untracked.int32(-900),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string('file:step3_inMINIAODSIM.root'),
    outputCommands = process.MINIAODSIMEventContent.outputCommands,
    overrideBranchesSplitLevel = cms.untracked.VPSet(
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedCandidates_packedPFCandidates__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenParticles_prunedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patTriggerObjectStandAlones_slimmedPatTrigger__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedGenParticles_packedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJets__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoVertexs_offlineSlimmedPrimaryVertices__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoCaloClusters_reducedEgamma_reducedESClusters_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEBRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEERecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenJets_slimmedGenJets__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJetsPuppi__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedESRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        )
    ),
    overrideInputFileSplitLevels = cms.untracked.bool(True),
    splitLevel = cms.untracked.int32(0)
)

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:step3_inDQM.root'),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.mix.input.nbPileupEvents.averageNumber = cms.double(35.000000)
process.mix.bunchspace = cms.int32(25)
process.mix.minBunch = cms.int32(-12)
process.mix.maxBunch = cms.int32(3)
#process.mix.input.fileNames = cms.untracked.vstring(['/store/relval/CMSSW_10_5_0_pre1/RelValMinBias_13/GEN-SIM/103X_mcRun2_asymptotic_v3-v1/20000/C54EE4FC-6249-F94E-B762-364EB260AAE0.root', '/store/relval/CMSSW_10_5_0_pre1/RelValMinBias_13/GEN-SIM/103X_mcRun2_asymptotic_v3-v1/20000/D7FA26FE-AD11-874E-A4AC-DAA9F93C908F.root', '/store/relval/CMSSW_10_5_0_pre1/RelValMinBias_13/GEN-SIM/103X_mcRun2_asymptotic_v3-v1/20000/4C2DA47A-296A-8143-B45A-E014D929A02F.root', '/store/relval/CMSSW_10_5_0_pre1/RelValMinBias_13/GEN-SIM/103X_mcRun2_asymptotic_v3-v1/20000/DEC569F7-1C6D-2D4B-B9DB-AC1AD26E992A.root', '/store/relval/CMSSW_10_5_0_pre1/RelValMinBias_13/GEN-SIM/103X_mcRun2_asymptotic_v3-v1/20000/AD61A7A9-F8BA-BB46-BCE2-C01604771D4F.root', '/store/relval/CMSSW_10_5_0_pre1/RelValMinBias_13/GEN-SIM/103X_mcRun2_asymptotic_v3-v1/20000/B6C6E680-3A2A-6D44-AA62-2ED6CD03DD83.root'])
process.mix.input.fileNames = cms.untracked.vstring(['/store/relval/CMSSW_10_5_0_pre1/RelValMinBias_13/GEN-SIM/103X_upgrade2018_realistic_v8-v1/20000/6D675869-EBCE-E14F-8D4F-676D96689255.root', '/store/relval/CMSSW_10_5_0_pre1/RelValMinBias_13/GEN-SIM/103X_upgrade2018_realistic_v8-v1/20000/264C7F72-1A79-C647-A302-B41F2987D828.root', '/store/relval/CMSSW_10_5_0_pre1/RelValMinBias_13/GEN-SIM/103X_upgrade2018_realistic_v8-v1/20000/70BED62C-5DA7-3642-9FA1-DBF9C95609CF.root', '/store/relval/CMSSW_10_5_0_pre1/RelValMinBias_13/GEN-SIM/103X_upgrade2018_realistic_v8-v1/20000/F75D457E-7A3D-6A40-8838-0E2B065984EA.root', '/store/relval/CMSSW_10_5_0_pre1/RelValMinBias_13/GEN-SIM/103X_upgrade2018_realistic_v8-v1/20000/30663F49-0B0C-9E46-B3DC-03318E33D367.root', '/store/relval/CMSSW_10_5_0_pre1/RelValMinBias_13/GEN-SIM/103X_upgrade2018_realistic_v8-v1/20000/9CBFA0CE-63F1-F948-966C-271D056CA5B3.root', '/store/relval/CMSSW_10_5_0_pre1/RelValMinBias_13/GEN-SIM/103X_upgrade2018_realistic_v8-v1/20000/8E775F4B-50AD-BF4B-AF31-696F37D4FB6A.root', '/store/relval/CMSSW_10_5_0_pre1/RelValMinBias_13/GEN-SIM/103X_upgrade2018_realistic_v8-v1/20000/1BBD8C59-3D17-6945-AA58-FC97CB3F72D5.root'])

process.mix.playback = True
process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)
process.RandomNumberGeneratorService.restoreStateLabel=cms.untracked.string("randomEngineStateProducer")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic', '')


############### ADDING NTUPLE in output
###############################
process.TFileService = cms.Service("TFileService", fileName = cms.string("trackingNTuple3.root") )


process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")

from CommonTools.RecoAlgos.trackingParticleRefSelector_cfi import trackingParticleRefSelector as _trackingParticleRefSelector
process.trackingParticlesIntime = _trackingParticleRefSelector.clone(
    signalOnly = False,
    intimeOnly = False,
    chargedOnly = True,
    tip = 1e5,
    lip = 1e5,
    minRapidity = -2.4,
    maxRapidity =  2.4,
    ptMin = 1,
)


process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("SimGeneral.MixingModule.trackingTruthProducerSelection_cfi")

process.trackingParticles.simHitCollections = cms.PSet( )

############## CHANGED THIS BECAUSE ERROR MESSAGE ABOUT MIXING  
#process.mix.playback = cms.untracked.bool(False)
#process.mix.digitizers = cms.PSet(
#     mergedtruth = cms.PSet(process.trackingParticles)
#)
for a in process.aliases: delattr(process, a)

process.trackingPerf = cms.EDAnalyzer('TrackingPerf',
      tracks                   = cms.untracked.InputTag('generalTracks'),
      trackLabel               = cms.InputTag('generalTracks'),
      trackingParticles        = cms.untracked.InputTag('trackingParticlesIntime'),
      trackingParticlesRef     = cms.untracked.bool(True),
      #trackAssociator         = cms.untracked.InputTag('trackingParticleRecoTrackAsssociation'),
      trackAssociator          = cms.untracked.InputTag('quickTrackAssociatorByHits'),
      beamSpot                 = cms.untracked.InputTag('offlineBeamSpot'),
      vertices                 = cms.untracked.InputTag('offlinePrimaryVertices'),
      jetInput                 = cms.InputTag('ak4PFJets'),
      pfJetCollection          = cms.InputTag('ak4PFJets'),
      #pfmetInput               = cms.InputTag('pfMet'),
      genParticles             = cms.InputTag('genParticles'),
      genJetInput              = cms.InputTag("slimmedGenJets"),
      genEventInfoInput        = cms.InputTag("generator"),
      LHEEventProductInput     = cms.InputTag("externalLHEProducer"),
      pfcands                  = cms.InputTag("packedPFCandidates"),
      parametersDefiner        = cms.untracked.string('LhcParametersDefinerForTP'),      
      #electronInput           = cms.untracked.InputTag("slimmedElectrons"),
      #electronInput            = cms.InputTag("patElectronsPFlow"),
      muonInput                = cms.InputTag("slimmedMuons"),
      TTRHBuilder              = cms.string('WithTrackAngle'),
      useCluster               = cms.untracked.bool(False)
)


process.TrackingPerfNtuple = cms.Path( 
	#process.LHCTransport*
	#process.g4SimHits*
	process.mix *
	#process.simHitTPAssocProducer*
	process.tpClusterProducer*
	process.quickTrackAssociatorByHits*
	process.trackingParticleRecoTrackAsssociation*
	process.trackingParticlesIntime*
	process.trackingPerf)


# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.recosim_step = cms.Path(process.recosim)
process.eventinterpretaion_step = cms.Path(process.EIsequence)
process.Flag_trackingFailureFilter = cms.Path(process.goodVertices+process.trackingFailureFilter)
process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)
process.Flag_trkPOGFilters = cms.Path(process.trkPOGFilters)
process.Flag_HcalStripHaloFilter = cms.Path(process.HcalStripHaloFilter)
process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~process.logErrorTooManyClusters)
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)
process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)
process.Flag_eeBadScFilter = cms.Path(process.eeBadScFilter)
process.Flag_METFilters = cms.Path(process.metFilters)
process.Flag_chargedHadronTrackResolutionFilter = cms.Path(process.chargedHadronTrackResolutionFilter)
process.Flag_globalTightHalo2016Filter = cms.Path(process.globalTightHalo2016Filter)
process.Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(process.CSCTightHaloTrkMuUnvetoFilter)
process.Flag_HBHENoiseIsoFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseIsoFilter)
process.Flag_BadChargedCandidateSummer16Filter = cms.Path(process.BadChargedCandidateSummer16Filter)
process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)
process.Flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)
process.Flag_ecalBadCalibFilter = cms.Path(process.ecalBadCalibFilter)
process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseFilter)
process.Flag_trkPOG_toomanystripclus53X = cms.Path(~process.toomanystripclus53X)
process.Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(process.EcalDeadCellBoundaryEnergyFilter)
process.Flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter)
process.Flag_trkPOG_manystripclus53X = cms.Path(~process.manystripclus53X)
process.Flag_BadPFMuonSummer16Filter = cms.Path(process.BadPFMuonSummer16Filter)
process.Flag_muonBadTrackFilter = cms.Path(process.muonBadTrackFilter)
process.Flag_CSCTightHalo2015Filter = cms.Path(process.CSCTightHalo2015Filter)
process.prevalidation_step = cms.Path(process.prevalidation)
process.prevalidation_step1 = cms.Path(process.prevalidationMiniAOD)
process.validation_step = cms.EndPath(process.validation)
process.validation_step1 = cms.EndPath(process.validationMiniAOD)
process.dqmoffline_step = cms.EndPath(process.DQMOffline)
process.dqmoffline_1_step = cms.EndPath(process.DQMOfflineExtraHLT)
process.dqmoffline_2_step = cms.EndPath(process.DQMOfflineMiniAOD)
process.dqmofflineOnPAT_step = cms.EndPath(process.PostDQMOffline)
process.dqmofflineOnPAT_1_step = cms.EndPath(process.PostDQMOffline)
process.dqmofflineOnPAT_2_step = cms.EndPath(process.PostDQMOfflineMiniAOD)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)
process.MINIAODSIMoutput_step = cms.EndPath(process.MINIAODSIMoutput)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.recosim_step,process.eventinterpretaion_step,process.Flag_HBHENoiseFilter,process.Flag_HBHENoiseIsoFilter,process.Flag_CSCTightHaloFilter,process.Flag_CSCTightHaloTrkMuUnvetoFilter,process.Flag_CSCTightHalo2015Filter,process.Flag_globalTightHalo2016Filter,process.Flag_globalSuperTightHalo2016Filter,process.Flag_HcalStripHaloFilter,process.Flag_hcalLaserEventFilter,process.Flag_EcalDeadCellTriggerPrimitiveFilter,process.Flag_EcalDeadCellBoundaryEnergyFilter,process.Flag_ecalBadCalibFilter,process.Flag_goodVertices,process.Flag_eeBadScFilter,process.Flag_ecalLaserCorrFilter,process.Flag_trkPOGFilters,process.Flag_chargedHadronTrackResolutionFilter,process.Flag_muonBadTrackFilter,process.Flag_BadChargedCandidateFilter,process.Flag_BadPFMuonFilter,process.Flag_BadChargedCandidateSummer16Filter,process.Flag_BadPFMuonSummer16Filter,process.Flag_trkPOG_manystripclus53X,process.Flag_trkPOG_toomanystripclus53X,process.Flag_trkPOG_logErrorTooManyClusters,process.Flag_METFilters,process.prevalidation_step,process.prevalidation_step1,process.validation_step,process.validation_step1,process.dqmoffline_step,process.dqmoffline_1_step,process.dqmoffline_2_step,process.dqmofflineOnPAT_step,process.dqmofflineOnPAT_1_step,process.dqmofflineOnPAT_2_step,process.RECOSIMoutput_step,process.MINIAODSIMoutput_step,process.TrackingPerfNtuple,process.DQMoutput_step)
process.schedule.associate(process.patTask)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from SimGeneral.MixingModule.fullMixCustomize_cff
from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn 

#call to customisation function setCrossingFrameOn imported from SimGeneral.MixingModule.fullMixCustomize_cff
process = setCrossingFrameOn(process)

# End of customisation functions
#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC 

#call to customisation function miniAOD_customizeAllMC imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
process = miniAOD_customizeAllMC(process)

process.options.numberOfThreads=cms.untracked.uint32(8)

# End of customisation functions

# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
