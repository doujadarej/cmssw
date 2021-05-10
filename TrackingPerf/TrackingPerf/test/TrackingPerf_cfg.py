import FWCore.ParameterSet.Config as cms

from RecoLocalTracker.Configuration.RecoLocalTracker_cff import *
from SimGeneral.TrackingAnalysis.simHitTPAssociation_cfi import *
from SimTracker.TrackerHitAssociation.tpClusterProducer_cfi import *
from SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi import *

from RecoTracker.TransientTrackingRecHit.TTRHBuilders_cff import *
from RecoLocalTracker.SiPixelRecHits.PixelCPEGeneric_cfi import *
from SimTracker.TrackAssociation.LhcParametersDefinerForTP_cfi import *

# Track Associators
from SimTracker.TrackAssociatorProducers.trackAssociatorByHits_cfi import *



process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

#
#process.load('Configuration.StandardSequences.Reconstruction_cff')
#
#process.load('Configuration.Geometry.GeometryRecoDB_cff')
#process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#
#process.load("RecoTracker.TrackProducer.TrackRefitters_cff") 
#process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = '92X_upgrade2017_realistic_v7'  
#

#process.GlobalTag.globaltag = cms.string('auto')

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
#process.GlobalTag.globaltag = '92X_upgrade2017_realistic_v7' 
from Configuration.AlCa.GlobalTag import GlobalTag 
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic', '')



process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2) )



process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:step3.root'),	 
)


process.TFileService = cms.Service("TFileService", fileName = cms.string("trackingNTuple.root") )



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
#process.RandomNumberGeneratorService.restoreStateLabel=cms.untracked.string("randomEngineStateProducer")


process.load("SimGeneral.TrackingAnalysis.simHitTPAssociation_cfi")

process.load("SimTracker.TrackerHitAssociation.tpClusterProducer_cfi")
process.load("SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi")
#process.load("SimTracker.TrackAssociatorProducers.trackAssociatorByHits_cfi")
process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")
process.load("SimG4Core.Application.g4SimHits_cfi")


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
      #electronInput            = cms.untracked.InputTag("slimmedElectrons"),
      #electronInput            = cms.InputTag("patElectronsPFlow"),
      muonInput                = cms.InputTag("slimmedMuons"),
      TTRHBuilder              = cms.string('WithTrackAngle'),
      useCluster               = cms.untracked.bool(False)
)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("Configuration.StandardSequences.Simulation_cff")


process.load("Geometry.CMSCommonData.cmsExtendedGeometryXML_cfi")
#process.load("SimTransport.HectorProducer.HectorTransport_cfi"



process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
	g4SimHits = cms.PSet(
   		 initialSeed = cms.untracked.uint32(123456789),
    		 engineName = cms.untracked.string('TRandom3')
  	),
  	LHCTransport = cms.PSet(
    		initialSeed = cms.untracked.uint32(321456789),
    		engineName = cms.untracked.string('TRandom3')
  	)

)

process.simHitTPAssocProducer.simHitSrc = cms.VInputTag( 
	cms.InputTag('g4SimHits','TrackerHitsTIBLowTof'),
	cms.InputTag('g4SimHits','TrackerHitsTIBHighTof'),
	cms.InputTag('g4SimHits','TrackerHitsTIDLowTof'),
	cms.InputTag('g4SimHits','TrackerHitsTIDHighTof'),
	cms.InputTag('g4SimHits','TrackerHitsTOBLowTof'),
	cms.InputTag('g4SimHits','TrackerHitsTOBHighTof'),
	cms.InputTag('g4SimHits','TrackerHitsTECLowTof'),
	cms.InputTag('g4SimHits','TrackerHitsTECHighTof'),
	cms.InputTag( 'g4SimHits','TrackerHitsPixelBarrelLowTof'),
	cms.InputTag('g4SimHits','TrackerHitsPixelBarrelHighTof'),
	cms.InputTag('g4SimHits','TrackerHitsPixelEndcapLowTof'),
	cms.InputTag('g4SimHits','TrackerHitsPixelEndcapHighTof') 
	)


#process.p = cms.Path(trackingParticlesIntime*simHitTPAssocProducer*tpClusterProducer*process.trackingPerf)
process.p = cms.Path( 
	#process.LHCTransport*
	#process.g4SimHits*
	process.mix *
	#process.simHitTPAssocProducer*
	process.tpClusterProducer*
	process.quickTrackAssociatorByHits*
	process.trackingParticleRecoTrackAsssociation*
	process.trackingParticlesIntime*
	process.trackingPerf)











