import FWCore.ParameterSet.Config as cms

from RecoTracker.TkSeedingLayers.seedingLayersEDProducer_cfi import *


DisplacedGeneralLayerTriplet = seedingLayersEDProducer.clone()


DisplacedGeneralLayerTriplet.layerList = cms.vstring(    
    #----------
    #TIB
    #----------

    'TIB1+TIB2+MTIB3',
    'TIB1+TIB2+MTIB4',
    'TIB1+MTIB3+MTIB4',
    'TIB2+MTIB3+MTIB4',

    #----------
    #TOB
    #----------
    'TOB1+TOB2+MTOB3',
    'TOB2+MTOB3+MTOB4',
    'MTOB3+MTOB4+MTOB5',
    'MTOB4+MTOB5+MTOB6',

    #----------
    #TIB+TOB
    #----------
    'MTIB4+TOB1+TOB2',
    'MTIB4+TOB2+MTOB3',
    'MTIB3+TOB1+TOB2',

    #----------
    #TID+TOB
    #----------

    'MTID1_pos+TOB1+TOB2','MTID1_neg+TOB1+TOB2',
    'MTID1_pos+TOB1+TOB2','MTID1_neg+TOB1+TOB2',
    'MTID2_pos+TOB1+TOB2','MTID2_neg+TOB1+TOB2',
    'MTID3_pos+TOB1+TOB2','MTID3_neg+TOB1+TOB2',

    #TOB+MTEC
    'TOB1+TOB2+MTEC1_pos','TOB1+TOB2+MTEC1_neg',

    #TID+TEC
    'TID1+TID2+TEC1_pos', 'TID1+TID2+TEC1_neg', 
    'TID2+TID3+TEC1_pos', #'TID2+TID3+TEC1_neg',
    'TID2+MTID3+TEC1_pos', 'TID2+MTID3+TEC1_neg',
    'MTID3+TEC1_pos+MTEC2_pos', 'MTID3+TEC1_neg+MTEC2_neg'

)





DisplacedGeneralLayerTriplet.TOB = cms.PSet(
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
        skipClusters   = cms.InputTag('displacedGeneralStepClusters')
)


DisplacedGeneralLayerTriplet.MTOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        rphiRecHits    = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
        skipClusters   = cms.InputTag('displacedGeneralStepClusters')
)


DisplacedGeneralLayerTriplet.TIB = cms.PSet(
        TTRHBuilder    = cms.string('WithTrackAngle'), 
	clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        skipClusters   = cms.InputTag('displacedGeneralStepClusters')
)

DisplacedGeneralLayerTriplet.MTIB = cms.PSet(
        TTRHBuilder    = cms.string('WithTrackAngle'), 
	clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
        rphiRecHits    = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        skipClusters   = cms.InputTag('displacedGeneralStepClusters')
)

DisplacedGeneralLayerTriplet.TID = cms.PSet(
        matchedRLayerecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        useRingSlector = cms.bool(True),
        TTRHBuilder = cms.string('WithTrackAngle'), 
	clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
        minRing = cms.int32(1),
        maxRing = cms.int32(2),
        skipClusters   = cms.InputTag('displacedGeneralStepClusters')

)

DisplacedGeneralLayerTriplet.MTID = cms.PSet(
        rphiRecHLayerits    = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        useRingSlector = cms.bool(True),
        TTRHBuilder = cms.string('WithTrackAngle'), 
        clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
        minRing = cms.int32(3),
        maxRing = cms.int32(3),
        skipClusters   = cms.InputTag('displacedGeneralStepClusters')

)

DisplacedGeneralLayerTriplet.TEC = cms.PSet(
        matchedRLayerecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        useRingSlector = cms.bool(True),
        TTRHBuilder = cms.string('WithTrackAngle'), 
        clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
        minRing = cms.int32(5),
        maxRing = cms.int32(5),
        skipClusters   = cms.InputTag('displacedGeneralStepClusters')

)

DisplacedGeneralLayerTriplet.MTEC = cms.PSet(
        rphiRecHLayerits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        useRingSlector = cms.bool(True),
        TTRHBuilder = cms.string('WithTrackAngle'), 
        clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
        minRing = cms.int32(6),
        maxRing = cms.int32(7),
        skipClusters   = cms.InputTag('displacedGeneralStepClusters')

)