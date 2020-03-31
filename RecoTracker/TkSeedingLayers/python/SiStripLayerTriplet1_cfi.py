import FWCore.ParameterSet.Config as cms

from RecoTracker.TkSeedingLayers.seedingLayersEDProducer_cfi import *


SiStripLayerTriplet1 = seedingLayersEDProducer.clone()


SiStripLayerTriplet1.layerList = cms.vstring(    
    #----------
    #TIB+TOB
    #----------

    'TIB1+TIB2+MTIB3+MTIB4',
    "TIB1+MTIB3+MTIB4+TOB1", 
    "TIB1+TIB2+MTIB3+TOB1", 
    "TIB2+MTIB3+MTIB4+TOB1", 
    "TIB2+MTIB4+TOB1+TOB2", 
    "TIB2+MTIB3+TOB1+TOB2", 
    "MTIB3+MTIB4+TOB1+TOB2", 
    "MTIB3+TOB1+TOB2+MTOB3", 
    "MTIB3+MTIB4+TOB2+MTOB3", 
    "MTIB4+TOB1+TOB2+MTOB3", 
    "MTIB4+TOB2+MTOB3+MTOB4", 
    "MTIB4+TOB1+MTOB3+MTOB4",

    "TIB2+MTIB3+MTID1+TOB1", 




    #----------
    #TIB+TID
    #----------

    "TIB1+TIB2+MTID1_pos+MTID2_pos", "TIB1+TIB2+MTID1_neg+MTID2_neg",
    "TIB1+TIB2+MTID1_pos+TOB1", "TIB1+TIB2+MTID1_neg+TOB1",
    "TIB1+MTID1_pos+TOB1+TOB2", "TIB1+MTID1_neg+TOB1+TOB2", 
    
    #"TIB1+MTID1_pos+MTID2_pos+TID3_pos", "TIB1+TID1_neg+TID2_neg+TID3_neg",
    #"TIB1+TIB2+TID1_pos+TID3_pos", "TIB1+TIB2+TID1_neg+TID3_neg",
    #"TIB1+TID1_pos+TID2_pos+TID3_pos", "TIB1+TID1_neg+TID2_neg+TID3_neg", 
    #"TIB2+TID1_pos+TID2_pos+TID3_pos", "TIB2+TID1_neg+TID2_neg+TID3_neg"

    #----------
    #TID+TEC
    #----------

    "TID1_pos+TID2_pos+TID3_pos+TEC1_pos", "TID1_neg+TID2_neg+TID3_neg+TEC1_neg",
    "TID1_pos+TID3_pos+TEC1_pos+TEC2_pos", "TID1_neg+TID3_neg+TEC1_neg+TEC2_neg",
    "TID1_pos+TID2_pos+TEC1_pos+TEC2_pos", "TID1_neg+TID2_neg+TEC1_neg+TEC2_neg",

    "TID2_pos+TID3_pos+TEC1_pos+TEC2_pos", "TID2_neg+TID3_neg+TEC1_neg+TEC2_neg",
    "TID2_pos+TEC1_pos+TEC2_pos+TEC3_pos", "TID2_neg+TEC1_neg+TEC2_neg+TEC3_neg",
    "TID2_pos+TID3_pos+TEC2_pos+TEC3_pos", "TID2_neg+TEC1_neg+TEC2_neg+TEC3_neg",

    "TID3_pos+TEC1_pos+TEC2_pos+TEC3_pos", "TID3_neg+TEC1_neg+TEC2_neg+TEC3_neg",

    "TID1_pos+TID2_pos+TID3_pos+MTEC1_pos", "TID1_neg+TID2_neg+TID3_neg+MTEC1_neg",
    "TID1_pos+TID2_pos+MTID3_pos+MTEC1_pos", "TID1_neg+TID2_neg+MTID3_neg+MTEC1_neg",
    "TID2_pos+TID3_pos+MTEC1_pos+MTEC2_pos", "TID2_neg+TID3_neg+MTEC1_neg+MTEC2_neg", 
    "TID3_pos+MTEC1_pos+MTEC2_pos+MTEC3_pos", "TID3_neg+MTEC1_neg+MTEC2_neg+MTEC3_neg",




    "MTID1_pos+MTID2_pos+MTID3_pos+TEC1_pos", "MTID1_neg+MTID2_neg+MTID3_neg+TEC1_neg",
    "MTID1_pos+MTID3_pos+TEC1_pos+TEC2_pos", "MTID1_neg+MTID3_neg+TEC1_neg+TEC2_neg",
    "MTID1_pos+MTID2_pos+TEC1_pos+TEC2_pos", "MTID1_neg+MTID2_neg+TEC1_neg+TEC2_neg",

    "MTID1_pos+MTID2_pos+MTID3_pos+MTEC1_pos", "MTID1_neg+MTID2_neg+MTID3_neg+MTEC1_neg",
    "MTID1_pos+MTID3_pos+MTEC1_pos+MTEC2_pos", "MTID1_neg+MTID3_neg+MTEC1_neg+MTEC2_neg",
    "MTID1_pos+MTID2_pos+MTEC1_pos+MTEC2_pos", "MTID1_neg+MTID2_neg+MTEC1_neg+MTEC2_neg", 

    "MTID3_pos+MTEC1_pos+MTEC2_pos+MTEC3_pos", "MTID3_neg+MTEC1_neg+MTEC2_neg+MTEC3_neg", 
    "MTID3_pos+MTEC2_pos+MTEC3_pos+MTEC4_pos", "MTID3_neg+MTEC2_neg+MTEC3_neg+MTEC4_neg",
    "MTID3_pos+MTEC1_pos+MTEC3_pos+MTEC4_pos", "MTID3_neg+MTEC1_neg+MTEC3_neg+MTEC4_neg",





    #----------
    #TID+TOB
    #----------
    
    "MTID1_pos+TOB1+TOB2+MTOB3", "MTID1_neg+TOB1+TOB2+MTOB3", 
    "MTID2_pos+TOB1+TOB2+MTOB3", "MTID2_neg+TOB1+TOB2+MTOB3",
    "MTID3_pos+TOB1+TOB2+MTOB3", "MTID3_neg+TOB1+TOB2+MTOB3",








)





SiStripLayerTriplet1.TOB = cms.PSet(
    matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
    TTRHBuilder = cms.string('WithTrackAngle'),
    clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
         skipClusters   = cms.InputTag('siStripTripletStep1Clusters')
)


SiStripLayerTriplet1.MTOB = cms.PSet(
    TTRHBuilder = cms.string('WithTrackAngle'),
    rphiRecHits    = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
    clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
         skipClusters   = cms.InputTag('siStripTripletStep1Clusters')
)


SiStripLayerTriplet1.TIB = cms.PSet(
         TTRHBuilder    = cms.string('WithTrackAngle'), 
	 clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
         matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
         skipClusters   = cms.InputTag('siStripTripletStep1Clusters')
)

SiStripLayerTriplet1.MTIB = cms.PSet(
         TTRHBuilder    = cms.string('WithTrackAngle'), 
	 clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
         rphiRecHits    = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
         skipClusters   = cms.InputTag('siStripTripletStep1Clusters')
)

SiStripLayerTriplet1.TID = cms.PSet(
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        useRingSlector = cms.bool(True),
        TTRHBuilder = cms.string('WithTrackAngle'), 
	clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
        minRing = cms.int32(1),
        maxRing = cms.int32(2),
         skipClusters   = cms.InputTag('siStripTripletStep1Clusters')

)

SiStripLayerTriplet1.MTID = cms.PSet(
        rphiRecHits    = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        useRingSlector = cms.bool(True),
        TTRHBuilder = cms.string('WithTrackAngle'), 
	clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
        minRing = cms.int32(3),
        maxRing = cms.int32(3),
         skipClusters   = cms.InputTag('siStripTripletStep1Clusters')

)

SiStripLayerTriplet1.TEC = cms.PSet(
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        useRingSlector = cms.bool(True),
        TTRHBuilder = cms.string('WithTrackAngle'), 
	clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
        minRing = cms.int32(1),
        maxRing = cms.int32(2),
         skipClusters   = cms.InputTag('siStripTripletStep1Clusters')

)

SiStripLayerTriplet1.MTEC = cms.PSet(
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        useRingSlector = cms.bool(True),
        TTRHBuilder = cms.string('WithTrackAngle'), 
	clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
        minRing = cms.int32(3),
        maxRing = cms.int32(4),
         skipClusters   = cms.InputTag('siStripTripletStep1Clusters')

)
