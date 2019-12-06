import FWCore.ParameterSet.Config as cms

from RecoTracker.TkSeedingLayers.seedingLayersEDProducer_cfi import *


SiStripLayerTriplet2 = seedingLayersEDProducer.clone()


SiStripLayerTriplet2.layerList = cms.vstring(    
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
    #TID
    #----------

    'TID1_pos+TID2_pos+TID3_pos','TID1_neg+TID2_neg+TID3_neg',#ring 1-2 (matched)
    'TID1_pos+TID2_pos+MTID3_pos','TID1_neg+TID2_neg+MTID3_neg',#ring 3 (mono)
    'TID1_pos+TID2_pos+MTEC1_pos','TID1_neg+TID2_neg+MTEC1_neg',
    'MTID1_pos+MTID2_pos+MTID3_pos','MTID1_neg+MTID2_neg+MTID3_neg',
    'TID1_pos+MTID2_pos+MTID3_pos','TID1_neg+MTID2_neg+MTID3_neg',
    'TID1_pos+TID2_pos+TEC1_pos','TID1_neg+TID2_neg+TEC1_neg',
    'MTID1_pos+MTID2_pos+MTEC1_pos','MTID1_neg+MTID2_neg+MTEC1_neg',
    'TID1_pos+MTID2_pos+MTEC1_pos','TID1_neg+MTID2_neg+MTEC1_neg',
    'TID1_pos+MTID2_pos+TEC1_pos','TID1_neg+MTID2_neg+TEC1_neg',
    'TID1_pos+TID2_pos+TEC2_pos','TID1_neg+TID2_neg+TEC2_neg',
    'TID1_pos+MTID2_pos+MTEC2_pos','TID1_neg+MTID2_neg+MTEC2_neg',
    'TID1_pos+MTID2_pos+TEC2_pos','TID1_neg+MTID2_neg+TEC2_neg',

    #----------
    #TID+TOB
    #----------

    'MTID1_pos+TOB1+TOB2','MTID1_neg+TOB1+TOB2',
    'MTID1_pos+TOB1+TOB2','MTID1_neg+TOB1+TOB2',
    'MTID2_pos+TOB1+TOB2','MTID2_neg+TOB1+TOB2',
    'MTID3_pos+TOB1+TOB2','MTID3_neg+TOB1+TOB2',

    #----------
    #TIB+TID
    #----------

    'TIB1+TIB2+MTID1_pos', 'TIB1+TIB2+MTID1_neg',
    'TIB1+MTID1_pos+MTID2_pos','TIB1+MTID1_neg+MTID2_neg',
    'TIB1+TIB2+TID1_pos','TIB1+TIB2+TID1_neg',
    'MTIB3+MTIB4+MTID1_pos','MTIB3+MTIB4+MTID1_neg',
    'TIB2+MTIB3+MTID1_pos','TIB2+MTIB3+MTID1_neg',
    'MTIB3+MTIB4+MTID1_pos','MTIB2+MTIB4+MTID1_neg',


    #----------
    #TID+TEC RING 1-4
    #----------

    'TID2_pos+TID3_pos+TEC1_pos' ,'TID2_neg+TID3_neg+TEC1_neg',
    'TID2_pos+TID3_pos+MTEC1_pos','TID2_neg+TID3_neg+MTEC1_neg',
    'MTID2_pos+MTID3_pos+TEC1_pos','MTID2_neg+MTID3_neg+TEC1_neg',
    'TID2_pos+TID3_pos+MTEC1_pos','TID2_neg+TID3_neg+MTEC1_neg',
    'TID2_pos+MTID3_pos+TEC1_pos' ,'TID2_neg+MTID3_neg+TEC1_neg',


    #---------
    #TEC RING 1-4
    #---------

    'TEC1_pos+TEC2_pos+TEC3_pos', 'TEC1_neg+TEC2_neg+TEC3_neg',
    'MTEC1_pos+MTEC2_pos+MTEC3_pos', 'MTEC1_neg+MTEC2_neg+MTEC3_neg',
    'MTEC1_pos+TEC2_pos+TEC3_pos', 'MTEC1_neg+TEC2_neg+TEC3_neg',
    'MTEC1_pos+MTEC2_pos+TEC3_pos', 'MTEC1_neg+MTEC2_neg+TEC3_neg',
    'TEC1_pos+MTEC2_pos+TEC3_pos', 'TEC1_neg+MTEC2_neg+TEC3_neg',
    'MTEC1_pos+MTEC2_pos+MTEC3_pos', 'MTEC1_neg+MTEC2_neg+MTEC3_neg',
    'MTEC1_pos+TEC2_pos+MTEC3_pos', 'MTEC1_neg+TEC2_neg+MTEC3_neg',

    'TEC2_pos+TEC3_pos+TEC4_pos', 'TEC2_neg+TEC3_neg+TEC4_neg',
    'MTEC2_pos+MTEC3_pos+MTEC4_pos', 'MTEC2_neg+MTEC3_neg+MTEC4_neg',
    'MTEC2_pos+TEC3_pos+TEC4_pos', 'MTEC2_neg+TEC3_neg+TEC4_neg',
    'MTEC2_pos+MTEC3_pos+TEC4_pos', 'MTEC2_neg+MTEC3_neg+TEC4_neg',
    'TEC2_pos+MTEC3_pos+TEC4_pos', 'TEC2_neg+MTEC3_neg+TEC4_neg',
    'MTEC2_pos+MTEC3_pos+MTEC4_pos', 'MTEC2_neg+MTEC3_neg+MTEC4_neg',
    'MTEC2_pos+TEC3_pos+MTEC4_pos', 'MTEC2_neg+TEC3_neg+MTEC4_neg',

    'TEC3_pos+TEC4_pos+TEC5_pos', 'TEC3_neg+TEC4_neg+TEC5_neg',
    'MTEC3_pos+MTEC4_pos+MTEC5_pos', 'MTEC3_neg+MTEC4_neg+MTEC5_neg',
    'MTEC3_pos+TEC4_pos+TEC5_pos', 'MTEC3_neg+TEC4_neg+TEC5_neg',
    'MTEC3_pos+MTEC4_pos+TEC5_pos', 'MTEC3_neg+MTEC4_neg+TEC5_neg',
    'TEC3_pos+MTEC4_pos+TEC5_pos', 'TEC3_neg+MTEC4_neg+TEC5_neg',
    'MTEC3_pos+MTEC4_pos+MTEC5_pos', 'MTEC3_neg+MTEC4_neg+MTEC5_neg',
    'MTEC3_pos+TEC4_pos+MTEC5_pos', 'MTEC3_neg+TEC4_neg+MTEC5_neg',

    'TEC4_pos+TEC5_pos+TEC6_pos', 'TEC4_neg+TEC5_neg+TEC6_neg',
    'MTEC4_pos+MTEC5_pos+MTEC6_pos', 'MTEC4_neg+MTEC5_neg+MTEC6_neg',
    'MTEC4_pos+TEC5_pos+TEC6_pos', 'MTEC4_neg+TEC5_neg+TEC6_neg',
    'MTEC4_pos+MTEC5_pos+TEC6_pos', 'MTEC4_neg+MTEC5_neg+TEC6_neg',
    'TEC4_pos+MTEC5_pos+TEC6_pos', 'TEC4_neg+MTEC5_neg+TEC6_neg',
    'MTEC4_pos+MTEC5_pos+MTEC6_pos', 'MTEC4_neg+MTEC5_neg+MTEC6_neg',
    'MTEC4_pos+TEC5_pos+MTEC6_pos', 'MTEC4_neg+TEC5_neg+MTEC6_neg'

)



SiStripLayerTriplet2.TOB = cms.PSet(
    matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
    TTRHBuilder = cms.string('WithTrackAngle'),
    clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
         skipClusters   = cms.InputTag('siStripTripletStep2Clusters')
)


SiStripLayerTriplet2.MTOB = cms.PSet(
    TTRHBuilder = cms.string('WithTrackAngle'),
    rphiRecHits    = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
    clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
         skipClusters   = cms.InputTag('siStripTripletStep2Clusters')
)


SiStripLayerTriplet2.TIB = cms.PSet(
         TTRHBuilder    = cms.string('WithTrackAngle'), 
	 clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
         matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
         skipClusters   = cms.InputTag('siStripTripletStep2Clusters')
)

SiStripLayerTriplet2.MTIB = cms.PSet(
         TTRHBuilder    = cms.string('WithTrackAngle'), 
	 clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
         rphiRecHits    = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
         skipClusters   = cms.InputTag('siStripTripletStep2Clusters')
)

SiStripLayerTriplet2.TID = cms.PSet(
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        useRingSlector = cms.bool(True),
        TTRHBuilder = cms.string('WithTrackAngle'), 
	clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
        minRing = cms.int32(1),
        maxRing = cms.int32(2),
         skipClusters   = cms.InputTag('siStripTripletStep2Clusters')

)

SiStripLayerTriplet2.MTID = cms.PSet(
        rphiRecHits    = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        useRingSlector = cms.bool(True),
        TTRHBuilder = cms.string('WithTrackAngle'), 
	clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
        minRing = cms.int32(3),
        maxRing = cms.int32(3),
         skipClusters   = cms.InputTag('siStripTripletStep2Clusters')

)

SiStripLayerTriplet2.TEC = cms.PSet(
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        useRingSlector = cms.bool(True),
        TTRHBuilder = cms.string('WithTrackAngle'), 
	clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
        minRing = cms.int32(5),
        maxRing = cms.int32(5),
         skipClusters   = cms.InputTag('siStripTripletStep2Clusters')

)

SiStripLayerTriplet2.MTEC = cms.PSet(
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        useRingSlector = cms.bool(True),
        TTRHBuilder = cms.string('WithTrackAngle'), 
	clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
        minRing = cms.int32(6),
        maxRing = cms.int32(7),
         skipClusters   = cms.InputTag('siStripTripletStep2Clusters')

)


