
import FWCore.ParameterSet.Config as cms

from RecoTracker.TkSeedingLayers.seedingLayersEDProducer_cfi import *


SiStripLayerTriplets = seedingLayersEDProducer.clone()


SiStripLayerTriplets.layerList = cms.vstring(    
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

    #----------
    #TIB+TID
    #----------
    
    'TIB1+TIB2+MTID1_pos',
    'TIB1+TIB2+MTID1_neg',


    #----------
    #TID+TEC RING 1-3
    #----------

    'TID2_pos+TID3_pos+TEC1_pos' ,
    'TID2_neg+TID3_neg+TEC1_neg',#ring 1-2 (matched)
    'TID2_pos+TID3_pos+MTEC1_pos',
    'TID2_neg+TID3_neg+MTEC1_neg',#ring 3 (mono)

    #TEC RING 1-3

    'TEC1_pos+TEC2_pos+TEC3_pos', 'TEC1_neg+TEC2_neg+TEC3_neg',
    'TEC1_pos+TEC2_pos+MTEC3_pos','TEC1_neg+TEC2_neg+MTEC3_neg',
    'TEC1_pos+TEC2_pos+TEC4_pos', 'TEC1_neg+TEC2_neg+TEC4_neg',
    'TEC1_pos+TEC2_pos+MTEC4_pos','TEC1_neg+TEC2_neg+MTEC4_neg',
    'TEC2_pos+TEC3_pos+TEC4_pos', 'TEC2_neg+TEC3_neg+TEC4_neg',
    'TEC2_pos+TEC3_pos+MTEC4_pos','TEC2_neg+TEC3_neg+MTEC4_neg',
    'TEC2_pos+TEC3_pos+TEC5_pos', 'TEC2_neg+TEC3_neg+TEC5_neg',
    'TEC2_pos+TEC3_pos+TEC6_pos', 'TEC2_neg+TEC3_neg+TEC6_neg',
    'TEC3_pos+TEC4_pos+TEC5_pos', 'TEC3_neg+TEC4_neg+TEC5_neg',
    'TEC3_pos+TEC4_pos+MTEC5_pos','TEC3_neg+TEC4_neg+MTEC5_neg',
    'TEC3_pos+TEC5_pos+TEC6_pos', 'TEC3_neg+TEC5_neg+TEC6_neg',
    'TEC4_pos+TEC5_pos+TEC6_pos', 'TEC4_neg+TEC5_neg+TEC6_neg', 
    'TEC3_pos+TEC5_pos+TEC6_pos', 'TEC3_neg+TEC5_neg+TEC6_neg',
    'TEC4_pos+TEC5_pos+TEC6_pos', 'TEC4_neg+TEC5_neg+TEC6_neg',  
    'TEC3_pos+MTEC5_pos+MTEC6_pos', 'TEC3_neg+MTEC5_neg+MTEC6_neg',
    'TEC4_pos+MTEC5_pos+MTEC6_pos', 'TEC4_neg+MTEC5_neg+MTEC6_neg'   

)



SiStripLayerTriplets.TOB = cms.PSet(
    matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
    TTRHBuilder = cms.string('WithTrackAngle'),
    clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone'))
)


SiStripLayerTriplets.MTOB = cms.PSet(
    TTRHBuilder = cms.string('WithTrackAngle'),
    rphiRecHits    = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
    clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone'))
)


SiStripLayerTriplets.TIB = cms.PSet(
         TTRHBuilder    = cms.string('WithTrackAngle'), 
	 clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
         matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
         #skipClusters   = cms.InputTag('pixelLessStepClusters')
)

SiStripLayerTriplets.MTIB = cms.PSet(
         TTRHBuilder    = cms.string('WithTrackAngle'), 
	 clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
         #skipClusters   = cms.InputTag('pixelLessStepClusters'),
         rphiRecHits    = cms.InputTag("siStripMatchedRecHits","rphiRecHit")
)

SiStripLayerTriplets.TID = cms.PSet(
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        #skipClusters = cms.InputTag('pixelLessStepClusters'),
        useRingSlector = cms.bool(True),
        TTRHBuilder = cms.string('WithTrackAngle'), 
	clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
        minRing = cms.int32(1),
        maxRing = cms.int32(2)

)

SiStripLayerTriplets.MTID = cms.PSet(
        rphiRecHits    = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        #skipClusters = cms.InputTag('pixelLessStepClusters'),
        useRingSlector = cms.bool(True),
        TTRHBuilder = cms.string('WithTrackAngle'), 
	clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
        minRing = cms.int32(3),
        maxRing = cms.int32(3)

)

SiStripLayerTriplets.TEC = cms.PSet(
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        #skipClusters = cms.InputTag('pixelLessStepClusters'),
        useRingSlector = cms.bool(True),
        TTRHBuilder = cms.string('WithTrackAngle'), 
	clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
        minRing = cms.int32(1),
        maxRing = cms.int32(2)

)

SiStripLayerTriplets.MTEC = cms.PSet(
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        #skipClusters = cms.InputTag('pixelLessStepClusters'),
        useRingSlector = cms.bool(True),
        TTRHBuilder = cms.string('WithTrackAngle'), 
	clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone')),
        minRing = cms.int32(3),
        maxRing = cms.int32(3)

)






#import FWCore.ParameterSet.Config as cms
#
#from RecoTracker.TkSeedingLayers.seedingLayersEDProducer_cfi import *
#
#SiStripLayerTriplets = seedingLayersEDProducer.clone()
#SiStripLayerTriplets.layerList = cms.vstring(
#	#for barrel SiStrip seeding
#	'TIB1+TIB2+TIB3',
#	'TIB2+TIB3+TIB4',
#	'TIB3+TIB4+TOB1',
#	'TIB4+TOB1+TOB2',
#	'TOB1+TOB2+TOB3',
#	'TOB2+TOB3+TOB4',
#	'TOB3+TOB4+TOB5',
#	'TOB4+TOB5+TOB6',
#	
#	#for endcap SiStrip seeding
#	'TID1_pos+TID2_pos+TID3_pos',
#	'TID2_pos+TID3_pos+TEC1_pos',
#	'TID3_pos+TEC1_pos+TEC2_pos',
#	'TEC1_pos+TEC2_pos+TEC3_pos',
#	'TEC2_pos+TEC3_pos+TEC4_pos',
#	'TEC4_pos+TEC5_pos+TEC6_pos',
#	'TEC6_pos+TEC7_pos+TEC8_pos',
#	'TEC7_pos+TEC8_pos+TEC9_pos',
#	
#	'TID1_neg+TID2_neg+TID3_neg',
#	'TID2_neg+TID3_neg+TEC1_neg',
#	'TID3_neg+TEC1_neg+TEC2_neg',
#	'TEC1_neg+TEC2_neg+TEC3_neg',
#	'TEC2_neg+TEC3_neg+TEC4_neg',
#	'TEC4_neg+TEC5_neg+TEC6_neg',
#	'TEC6_neg+TEC7_neg+TEC8_neg',
#	'TEC7_neg+TEC8_neg+TEC9_neg',
#	
#	
#	#mixed barel and endcap SiStrip seeding
#	
#	
#	'TIB1+TID1_pos+TID2_pos',
#	'TIB2+TID1_pos+TID2_pos',
#	'TIB3+TID1_pos+TID2_pos',
#	'TIB4+TID1_pos+TID2_pos',
#	
#	'TID2_pos+TID3_pos+TEC1_pos',
#	'TID3_pos+TEC1_pos+TEC2_pos',
#	
#	'TOB1+TEC1_pos+TEC2_pos',
#	'TOB2+TEC1_pos+TEC2_pos',
#	'TOB3+TEC1_pos+TEC2_pos',
#	'TOB4+TEC1_pos+TEC2_pos',
#	'TOB5+TEC1_pos+TEC2_pos',
#	'TOB6+TEC1_pos+TEC2_pos',
#	
#	'TIB1+TID1_neg+TID2_neg',
#	'TIB2+TID1_neg+TID2_neg',
#	'TIB3+TID1_neg+TID2_neg',
#	'TIB4+TID1_neg+TID2_neg',
#	
#	'TID2_neg+TID3_neg+TEC1_neg',
#	'TID3_neg+TEC1_neg+TEC2_neg',
#	
#	'TOB1+TEC1_neg+TEC2_neg',
#	'TOB2+TEC1_neg+TEC2_neg',
#	'TOB3+TEC1_neg+TEC2_neg',
#	'TOB4+TEC1_neg+TEC2_neg',
#	'TOB5+TEC1_neg+TEC2_neg',
#	'TOB6+TEC1_neg+TEC2_neg'
#)
#
#
#
#
#SiStripLayerTriplets.TOB = cms.PSet(
#    matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
#    TTRHBuilder = cms.string('WithTrackAngle')
#    ,clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone'))
#)
#
#SiStripLayerTriplets.TIB = cms.PSet(
#        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
#        TTRHBuilder = cms.string('WithTrackAngle')
#        ,clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone'))
#)
#
#SiStripLayerTriplets.TID = cms.PSet(
#    matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
#    useRingSlector = cms.bool(True),
#    TTRHBuilder = cms.string('WithTrackAngle'),
#    minRing = cms.int32(1),
#    maxRing = cms.int32(2)
#   ,clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone'))
#)
#
#SiStripLayerTriplets.TEC = cms.PSet(
#    matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
#    useRingSlector = cms.bool(True),
#    TTRHBuilder = cms.string('WithTrackAngle'),
#    minRing = cms.int32(1),
#    maxRing = cms.int32(2)
#   ,clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutNone'))
#)
#
