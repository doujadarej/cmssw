import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Modifier_tracker_apv_vfp30_2016_cff import tracker_apv_vfp30_2016 as _tracker_apv_vfp30_2016
import RecoTracker.IterativeTracking.iterativeTkConfig as _cfg

#----------------------------------------- NEW CLUSTERS (remove previously used clusters)
siStripTripletStep2Clusters = _cfg.clusterRemoverForIter("SiStripTripletStep2")
for _eraName, _postfix, _era in _cfg.nonDefaultEras():
    _era.toReplaceWith(siStripTripletStep2Clusters, _cfg.clusterRemoverForIter("SiStripTripletStep2", _eraName, _postfix))

#----------------------------------------- SEEDING LAYERS
import RecoTracker.TkSeedingLayers.SiStripLayerTriplet2_cfi
siStripTripletStep2SeedLayers = RecoTracker.TkSeedingLayers.SiStripLayerTriplet2_cfi.SiStripLayerTriplet2.clone()


#----------------------------------------- TrackingRegion
from RecoTracker.TkTrackingRegions.globalTrackingRegion_cfi import globalTrackingRegion as _globalTrackingRegion
#from RecoTracker.TkTrackingRegions.GlobalTrackingRegion_cfi import GlobalTrackingRegion as _globalTrackingRegion
siStripTripletStep2TrackingRegions = _globalTrackingRegion.clone(
   RegionPSet = dict(
     precise = cms.bool(True),
     useMultipleScattering = cms.bool(True),
     originHalfLength = cms.double(60),
     originRadius = cms.double(20),
     ptMin = cms.double(1)
   )
)




#----------------------------------------- Triplet seeding

from RecoPixelVertexing.PixelLowPtUtilities.ClusterShapeHitFilterESProducer_cfi import ClusterShapeHitFilterESProducer as _ClusterShapeHitFilterESProducer
siStripTripletStep2ClusterShapeHitFilter = _ClusterShapeHitFilterESProducer.clone(
    ComponentName = 'siStripTripletStep2ClusterShapeHitFilter',
    doStripShapeCut = cms.bool(False),
    clusterChargeCut = dict(refToPSet_ = 'SiStripClusterChargeCutTight')
)

from RecoTracker.TkHitPairs.hitPairEDProducer_cfi import hitPairEDProducer as _hitPairEDProducer
siStripTripletStep2HitDoublets = _hitPairEDProducer.clone(
    seedingLayers = "siStripTripletStep2SeedLayers",
    trackingRegions = "siStripTripletStep2TrackingRegions",
    maxElement = 500000000,
    produceIntermediateHitDoublets = True,
)

from RecoTracker.TkSeedGenerator.multiHitFromChi2EDProducer_cfi import multiHitFromChi2EDProducer as _multiHitFromChi2EDProducer
siStripTripletStep2HitTriplets = _multiHitFromChi2EDProducer.clone(
    doublets = "siStripTripletStep2HitDoublets",
    extraPhiKDBox = 0.01,
)

from RecoTracker.TkSeedGenerator.seedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer_cff import seedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer as _seedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer
from RecoPixelVertexing.PixelLowPtUtilities.StripSubClusterShapeSeedFilter_cfi import StripSubClusterShapeSeedFilter as _StripSubClusterShapeSeedFilter
siStripTripletStep2Seeds = _seedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer.clone(
    seedingHitSets = "siStripTripletStep2HitTriplets",
    SeedComparitorPSet = dict(
        ComponentName = 'CombinedSeedComparitor',
        mode = cms.string("and"),
        comparitors = cms.VPSet(
            cms.PSet(# FIXME: is this defined in any cfi that could be imported instead of copy-paste?
                ComponentName = cms.string('PixelClusterShapeSeedComparitor'),
                FilterAtHelixStage = cms.bool(True),
                FilterPixelHits = cms.bool(False),
                FilterStripHits = cms.bool(True),
                ClusterShapeHitFilterName = cms.string('siStripTripletStep2ClusterShapeHitFilter'),
                ClusterShapeCacheSrc = cms.InputTag("siPixelClusterShapeCache") # not really needed here since FilterPixelHits=False
            ),
            _StripSubClusterShapeSeedFilter.clone()
        )
    )
)
#from RecoPixelVertexing.PixelTriplets.pixelTripletLargeTipEDProducer_cfi import pixelTripletLargeTipEDProducer as _pixelTripletLargeTipEDProducer
#from RecoPixelVertexing.PixelLowPtUtilities.ClusterShapeHitFilterESProducer_cfi import *
#siStripTripletStepHitTriplets = _pixelTripletHLTEDProducer.clone(
#    doublets = "siStripTripletStepHitDoublets",
#    produceSeedingHitSets = True
#)
#
#from RecoTracker.TkSeedGenerator.seedCreatorFromRegionConsecutiveHitsEDProducer_cff import seedCreatorFromRegionConsecutiveHitsEDProducer as _seedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer
#from RecoPixelVertexing.PixelLowPtUtilities.StripSubClusterShapeSeedFilter_cfi import StripSubClusterShapeSeedFilter as _StripSubClusterShapeSeedFilter


#_siStripTripletStep2SeedComparitorPSet = dict(
#    ComponentName = 'CombinedSeedComparitor',
#    mode = cms.string("and"),
#    comparitors = cms.VPSet(
#        cms.PSet(# FIXME: is this defined in any cfi that could be imported instead of copy-paste?
#            ComponentName = cms.string('PixelClusterShapeSeedComparitor'),
#            FilterAtHelixStage = cms.bool(True),
#            FilterPixelHits = cms.bool(False),
#            FilterStripHits = cms.bool(True),
#            ClusterShapeHitFilterName = cms.string('siStripTripletStep2ClusterShapeHitFilter'),
#            ClusterShapeCacheSrc = cms.InputTag("siPixelClusterShapeCache") # not really needed here since FilterPixelHits=False
#        ),
#        _StripSubClusterShapeSeedFilter.clone()
#    )
#)

#siStripTripletStep2Seeds = _seedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer.clone(#empirically better than 'SeedFromConsecutiveHitsTripletOnlyCreator'
#    seedingHitSets = "siStripTripletStep2HitTriplets",
#    SeedComparitorPSet = _siStripTripletStep2SeedComparitorPSet,
#)









#----------------------------------------- QUALITY CUTS DURING TRACK BUILDING
import TrackingTools.TrajectoryFiltering.TrajectoryFilter_cff as _TrajectoryFilter_cff
_siStripTripletStep2TrajectoryFilterBase = _TrajectoryFilter_cff.CkfBaseTrajectoryFilter_block.clone(
    maxLostHits = 1,
    minimumNumberOfHits = 4,
    minPt = 1,
)

siStripTripletStep2TrajectoryFilter = _siStripTripletStep2TrajectoryFilterBase.clone(
    seedPairPenalty = 1,
)


siStripTripletStep2TrajectoryFilterInOut = siStripTripletStep2TrajectoryFilter.clone(
    minimumNumberOfHits = 4,
)

import RecoTracker.MeasurementDet.Chi2ChargeMeasurementEstimator_cfi
siStripTripletStep2Chi2Est = RecoTracker.MeasurementDet.Chi2ChargeMeasurementEstimator_cfi.Chi2ChargeMeasurementEstimator.clone(
    ComponentName = cms.string('siStripTripletStep2Chi2Est'),
    nSigma = cms.double(3.0),
    MaxChi2 = cms.double(20.0),
    clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutTight')),
)



#----------------------------------------- TRACK BUILDING
import RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilder_cfi
siStripTripletStep2TrajectoryBuilder = RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilder_cfi.GroupedCkfTrajectoryBuilder.clone(
    MeasurementTrackerName = '',
    trajectoryFilter = cms.PSet(refToPSet_ = cms.string('siStripTripletStep2TrajectoryFilter')),
    inOutTrajectoryFilter = cms.PSet(refToPSet_ = cms.string('siStripTripletStep2TrajectoryFilterInOut')),
    useSameTrajFilter = False,
    minNrOfHitsForRebuild = 4,
    maxCand = 2,
    estimator = cms.string('siStripTripletStep2Chi2Est'),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    # 0.63 GeV is the maximum pT for a charged particle to loop within the 1.1m radius
    # of the outermost Tracker barrel layer (with B=3.8T)
    maxPtForLooperReconstruction = cms.double(0.7) 
    )



#----------------------------------------- MAKING OF TRACK CANDIDATES
import RecoTracker.CkfPattern.CkfTrackCandidates_cfi
siStripTripletStep2TrackCandidates = RecoTracker.CkfPattern.CkfTrackCandidates_cfi.ckfTrackCandidates.clone(
    src = cms.InputTag('siStripTripletStep2Seeds'),
    ### these two parameters are relevant only for the CachingSeedCleanerBySharedInput
    numHitsForSeedCleaner = cms.int32(50),
    onlyPixelHitsForSeedCleaner = cms.bool(True),

    TrajectoryBuilderPSet = cms.PSet(refToPSet_ = cms.string('siStripTripletStep2TrajectoryBuilder')),
    #clustersToSkip = cms.InputTag('siStripTripletStep2Clusters'),
    doSeedingRegionRebuilding = True,
    useHitsSplitting = True,
    cleanTrajectoryAfterInOut = True
)

from TrackingTools.TrajectoryCleaning.TrajectoryCleanerBySharedHits_cfi import trajectoryCleanerBySharedHits
siStripTripletStep2TrajectoryCleanerBySharedHits = trajectoryCleanerBySharedHits.clone(
    ComponentName = cms.string('siStripTripletStep2TrajectoryCleanerBySharedHits'),
    fractionShared = cms.double(0.25),
    allowSharedFirstHit = cms.bool(True)
    )
siStripTripletStep2TrackCandidates.TrajectoryCleaner = 'siStripTripletStep2TrajectoryCleanerBySharedHits'





# ----------------------------------------- TRACK FITTING AND SMOOTHING OPTIONS
import TrackingTools.TrackFitters.RungeKuttaFitters_cff
siStripTripletStep2FitterSmoother = TrackingTools.TrackFitters.RungeKuttaFitters_cff.KFFittingSmootherWithOutliersRejectionAndRK.clone(
    ComponentName = 'siStripTripletStep2FitterSmoother',
    EstimateCut = 30,
    MinNumberOfHits = 4,
    Fitter = cms.string('siStripTripletStep2RKFitter'),
    Smoother = cms.string('siStripTripletStep2RKSmoother')
    )


siStripTripletStep2FitterSmootherForLoopers = siStripTripletStep2FitterSmoother.clone(
    ComponentName = 'siStripTripletStep2FitterSmootherForLoopers',
    Fitter = cms.string('siStripTripletStep2RKFitterForLoopers'),
    Smoother = cms.string('siStripTripletStep2RKSmootherForLoopers')
)

# Also necessary to specify minimum number of hits after final track fit
siStripTripletStep2RKTrajectoryFitter = TrackingTools.TrackFitters.RungeKuttaFitters_cff.RKTrajectoryFitter.clone(
    ComponentName = cms.string('siStripTripletStep2RKFitter'),
    minHits = 4
)


siStripTripletStep2RKTrajectoryFitterForLoopers = siStripTripletStep2RKTrajectoryFitter.clone(
    ComponentName = cms.string('siStripTripletStep2RKFitterForLoopers'),
    Propagator = cms.string('PropagatorWithMaterialForLoopers'),
)

siStripTripletStep2RKTrajectorySmoother = TrackingTools.TrackFitters.RungeKuttaFitters_cff.RKTrajectorySmoother.clone(
    ComponentName = cms.string('siStripTripletStep2RKSmoother'),
    errorRescaling = 10.0,
    minHits = 4
)


siStripTripletStep2RKTrajectorySmootherForLoopers = siStripTripletStep2RKTrajectorySmoother.clone(
    ComponentName = cms.string('siStripTripletStep2RKSmootherForLoopers'),
    Propagator = cms.string('PropagatorWithMaterialForLoopers'),
)

import TrackingTools.TrackFitters.FlexibleKFFittingSmoother_cfi
siStripTripletFlexibleKFFittingSmoother = TrackingTools.TrackFitters.FlexibleKFFittingSmoother_cfi.FlexibleKFFittingSmoother.clone(
    ComponentName = cms.string('siStripTripletFlexibleKFFittingSmoother'),
    standardFitter = cms.string('siStripTripletStep2FitterSmoother'),
    looperFitter = cms.string('siStripTripletStep2FitterSmootherForLoopers'),
)


#----------------------------------------- TRACK FITTING
#import RecoTracker.TrackProducer.TrackProducer_cfi
#siStripTripletStepTracks = RecoTracker.TrackProducer.TrackProducer_cfi.TrackProducer.clone(
#    src = 'siStripTripletStepTrackCandidates',
#    AlgorithmName = cms.string('siStripTripletStep'),
#    Fitter = cms.string('FlexibleKFFittingSmoother')
#    )
#
#from TrackingTools.TrajectoryCleaning.TrajectoryCleanerBySharedHits_cfi import trajectoryCleanerBySharedHits
#siStripTripletStepTrajectoryCleanerBySharedHits = trajectoryCleanerBySharedHits.clone(
#        ComponentName = cms.string('siStripTripletStepTrajectoryCleanerBySharedHits'),
#            fractionShared = cms.double(0.16),
#            allowSharedFirstHit = cms.bool(True)
#            )
#siStripTripletStepTrackCandidates.TrajectoryCleaner = 'siStripTripletStepTrajectoryCleanerBySharedHits'

import RecoTracker.TrackProducer.TrackProducer_cfi
siStripTripletStep2Tracks = RecoTracker.TrackProducer.TrackProducer_cfi.TrackProducer.clone(
    src = 'siStripTripletStep2TrackCandidates',
    AlgorithmName = cms.string('siStripTripletStep2'),
    #Fitter = 'siStripTripletStep2FitterSmoother',
    Fitter = 'siStripTripletFlexibleKFFittingSmoother',
    )



#----------------------------------------- TRACK SELECTION AND QUALITY FLAG SETTING.
from RecoTracker.FinalTrackSelectors.TrackMVAClassifierPrompt_cfi import *
from RecoTracker.FinalTrackSelectors.TrackMVAClassifierDetached_cfi import *
siStripTripletStep2Classifier1 = TrackMVAClassifierDetached.clone()
siStripTripletStep2Classifier1.src = 'siStripTripletStep2Tracks'
siStripTripletStep2Classifier1.mva.GBRForestLabel = 'MVASelectorIter6_13TeV'
siStripTripletStep2Classifier1.qualityCuts = [-0.6,-0.45,-0.3]
siStripTripletStep2Classifier2 = TrackMVAClassifierPrompt.clone()
siStripTripletStep2Classifier2.src = 'siStripTripletStep2Tracks'
siStripTripletStep2Classifier2.mva.GBRForestLabel = 'MVASelectorIter0_13TeV'
siStripTripletStep2Classifier2.qualityCuts = [0.0,0.0,0.0]

from RecoTracker.FinalTrackSelectors.ClassifierMerger_cfi import *
siStripTripletStep2 = ClassifierMerger.clone()
siStripTripletStep2.inputClassifiers=['siStripTripletStep2Classifier1','siStripTripletStep2Classifier2']

from Configuration.Eras.Modifier_trackingPhase1_cff import trackingPhase1
#commented#from Configuration.Eras.Modifier_trackingPhase1QuadProp_cff import trackingPhase1QuadProp
trackingPhase1.toReplaceWith(siStripTripletStep2, siStripTripletStep2Classifier1.clone(
     mva = dict(GBRForestLabel = 'MVASelectorTobTecStep_Phase1'),
     #commented# GBRForestLabel = 'MVASelectorTobTecStep_Phase1',
     qualityCuts = [-0.6,-0.45,-0.3],
))
#commented# trackingPhase1QuadProp.toReplaceWith(siStripTripletStep2, siStripTripletStep2Classifier1.clone(
#commented#     GBRForestLabel = 'MVASelectorTobTecStep_Phase1',
#commented#      qualityCuts = [-0.6,-0.45,-0.3],
#commented# ))







SiStripTripletStep2Task = cms.Task(siStripTripletStep2Clusters,
                          siStripTripletStep2SeedLayers,
                          siStripTripletStep2TrackingRegions,
                          siStripTripletStep2HitDoublets,
                          siStripTripletStep2HitTriplets,
                          siStripTripletStep2Seeds,
                          #siStripTripletStep2SeedLayersPair,
                          #siStripTripletStep2TrackingRegionsPair,
                          #siStripTripletStep2HitDoubletsPair,
                          #siStripTripletStep2SeedsPair,
                          siStripTripletStep2TrackCandidates,
                          siStripTripletStep2Tracks,
                          siStripTripletStep2Classifier1,siStripTripletStep2Classifier2,
                          siStripTripletStep2)

SiStripTripletStep2 = cms.Sequence(SiStripTripletStep2Task)