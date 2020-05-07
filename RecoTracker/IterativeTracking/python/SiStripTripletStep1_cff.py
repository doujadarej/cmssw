import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Modifier_tracker_apv_vfp30_2016_cff import tracker_apv_vfp30_2016 as _tracker_apv_vfp30_2016
import RecoTracker.IterativeTracking.iterativeTkConfig as _cfg

#----------------------------------------- NEW CLUSTERS (remove previously used clusters)
siStripTripletStep1Clusters = _cfg.clusterRemoverForIter("SiStripTripletStep1")
for _eraName, _postfix, _era in _cfg.nonDefaultEras():
    _era.toReplaceWith(siStripTripletStep1Clusters, _cfg.clusterRemoverForIter("SiStripTripletStep1", _eraName, _postfix))

#----------------------------------------- SEEDING LAYERS
import RecoTracker.TkSeedingLayers.SiStripLayerTriplet1_cfi
siStripTripletStep1SeedLayers = RecoTracker.TkSeedingLayers.SiStripLayerTriplet1_cfi.SiStripLayerTriplet1.clone()


#----------------------------------------- TrackingRegion
from RecoTracker.TkTrackingRegions.globalTrackingRegion_cfi import globalTrackingRegion as _globalTrackingRegion
#from RecoTracker.TkTrackingRegions.GlobalTrackingRegion_cfi import GlobalTrackingRegion as _globalTrackingRegion
siStripTripletStep1TrackingRegions = _globalTrackingRegion.clone(
   RegionPSet = dict(
     precise = cms.bool(True),
     useMultipleScattering = cms.bool(True),
     originHalfLength = cms.double(55),
     originRadius = cms.double(10),
     ptMin = cms.double(1)
   )
)




#----------------------------------------- Triplet seeding

from RecoPixelVertexing.PixelLowPtUtilities.ClusterShapeHitFilterESProducer_cfi import ClusterShapeHitFilterESProducer as _ClusterShapeHitFilterESProducer
siStripTripletStep1ClusterShapeHitFilter = _ClusterShapeHitFilterESProducer.clone(
    ComponentName = 'siStripTripletStep1ClusterShapeHitFilter',
    doStripShapeCut = cms.bool(False),
    clusterChargeCut = dict(refToPSet_ = 'SiStripClusterChargeCutTight')
)

from RecoTracker.TkHitPairs.hitPairEDProducer_cfi import hitPairEDProducer as _hitPairEDProducer
siStripTripletStep1HitDoublets = _hitPairEDProducer.clone(
    seedingLayers = "siStripTripletStep1SeedLayers",
    trackingRegions = "siStripTripletStep1TrackingRegions",
    maxElement = 500000000,
    produceIntermediateHitDoublets = True,
)

from RecoTracker.TkSeedGenerator.multiHitFromChi2EDProducer_cfi import multiHitFromChi2EDProducer as _multiHitFromChi2EDProducer
siStripTripletStep1HitTriplets = _multiHitFromChi2EDProducer.clone(
    doublets = "siStripTripletStep1HitDoublets",
    extraPhiKDBox = 0.01,
)


from RecoTracker.TkSeedGenerator.seedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer_cff import seedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer as _seedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer
from RecoPixelVertexing.PixelLowPtUtilities.StripSubClusterShapeSeedFilter_cfi import StripSubClusterShapeSeedFilter as _StripSubClusterShapeSeedFilter
siStripTripletStep1Seeds = _seedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer.clone(
    seedingHitSets = "siStripTripletStep1HitTriplets",
    SeedComparitorPSet = dict(
        ComponentName = 'CombinedSeedComparitor',
        mode = cms.string("and"),
        comparitors = cms.VPSet(
            cms.PSet(# FIXME: is this defined in any cfi that could be imported instead of copy-paste?
                ComponentName = cms.string('PixelClusterShapeSeedComparitor'),
                FilterAtHelixStage = cms.bool(True),
                FilterPixelHits = cms.bool(False),
                FilterStripHits = cms.bool(True),
                ClusterShapeHitFilterName = cms.string('siStripTripletStep1ClusterShapeHitFilter'),
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


#_siStripTripletStepSeedComparitorPSet = dict(
#    ComponentName = 'CombinedSeedComparitor',
#    mode = cms.string("and"),
#    comparitors = cms.VPSet(
#        cms.PSet(# FIXME: is this defined in any cfi that could be imported instead of copy-paste?
#            ComponentName = cms.string('PixelClusterShapeSeedComparitor'),
#            FilterAtHelixStage = cms.bool(True),
#            FilterPixelHits = cms.bool(False),
#            FilterStripHits = cms.bool(True),
#            ClusterShapeHitFilterName = cms.string('siStripTripletStepClusterShapeHitFilter'),
#            ClusterShapeCacheSrc = cms.InputTag("siPixelClusterShapeCache") # not really needed here since FilterPixelHits=False
#        ),
#        _StripSubClusterShapeSeedFilter.clone()
#    )
#)

#siStripTripletStepSeeds = _seedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer.clone(#empirically better than 'SeedFromConsecutiveHitsTripletOnlyCreator'
#    seedingHitSets = "siStripTripletStepHitTriplets",
#    SeedComparitorPSet = _siStripTripletStepSeedComparitorPSet,
#)









#----------------------------------------- QUALITY CUTS DURING TRACK BUILDING
import TrackingTools.TrajectoryFiltering.TrajectoryFilter_cff as _TrajectoryFilter_cff
_siStripTripletStep1TrajectoryFilterBase = _TrajectoryFilter_cff.CkfBaseTrajectoryFilter_block.clone(
    maxLostHits = 1,
    minimumNumberOfHits = 4,
    minPt = 1,
)

siStripTripletStep1TrajectoryFilter = _siStripTripletStep1TrajectoryFilterBase.clone(
    seedPairPenalty = 1,
)


siStripTripletStep1TrajectoryFilterInOut = siStripTripletStep1TrajectoryFilter.clone(
    minimumNumberOfHits = 4,
)

import RecoTracker.MeasurementDet.Chi2ChargeMeasurementEstimator_cfi
siStripTripletStep1Chi2Est = RecoTracker.MeasurementDet.Chi2ChargeMeasurementEstimator_cfi.Chi2ChargeMeasurementEstimator.clone(
    ComponentName = cms.string('siStripTripletStep1Chi2Est'),
    nSigma = cms.double(3.0),
    MaxChi2 = cms.double(10.0),
    clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutTight')),
)



#----------------------------------------- TRACK BUILDING
import RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilder_cfi
siStripTripletStep1TrajectoryBuilder = RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilder_cfi.GroupedCkfTrajectoryBuilder.clone(
    MeasurementTrackerName = '',
    trajectoryFilter = cms.PSet(refToPSet_ = cms.string('siStripTripletStep1TrajectoryFilter')),
    inOutTrajectoryFilter = cms.PSet(refToPSet_ = cms.string('siStripTripletStep1TrajectoryFilterInOut')),
    useSameTrajFilter = False,
    minNrOfHitsForRebuild = 4,
    maxCand = 2,
    estimator = cms.string('siStripTripletStep1Chi2Est'),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    # 0.63 GeV is the maximum pT for a charged particle to loop within the 1.1m radius
    # of the outermost Tracker barrel layer (with B=3.8T)
    maxPtForLooperReconstruction = cms.double(0.7) 
    )



#----------------------------------------- MAKING OF TRACK CANDIDATES
import RecoTracker.CkfPattern.CkfTrackCandidates_cfi
siStripTripletStep1TrackCandidates = RecoTracker.CkfPattern.CkfTrackCandidates_cfi.ckfTrackCandidates.clone(
    src = cms.InputTag('siStripTripletStep1Seeds'),
    ### these two parameters are relevant only for the CachingSeedCleanerBySharedInput
    numHitsForSeedCleaner = cms.int32(50),
    onlyPixelHitsForSeedCleaner = cms.bool(True),

    TrajectoryBuilderPSet = cms.PSet(refToPSet_ = cms.string('siStripTripletStep1TrajectoryBuilder')),
    #clustersToSkip = cms.InputTag('siStripTripletStep1Clusters'),
    doSeedingRegionRebuilding = True,
    useHitsSplitting = True,
    cleanTrajectoryAfterInOut = True
)

from TrackingTools.TrajectoryCleaning.TrajectoryCleanerBySharedHits_cfi import trajectoryCleanerBySharedHits
siStripTripletStep1TrajectoryCleanerBySharedHits = trajectoryCleanerBySharedHits.clone(
    ComponentName = cms.string('siStripTripletStep1TrajectoryCleanerBySharedHits'),
    fractionShared = cms.double(0.25),
    allowSharedFirstHit = cms.bool(True)
    )
siStripTripletStep1TrackCandidates.TrajectoryCleaner = 'siStripTripletStep1TrajectoryCleanerBySharedHits'





# ----------------------------------------- TRACK FITTING AND SMOOTHING OPTIONS
import TrackingTools.TrackFitters.RungeKuttaFitters_cff
siStripTripletStep1FitterSmoother = TrackingTools.TrackFitters.RungeKuttaFitters_cff.KFFittingSmootherWithOutliersRejectionAndRK.clone(
    ComponentName = 'siStripTripletStep1FitterSmoother',
    EstimateCut = 30,
    MinNumberOfHits = 8,
    Fitter = cms.string('siStripTripletStep1RKFitter'),
    Smoother = cms.string('siStripTripletStep1RKSmoother')
    )


siStripTripletStep1FitterSmootherForLoopers = siStripTripletStep1FitterSmoother.clone(
    ComponentName = 'siStripTripletStep1FitterSmootherForLoopers',
    Fitter = cms.string('siStripTripletStep1RKFitterForLoopers'),
    Smoother = cms.string('siStripTripletStep1RKSmootherForLoopers')
)

# Also necessary to specify minimum number of hits after final track fit
siStripTripletStep1RKTrajectoryFitter = TrackingTools.TrackFitters.RungeKuttaFitters_cff.RKTrajectoryFitter.clone(
    ComponentName = cms.string('siStripTripletStep1RKFitter'),
    minHits = 8
)


siStripTripletStep1RKTrajectoryFitterForLoopers = siStripTripletStep1RKTrajectoryFitter.clone(
    ComponentName = cms.string('siStripTripletStep1RKFitterForLoopers'),
    Propagator = cms.string('PropagatorWithMaterialForLoopers'),
)

siStripTripletStep1RKTrajectorySmoother = TrackingTools.TrackFitters.RungeKuttaFitters_cff.RKTrajectorySmoother.clone(
    ComponentName = cms.string('siStripTripletStep1RKSmoother'),
    errorRescaling = 10.0,
    minHits = 8
)


siStripTripletStep1RKTrajectorySmootherForLoopers = siStripTripletStep1RKTrajectorySmoother.clone(
    ComponentName = cms.string('siStripTripletStep1RKSmootherForLoopers'),
    Propagator = cms.string('PropagatorWithMaterialForLoopers'),
)

import TrackingTools.TrackFitters.FlexibleKFFittingSmoother_cfi
siStripTripletFlexibleKFFittingSmoother = TrackingTools.TrackFitters.FlexibleKFFittingSmoother_cfi.FlexibleKFFittingSmoother.clone(
    ComponentName = cms.string('siStripTripletFlexibleKFFittingSmoother'),
    standardFitter = cms.string('siStripTripletStep1FitterSmoother'),
    looperFitter = cms.string('siStripTripletStep1FitterSmootherForLoopers'),
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
siStripTripletStep1Tracks = RecoTracker.TrackProducer.TrackProducer_cfi.TrackProducer.clone(
    src = 'siStripTripletStep1TrackCandidates',
    AlgorithmName = cms.string('siStripTripletStep1'),
    #Fitter = 'siStripTripletStep1FitterSmoother',
    Fitter = 'siStripTripletFlexibleKFFittingSmoother',
    )


#---------------------------------------- TRACK SELECTION AND QUALITY FLAG SETTING.

from RecoTracker.FinalTrackSelectors.TrackMVAClassifierPrompt_cfi import *
from RecoTracker.FinalTrackSelectors.TrackMVAClassifierDetached_cfi import *
siStripTripletStep1Classifier1 = TrackMVAClassifierDetached.clone()
siStripTripletStep1Classifier1.src = 'siStripTripletStep1Tracks'
siStripTripletStep1Classifier1.mva.GBRForestLabel = 'MVASelectorIter6_13TeV'
siStripTripletStep1Classifier1.qualityCuts = [-0.6,-0.45,-0.3]
siStripTripletStep1Classifier2 = TrackMVAClassifierPrompt.clone()
siStripTripletStep1Classifier2.src = 'siStripTripletStep1Tracks'
siStripTripletStep1Classifier2.mva.GBRForestLabel = 'MVASelectorIter0_13TeV'
siStripTripletStep1Classifier2.qualityCuts = [0.0,0.0,0.0]
from RecoTracker.FinalTrackSelectors.ClassifierMerger_cfi import *
siStripTripletStep1 = ClassifierMerger.clone()
siStripTripletStep1.inputClassifiers=['siStripTripletStep1Classifier1','siStripTripletStep1Classifier2']
from Configuration.Eras.Modifier_trackingPhase1_cff import trackingPhase1
#commented# from Configuration.Eras.Modifier_trackingPhase1QuadProp_cff import trackingPhase1QuadProp
trackingPhase1.toReplaceWith(siStripTripletStep1, siStripTripletStep1Classifier1.clone(
    mva = dict(GBRForestLabel = 'MVASelectorTobTecStep_Phase1'),
    #commented# GBRForestLabel = 'MVASelectorTobTecStep_Phase1',
    qualityCuts = [-0.6,-0.45,-0.3],
))
#commented# trackingPhase1QuadProp.toReplaceWith(siStripTripletStep1, siStripTripletStep1Classifier1.clone(
#commented#      GBRForestLabel = 'MVASelectorTobTecStep_Phase1',
#commented#      qualityCuts = [-0.6,-0.45,-0.3],
#commented# ))


SiStripTripletStep1Task = cms.Task(siStripTripletStep1Clusters,
                          siStripTripletStep1SeedLayers,
                          siStripTripletStep1TrackingRegions,
                          siStripTripletStep1HitDoublets,
                          siStripTripletStep1HitTriplets,
                          siStripTripletStep1Seeds,
                          #siStripTripletStep1SeedLayersPair,
                          #siStripTripletStep1TrackingRegionsPair,
                          #siStripTripletStep1HitDoubletsPair,
                          #siStripTripletStep1SeedsPair,
                          siStripTripletStep1TrackCandidates,
                          siStripTripletStep1Tracks,
                          siStripTripletStep1Classifier1,siStripTripletStep1Classifier2,
                          siStripTripletStep1)

SiStripTripletStep1 = cms.Sequence(SiStripTripletStep1Task)


