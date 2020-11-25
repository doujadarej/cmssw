import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Modifier_tracker_apv_vfp30_2016_cff import tracker_apv_vfp30_2016 as _tracker_apv_vfp30_2016
import RecoTracker.IterativeTracking.iterativeTkConfig as _cfg

#----------------------------------------- NEW CLUSTERS (remove previously used clusters)
displacedGeneralStepClusters = _cfg.clusterRemoverForIter("DisplacedGeneralStep")
for _eraName, _postfix, _era in _cfg.nonDefaultEras():
    _era.toReplaceWith(displacedGeneralStepClusters, _cfg.clusterRemoverForIter("DisplacedGeneralStep", _eraName, _postfix))

#----------------------------------------- SEEDING LAYERS
import RecoTracker.TkSeedingLayers.DisplacedGeneralLayerTriplet_cfi
displacedGeneralStepSeedLayers = RecoTracker.TkSeedingLayers.DisplacedGeneralLayerTriplet_cfi.DisplacedGeneralLayerTriplet.clone()


#----------------------------------------- TrackingRegion
from RecoTracker.TkTrackingRegions.globalTrackingRegion_cfi import globalTrackingRegion as _globalTrackingRegion
displacedGeneralStepTrackingRegions = _globalTrackingRegion.clone(
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
displacedGeneralStepClusterShapeHitFilter = _ClusterShapeHitFilterESProducer.clone(
    ComponentName = 'displacedGeneralStepClusterShapeHitFilter',
    doStripShapeCut = cms.bool(False),
    clusterChargeCut = dict(refToPSet_ = 'SiStripClusterChargeCutTight')
)

from RecoTracker.TkHitPairs.hitPairEDProducer_cfi import hitPairEDProducer as _hitPairEDProducer
displacedGeneralStepHitDoublets = _hitPairEDProducer.clone(
    seedingLayers = "displacedGeneralStepSeedLayers",
    trackingRegions = "displacedGeneralStepTrackingRegions",
    maxElement = 500000000,
    produceIntermediateHitDoublets = True,
)

from RecoTracker.TkSeedGenerator.multiHitFromChi2EDProducer_cfi import multiHitFromChi2EDProducer as _multiHitFromChi2EDProducer
displacedGeneralStepHitTriplets = _multiHitFromChi2EDProducer.clone(
    doublets = "displacedGeneralStepHitDoublets",
    extraPhiKDBox = 0.01,
)


from RecoTracker.TkSeedGenerator.seedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer_cff import seedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer as _seedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer
from RecoPixelVertexing.PixelLowPtUtilities.StripSubClusterShapeSeedFilter_cfi import StripSubClusterShapeSeedFilter as _StripSubClusterShapeSeedFilter
displacedGeneralStepSeeds = _seedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer.clone(
    seedingHitSets = "displacedGeneralStepHitTriplets",
    SeedComparitorPSet = dict(
        ComponentName = 'CombinedSeedComparitor',
        mode = cms.string("and"),
        comparitors = cms.VPSet(
            cms.PSet(# FIXME: is this defined in any cfi that could be imported instead of copy-paste?
                ComponentName = cms.string('PixelClusterShapeSeedComparitor'),
                FilterAtHelixStage = cms.bool(True),
                FilterPixelHits = cms.bool(False),
                FilterStripHits = cms.bool(True),
                ClusterShapeHitFilterName = cms.string('displacedGeneralStepClusterShapeHitFilter'),
                ClusterShapeCacheSrc = cms.InputTag("siPixelClusterShapeCache") # not really needed here since FilterPixelHits=False
            ), 
            _StripSubClusterShapeSeedFilter.clone()
        )
    )
)



#----------------------------------------- QUALITY CUTS DURING TRACK BUILDING
import TrackingTools.TrajectoryFiltering.TrajectoryFilter_cff as _TrajectoryFilter_cff
_displacedGeneralStepTrajectoryFilterBase = _TrajectoryFilter_cff.CkfBaseTrajectoryFilter_block.clone(
    maxLostHits = 1,
    minimumNumberOfHits = 4,
    minPt = 1,
)

displacedGeneralStepTrajectoryFilter = _displacedGeneralStepTrajectoryFilterBase.clone(
    seedPairPenalty = 1,
)


displacedGeneralStepTrajectoryFilterInOut = displacedGeneralStepTrajectoryFilter.clone(
    minimumNumberOfHits = 4,
)

import RecoTracker.MeasurementDet.Chi2ChargeMeasurementEstimator_cfi
displacedGeneralStepChi2Est = RecoTracker.MeasurementDet.Chi2ChargeMeasurementEstimator_cfi.Chi2ChargeMeasurementEstimator.clone(
    ComponentName = cms.string('displacedGeneralStepChi2Est'),
    nSigma = cms.double(3.0),
    MaxChi2 = cms.double(10.0),
    clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutTight')),
)



#----------------------------------------- TRACK BUILDING
import RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilder_cfi
displacedGeneralStepTrajectoryBuilder = RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilder_cfi.GroupedCkfTrajectoryBuilder.clone(
    MeasurementTrackerName = '',
    trajectoryFilter = cms.PSet(refToPSet_ = cms.string('displacedGeneralStepTrajectoryFilter')),
    inOutTrajectoryFilter = cms.PSet(refToPSet_ = cms.string('displacedGeneralStepTrajectoryFilterInOut')),
    useSameTrajFilter = False,
    minNrOfHitsForRebuild = 4,
    maxCand = 2,
    estimator = cms.string('displacedGeneralStepChi2Est'),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    # 0.63 GeV is the maximum pT for a charged particle to loop within the 1.1m radius
    # of the outermost Tracker barrel layer (with B=3.8T)
    maxPtForLooperReconstruction = cms.double(0.7) 
    )



#----------------------------------------- MAKING OF TRACK CANDIDATES
import RecoTracker.CkfPattern.CkfTrackCandidates_cfi
displacedGeneralStepTrackCandidates = RecoTracker.CkfPattern.CkfTrackCandidates_cfi.ckfTrackCandidates.clone(
    src = cms.InputTag('displacedGeneralStepSeeds'),
    ### these two parameters are relevant only for the CachingSeedCleanerBySharedInput
    numHitsForSeedCleaner = cms.int32(50),
    onlyPixelHitsForSeedCleaner = cms.bool(True),

    TrajectoryBuilderPSet = cms.PSet(refToPSet_ = cms.string('displacedGeneralStepTrajectoryBuilder')),
    clustersToSkip = cms.InputTag('displacedGeneralStepClusters'),
    doSeedingRegionRebuilding = True,
    useHitsSplitting = True,
    cleanTrajectoryAfterInOut = True
)

from TrackingTools.TrajectoryCleaning.TrajectoryCleanerBySharedHits_cfi import trajectoryCleanerBySharedHits
displacedGeneralStepTrajectoryCleanerBySharedHits = trajectoryCleanerBySharedHits.clone(
    ComponentName = cms.string('displacedGeneralStepTrajectoryCleanerBySharedHits'),
    fractionShared = cms.double(0.25),
    allowSharedFirstHit = cms.bool(True)
    )
displacedGeneralStepTrackCandidates.TrajectoryCleaner = 'displacedGeneralStepTrajectoryCleanerBySharedHits'





# ----------------------------------------- TRACK FITTING AND SMOOTHING OPTIONS
import TrackingTools.TrackFitters.RungeKuttaFitters_cff
displacedGeneralStepFitterSmoother = TrackingTools.TrackFitters.RungeKuttaFitters_cff.KFFittingSmootherWithOutliersRejectionAndRK.clone(
    ComponentName = 'displacedGeneralStepFitterSmoother',
    EstimateCut = 30,
    MinNumberOfHits = 8,
    Fitter = cms.string('displacedGeneralStepRKFitter'),
    Smoother = cms.string('displacedGeneralStepRKSmoother')
    )


displacedGeneralStepFitterSmootherForLoopers = displacedGeneralStepFitterSmoother.clone(
    ComponentName = 'displacedGeneralStepFitterSmootherForLoopers',
    Fitter = cms.string('displacedGeneralStepRKFitterForLoopers'),
    Smoother = cms.string('displacedGeneralStepRKSmootherForLoopers')
)

# Also necessary to specify minimum number of hits after final track fit
displacedGeneralStepRKTrajectoryFitter = TrackingTools.TrackFitters.RungeKuttaFitters_cff.RKTrajectoryFitter.clone(
    ComponentName = cms.string('displacedGeneralStepRKFitter'),
    minHits = 8
)


displacedGeneralStepRKTrajectoryFitterForLoopers = displacedGeneralStepRKTrajectoryFitter.clone(
    ComponentName = cms.string('displacedGeneralStepRKFitterForLoopers'),
    Propagator = cms.string('PropagatorWithMaterialForLoopers'),
)

displacedGeneralStepRKTrajectorySmoother = TrackingTools.TrackFitters.RungeKuttaFitters_cff.RKTrajectorySmoother.clone(
    ComponentName = cms.string('displacedGeneralStepRKSmoother'),
    errorRescaling = 10.0,
    minHits = 8
)


displacedGeneralStepRKTrajectorySmootherForLoopers = displacedGeneralStepRKTrajectorySmoother.clone(
    ComponentName = cms.string('displacedGeneralStepRKSmootherForLoopers'),
    Propagator = cms.string('PropagatorWithMaterialForLoopers'),
)

import TrackingTools.TrackFitters.FlexibleKFFittingSmoother_cfi
generalDisplacedFlexibleKFFittingSmoother = TrackingTools.TrackFitters.FlexibleKFFittingSmoother_cfi.FlexibleKFFittingSmoother.clone(
    ComponentName = cms.string('generalDisplacedFlexibleKFFittingSmoother'),
    standardFitter = cms.string('displacedGeneralStepFitterSmoother'),
    looperFitter = cms.string('displacedGeneralStepFitterSmootherForLoopers'),
)




import RecoTracker.TrackProducer.TrackProducer_cfi
displacedGeneralStepTracks = RecoTracker.TrackProducer.TrackProducer_cfi.TrackProducer.clone(
    src = 'displacedGeneralStepTrackCandidates',
    AlgorithmName = cms.string('displacedGeneralStep'),
    Fitter = 'generalDisplacedFlexibleKFFittingSmoother',
    )


#---------------------------------------- TRACK SELECTION AND QUALITY FLAG SETTING.

from RecoTracker.FinalTrackSelectors.TrackMVAClassifierPrompt_cfi import *
from RecoTracker.FinalTrackSelectors.TrackMVAClassifierDetached_cfi import *
displacedGeneralStepClassifier1 = TrackMVAClassifierDetached.clone()
displacedGeneralStepClassifier1.src = 'displacedGeneralStepTracks'
displacedGeneralStepClassifier1.mva.GBRForestLabel = 'MVASelectorIter6_13TeV'
displacedGeneralStepClassifier1.qualityCuts = [-0.6,-0.45,-0.3]
displacedGeneralStepClassifier2 = TrackMVAClassifierPrompt.clone()
displacedGeneralStepClassifier2.src = 'displacedGeneralStepTracks'
displacedGeneralStepClassifier2.mva.GBRForestLabel = 'MVASelectorIter0_13TeV'
displacedGeneralStepClassifier2.qualityCuts = [0.0,0.0,0.0]
from RecoTracker.FinalTrackSelectors.ClassifierMerger_cfi import *
displacedGeneralStep = ClassifierMerger.clone()
displacedGeneralStep.inputClassifiers=['displacedGeneralStepClassifier1','displacedGeneralStepClassifier2']
from Configuration.Eras.Modifier_trackingPhase1_cff import trackingPhase1
#commented# from Configuration.Eras.Modifier_trackingPhase1QuadProp_cff import trackingPhase1QuadProp
trackingPhase1.toReplaceWith(displacedGeneralStep, displacedGeneralStepClassifier1.clone(
    mva = dict(GBRForestLabel = 'MVASelectorTobTecStep_Phase1'),
    qualityCuts = [-0.6,-0.45,-0.3],
))



DisplacedGeneralStepTask = cms.Task(displacedGeneralStepClusters,
                          displacedGeneralStepSeedLayers,
                          displacedGeneralStepTrackingRegions,
                          displacedGeneralStepHitDoublets,
                          displacedGeneralStepHitTriplets,
                          displacedGeneralStepSeeds,
                          displacedGeneralStepTrackCandidates,
                          displacedGeneralStepTracks,
                          displacedGeneralStepClassifier1,displacedGeneralStepClassifier2,
                          displacedGeneralStep)

DisplacedGeneralStep = cms.Sequence(DisplacedGeneralStepTask)
