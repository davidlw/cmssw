import FWCore.ParameterSet.Config as cms

################################################################################### 
# pp iterative tracking modified for hiOffline reco (the vertex is the one reconstructed in HI)
################################### 3rd step: low-pT and displaced tracks from pixel triplets

from RecoHI.HiTracking.HITrackingRegionProducer_cfi import *

###################################
from RecoTracker.IterativeTracking.DetachedTripletStep_cff import *

# NEW CLUSTERS (remove previously used clusters)
hiDetachedTripletStepClusters = cms.EDProducer("TrackClusterRemover",
                                                clusterLessSolution= cms.bool(True),
                                                oldClusterRemovalInfo = cms.InputTag("hiPixelPairClusters"),
                                                trajectories = cms.InputTag("hiPixelPairGlobalPrimTracks"),
                                                overrideTrkQuals = cms.InputTag('hiPixelPairStepSelector','hiPixelPairStep'),
                                                TrackQuality = cms.string('highPurity'),
                                                pixelClusters = cms.InputTag("siPixelClusters"),
                                                stripClusters = cms.InputTag("siStripClusters"),
                                                Common = cms.PSet(
    maxChi2 = cms.double(9.0),
    ),
                                                Strip = cms.PSet(
    maxChi2 = cms.double(9.0),
    #Yen-Jie's mod to preserve merged clusters
    maxSize = cms.uint32(2)
    )
                                                )



# SEEDING LAYERS
hiDetachedTripletStepSeedLayers =  RecoTracker.IterativeTracking.DetachedTripletStep_cff.detachedTripletStepSeedLayers.clone(
    ComponentName = 'hiDetachedTripletStepSeedLayers'
    )
hiDetachedTripletStepSeedLayers.BPix.skipClusters = cms.InputTag('hiDetachedTripletStepClusters')
hiDetachedTripletStepSeedLayers.FPix.skipClusters = cms.InputTag('hiDetachedTripletStepClusters')

# seeding
#hiDetachedTripletStepSeeds     = RecoTracker.IterativeTracking.DetachedTripletStep_cff.detachedTripletStepSeeds.clone()
#hiDetachedTripletStepSeeds.ClusterCheckPSet.doClusterCheck                             = False # do not check for max number of clusters pixel or strips
#hiDetachedTripletStepSeeds.OrderedHitsFactoryPSet.SeedingLayers = 'hiDetachedTripletStepSeedLayers'
#from RecoPixelVertexing.PixelLowPtUtilities.ClusterShapeHitFilterESProducer_cfi import *
#hiDetachedTripletStepSeeds.OrderedHitsFactoryPSet.GeneratorPSet.SeedComparitorPSet.ComponentName = 'LowPtClusterShapeSeedComparitor'
#hiDetachedTripletStepSeeds.RegionFactoryPSet.RegionPSet.ptMin = 0.3
#hiDetachedTripletStepSeeds.OrderedHitsFactoryPSet.GeneratorPSet.maxElement = 5000000
#hiDetachedTripletStepSeeds.ClusterCheckPSet.MaxNumberOfPixelClusters = 5000000
#hiDetachedTripletStepSeeds.ClusterCheckPSet.MaxNumberOfCosmicClusters = 50000000
#hiDetachedTripletStepSeeds.RegionFactoryPSet.RegionPSet.originHalfLength = 10.0
#hiDetachedTripletStepSeeds.RegionFactoryPSet.RegionPSet.originRadius = 1.5

from RecoTracker.TkTrackingRegions.GlobalTrackingRegionFromBeamSpot_cfi import RegionPsetFomBeamSpotBlock
hiDetachedTripletStepSeeds = RecoTracker.IterativeTracking.DetachedTripletStep_cff.detachedTripletStepSeeds.clone(
    RegionFactoryPSet = RegionPsetFomBeamSpotBlock.clone(
#    ComponentName = cms.string('GlobalRegionProducerFromBeamSpot'),
#    RegionPSet = RegionPsetFomBeamSpotBlock.RegionPSet.clone(
#    ptMin = 4.0,
#    originRadius = 0.005,
#    nSigmaZ = 4.0
#    )
      ComponentName = cms.string('GlobalTrackingRegionWithVerticesProducer'),
      RegionPSet = cms.PSet(
        precise = cms.bool(True),
        beamSpot = cms.InputTag("offlineBeamSpot"),
        useFixedError = cms.bool(True),
        nSigmaZ = cms.double(4.0),
        sigmaZVertex = cms.double(4.0),
        fixedError = cms.double(1.0),
        VertexCollection = cms.InputTag("hiSelectedVertex"),
        ptMin = cms.double(0.4),
        useFoundVertices = cms.bool(True),
        originRadius = cms.double(1.0)
      )
    )
)

hiDetachedTripletStepSeeds.ClusterCheckPSet.doClusterCheck                             = False # do not check for max number of clusters pixel or strips
hiDetachedTripletStepSeeds.OrderedHitsFactoryPSet.SeedingLayers = 'hiDetachedTripletStepSeedLayers'
from RecoPixelVertexing.PixelLowPtUtilities.ClusterShapeHitFilterESProducer_cfi import *
hiDetachedTripletStepSeeds.SeedComparitorPSet.ComponentName = 'LowPtClusterShapeSeedComparitor'
hiDetachedTripletStepSeeds.OrderedHitsFactoryPSet.GeneratorPSet.maxElement = 5000000
hiDetachedTripletStepSeeds.ClusterCheckPSet.MaxNumberOfPixelClusters = 5000000
hiDetachedTripletStepSeeds.ClusterCheckPSet.MaxNumberOfCosmicClusters = 50000000

# building: feed the new-named seeds
hiDetachedTripletStepTrajectoryFilter = RecoTracker.IterativeTracking.DetachedTripletStep_cff.detachedTripletStepTrajectoryFilter.clone(
    ComponentName    = 'hiDetachedTripletStepTrajectoryFilter'
    )

# TRACK BUILDING
import RecoTracker.CkfPattern.CkfTrajectoryBuilderESProducer_cfi
hiDetachedTripletStepTrajectoryBuilder = RecoTracker.CkfPattern.CkfTrajectoryBuilderESProducer_cfi.CkfTrajectoryBuilder.clone(
    ComponentName = 'hiDetachedTripletStepTrajectoryBuilder',
    MeasurementTrackerName = '',
    trajectoryFilterName = 'hiDetachedTripletStepTrajectoryFilter',
    clustersToSkip = cms.InputTag('hiDetachedTripletStepClusters'),
    maxCand = 3,
    estimator = cms.string('detachedTripletStepChi2Est'),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.7)
    ) 

#hiDetachedTripletStepTrajectoryBuilder = RecoTracker.IterativeTracking.DetachedTripletStep_cff.detachedTripletStepTrajectoryBuilder.clone(
#    ComponentName        = 'hiDetachedTripletStepTrajectoryBuilder',
#    trajectoryFilterName = 'hiDetachedTripletStepTrajectoryFilter',
#    clustersToSkip       = cms.InputTag('hiDetachedTripletStepClusters')
#)

hiDetachedTripletStepTrackCandidates        =  RecoTracker.IterativeTracking.DetachedTripletStep_cff.detachedTripletStepTrackCandidates.clone(
    src               = cms.InputTag('hiDetachedTripletStepSeeds'),
    TrajectoryBuilder = 'hiDetachedTripletStepTrajectoryBuilder',
    maxNSeeds=100000
    )

# fitting: feed new-names
hiDetachedTripletStepTracks                 = RecoTracker.IterativeTracking.DetachedTripletStep_cff.detachedTripletStepTracks.clone(
    src                 = 'hiDetachedTripletStepTrackCandidates',
    #AlgorithmName = cms.string('iter7'),
    AlgorithmName = cms.string('iter3')
    )


# Track selection
import RecoHI.HiTracking.hiMultiTrackSelector_cfi
hiDetachedTripletStepSelector = RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiMultiTrackSelector.clone(
    src='hiDetachedTripletStepTracks',
    trackSelectors= cms.VPSet(
    RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiLooseMTS.clone(
    name = 'hiDetachedTripletStepLoose',
    d0_par2 = [9999.0, 0.0],
    dz_par2 = [9999.0, 0.0],
    applyAdaptedPVCuts = False
    ), #end of pset
    RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiTightMTS.clone(
    name = 'hiDetachedTripletStepTight',
    preFilterName = 'hiDetachedTripletStepLoose',
    d0_par2 = [9999.0, 0.0],
    dz_par2 = [9999.0, 0.0],
    applyAdaptedPVCuts = False
    ),
    RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiHighpurityMTS.clone(
    name = 'hiDetachedTripletStep',
    preFilterName = 'hiDetachedTripletStepTight',
    d0_par2 = [9999.0, 0.0],
    dz_par2 = [9999.0, 0.0],
    applyAdaptedPVCuts = False
    ),
    ) #end of vpset
    ) #end of clone  

hiDetachedTripletStep = cms.Sequence(hiDetachedTripletStepClusters*
                                          hiDetachedTripletStepSeeds*
                                          hiDetachedTripletStepTrackCandidates*
                                          hiDetachedTripletStepTracks*
                                          hiDetachedTripletStepSelector
                                          )
