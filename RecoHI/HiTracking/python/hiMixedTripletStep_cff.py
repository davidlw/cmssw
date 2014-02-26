import FWCore.ParameterSet.Config as cms

################################################################################### 
# pp iterative tracking modified for hiOffline reco (the vertex is the one reconstructed in HI)
################################### 4th step: large impact parameter tracking using mixed-triplet seeding

from RecoHI.HiTracking.HITrackingRegionProducer_cfi import *

###################################
from RecoTracker.IterativeTracking.MixedTripletStep_cff import *

# NEW CLUSTERS (remove previously used clusters)
hiMixedTripletStepClusters = cms.EDProducer("TrackClusterRemover",
                                                clusterLessSolution= cms.bool(True),
                                                oldClusterRemovalInfo = cms.InputTag("hiDetachedTripletStepClusters"),
                                                trajectories = cms.InputTag("hiDetachedTripletStepTracks"),
                                                overrideTrkQuals = cms.InputTag('hiDetachedTripletStepSelector','hiDetachedTripletStep'),
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



# SEEDING LAYERS A
hiMixedTripletStepSeedLayersA =  RecoTracker.IterativeTracking.MixedTripletStep_cff.mixedTripletStepSeedLayersA.clone(
    ComponentName = 'hiMixedTripletStepSeedLayersA'
    )
hiMixedTripletStepSeedLayersA.BPix.skipClusters = cms.InputTag('hiMixedTripletStepClusters')
hiMixedTripletStepSeedLayersA.FPix.skipClusters = cms.InputTag('hiMixedTripletStepClusters')
hiMixedTripletStepSeedLayersA.TEC.skipClusters  = cms.InputTag('hiMixedTripletStepClusters')
hiMixedTripletStepSeedLayersA.layerList = cms.vstring('BPix1+BPix2+BPix3',
                                                           'BPix1+BPix2+FPix1_pos', 'BPix1+BPix2+FPix1_neg',
                                                           'BPix1+FPix1_pos+FPix2_pos', 'BPix1+FPix1_neg+FPix2_neg',
                                                           'BPix2+FPix1_pos+FPix2_pos', 'BPix2+FPix1_neg+FPix2_neg',)
#                                                           'FPix1_pos+FPix2_pos+TEC1_pos', 'FPix1_neg+FPix2_neg+TEC1_neg',)

# SEEDS A
from RecoTracker.TkTrackingRegions.GlobalTrackingRegionFromBeamSpot_cfi import RegionPsetFomBeamSpotBlock
hiMixedTripletStepSeedsA = RecoTracker.IterativeTracking.MixedTripletStep_cff.mixedTripletStepSeedsA.clone(
    RegionFactoryPSet = RegionPsetFomBeamSpotBlock.clone(
      ComponentName = cms.string('GlobalTrackingRegionWithVerticesProducer'),
      RegionPSet = cms.PSet(
        precise = cms.bool(True),
        beamSpot = cms.InputTag("offlineBeamSpot"),
        useFixedError = cms.bool(True),
        nSigmaZ = cms.double(4.0),
        sigmaZVertex = cms.double(4.0),
        fixedError = cms.double(0.5),
        VertexCollection = cms.InputTag("hiSelectedVertex"),
        ptMin = cms.double(0.5),
        useFoundVertices = cms.bool(True),
        originRadius = cms.double(1.)
      )
    )
)

hiMixedTripletStepSeedsA.ClusterCheckPSet.doClusterCheck                             = False # do not check for max number of clusters pixel or strips
hiMixedTripletStepSeedsA.OrderedHitsFactoryPSet.SeedingLayers = 'hiMixedTripletStepSeedLayersA'
from RecoPixelVertexing.PixelLowPtUtilities.ClusterShapeHitFilterESProducer_cfi import *
#hiMixedTripletStepSeedsA.OrderedHitsFactoryPSet.GeneratorPSet.SeedComparitorPSet.ComponentName = 'LowPtClusterShapeSeedComparitor'
hiMixedTripletStepSeedsA.OrderedHitsFactoryPSet.GeneratorPSet.maxElement = 5000000
hiMixedTripletStepSeedsA.ClusterCheckPSet.MaxNumberOfPixelClusters = 5000000
hiMixedTripletStepSeedsA.ClusterCheckPSet.MaxNumberOfCosmicClusters = 50000000

# SEEDING LAYERS B
hiMixedTripletStepSeedLayersB =  RecoTracker.IterativeTracking.MixedTripletStep_cff.mixedTripletStepSeedLayersB.clone(
    ComponentName = 'hiMixedTripletStepSeedLayersB',
    )
hiMixedTripletStepSeedLayersB.BPix.skipClusters = cms.InputTag('hiMixedTripletStepClusters')
hiMixedTripletStepSeedLayersB.TIB.skipClusters  = cms.InputTag('hiMixedTripletStepClusters')
#hiMixedTripletStepSeedLayersB.layerList = cms.vstring('BPix2+BPix3+TIB1','BPix2+BPix3+TIB2')
hiMixedTripletStepSeedLayersB.layerList = cms.vstring('BPix2+BPix3+TIB1')

hiMixedTripletStepSeedsB = RecoTracker.IterativeTracking.MixedTripletStep_cff.mixedTripletStepSeedsB.clone(
    RegionFactoryPSet = RegionPsetFomBeamSpotBlock.clone(
      ComponentName = cms.string('GlobalTrackingRegionWithVerticesProducer'),
      RegionPSet = cms.PSet(
        precise = cms.bool(True),
        beamSpot = cms.InputTag("offlineBeamSpot"),
        useFixedError = cms.bool(True),
        nSigmaZ = cms.double(4.0),
        sigmaZVertex = cms.double(4.0),
        fixedError = cms.double(0.5),
        VertexCollection = cms.InputTag("hiSelectedVertex"),
        ptMin = cms.double(0.7),
        useFoundVertices = cms.bool(True),
        originRadius = cms.double(1.0)
      )
    )
) 

hiMixedTripletStepSeedsB.ClusterCheckPSet.doClusterCheck                             = False # do not check for max number of clusters pixel or strips
hiMixedTripletStepSeedsB.OrderedHitsFactoryPSet.SeedingLayers = 'hiMixedTripletStepSeedLayersB'
#hiMixedTripletStepSeedsB.OrderedHitsFactoryPSet.GeneratorPSet.SeedComparitorPSet.ComponentName = 'LowPtClusterShapeSeedComparitor'
hiMixedTripletStepSeedsB.OrderedHitsFactoryPSet.GeneratorPSet.maxElement = 5000000
hiMixedTripletStepSeedsB.ClusterCheckPSet.MaxNumberOfPixelClusters = 5000000
hiMixedTripletStepSeedsB.ClusterCheckPSet.MaxNumberOfCosmicClusters = 50000000

# combine seeds
hiMixedTripletStepSeeds = RecoTracker.IterativeTracking.MixedTripletStep_cff.mixedTripletStepSeeds.clone(
    seedCollections = cms.VInputTag(
        cms.InputTag('hiMixedTripletStepSeedsA'),
        cms.InputTag('hiMixedTripletStepSeedsB'),
        )
    )

# track building
hiMixedTripletStepTrajectoryFilter = RecoTracker.IterativeTracking.MixedTripletStep_cff.mixedTripletStepTrajectoryFilter.clone(
    ComponentName        = 'hiMixedTripletStepTrajectoryFilter'
   )

hiMixedTripletStepTrajectoryBuilder = RecoTracker.IterativeTracking.MixedTripletStep_cff.mixedTripletStepTrajectoryBuilder.clone(
    ComponentName        = 'hiMixedTripletStepTrajectoryBuilder',
    trajectoryFilterName = 'hiMixedTripletStepTrajectoryFilter',
    clustersToSkip       = cms.InputTag('hiMixedTripletStepClusters'),
)

hiMixedTripletStepTrackCandidates        =  RecoTracker.IterativeTracking.MixedTripletStep_cff.mixedTripletStepTrackCandidates.clone(
    src               = cms.InputTag('hiMixedTripletStepSeeds'),
    TrajectoryBuilder = 'hiMixedTripletStepTrajectoryBuilder',
    maxNSeeds = 100000
    )

# fitting: feed new-names
hiMixedTripletStepTracks                 = RecoTracker.IterativeTracking.MixedTripletStep_cff.mixedTripletStepTracks.clone(
    src                 = 'hiMixedTripletStepTrackCandidates',
    #AlgorithmName = cms.string('iter8'),
    AlgorithmName = cms.string('iter4'),
    )

# Track selection
import RecoHI.HiTracking.hiMultiTrackSelector_cfi
hiMixedTripletStepSelector = RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiMultiTrackSelector.clone(
    src='hiMixedTripletStepTracks',
    trackSelectors= cms.VPSet(
    RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiLooseMTS.clone(
    name = 'hiMixedTripletStepLoose',
    d0_par2 = [9999.0, 0.0],
    dz_par2 = [9999.0, 0.0],
    applyAdaptedPVCuts = False
    ), #end of pset
    RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiTightMTS.clone(
    name = 'hiMixedTripletStepTight',
    preFilterName = 'hiMixedTripletStepLoose',
    d0_par2 = [9999.0, 0.0],
    dz_par2 = [9999.0, 0.0],
    applyAdaptedPVCuts = False
    ),
    RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiHighpurityMTS.clone(
    name = 'hiMixedTripletStep',
    preFilterName = 'hiMixedTripletStepTight',
    d0_par2 = [9999.0, 0.0],
    dz_par2 = [9999.0, 0.0],
    applyAdaptedPVCuts = False
    ),
    ) #end of vpset
    ) #end of clone  

hiMixedTripletStep = cms.Sequence(hiMixedTripletStepClusters*
                                       hiMixedTripletStepSeedsA*
                                       hiMixedTripletStepSeedsB*
                                       hiMixedTripletStepSeeds*
                                       hiMixedTripletStepTrackCandidates*
                                       hiMixedTripletStepTracks*
                                       hiMixedTripletStepSelector
                                       )
