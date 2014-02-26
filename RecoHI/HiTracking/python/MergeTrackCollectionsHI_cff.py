import FWCore.ParameterSet.Config as cms

import RecoTracker.FinalTrackSelectors.trackListMerger_cfi
#hiGeneralTracks = RecoTracker.FinalTrackSelectors.trackListMerger_cfi.trackListMerger.clone(
#    TrackProducers = (cms.InputTag('hiGlobalPrimTracks'),
#                      cms.InputTag('hiSecondPixelTripletGlobalPrimTracks'),
#                      cms.InputTag('hiPixelPairGlobalPrimTracks')),
#    hasSelector=cms.vint32(1,1,1),
#    selectedTrackQuals = cms.VInputTag(
#    cms.InputTag("hiInitialStepSelector","hiInitialStep"),
#    cms.InputTag("hiSecondPixelTripletStepSelector","hiSecondPixelTripletStep"),
#    cms.InputTag("hiPixelPairStepSelector","hiPixelPairStep"),
#    ),                    
#    setsToMerge = cms.VPSet( cms.PSet( tLists=cms.vint32(0,1,2), pQual=cms.bool(True)),  # should this be False?
#                             ),
#    copyExtras = True,
#    makeReKeyedSeeds = cms.untracked.bool(False)
#    )

#Wei's modification
hiGeneralTracks = RecoTracker.FinalTrackSelectors.trackListMerger_cfi.trackListMerger.clone(
    TrackProducers = (cms.InputTag('hiGlobalPrimTracks'),
                      cms.InputTag('hiSecondPixelTripletGlobalPrimTracks'),
                      cms.InputTag('hiPixelPairGlobalPrimTracks'),
                      cms.InputTag('hiDetachedTripletStepTracks'),
                      cms.InputTag('hiMixedTripletStepTracks')),

    hasSelector=cms.vint32(1,1,1,1,1),
    selectedTrackQuals = cms.VInputTag(
    cms.InputTag("hiInitialStepSelector","hiInitialStep"),
    cms.InputTag("hiSecondPixelTripletStepSelector","hiSecondPixelTripletStep"),
    cms.InputTag("hiPixelPairStepSelector","hiPixelPairStep"),
    cms.InputTag("hiDetachedTripletStepSelector","hiDetachedTripletStep"),
    cms.InputTag("hiMixedTripletStepSelector","hiMixedTripletStep"),
    ),
    setsToMerge = cms.VPSet( cms.PSet( tLists=cms.vint32(0,1,2,3,4), pQual=cms.bool(True)),  # should this be False?
                             ),
    copyExtras = True,
    makeReKeyedSeeds = cms.untracked.bool(False)
    )
