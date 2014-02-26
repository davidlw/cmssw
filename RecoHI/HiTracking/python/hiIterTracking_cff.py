from RecoHI.HiTracking.hiSecondPixelTripletStep_cff import *
from RecoHI.HiTracking.hiMixedTripletStep_cff import *
from RecoHI.HiTracking.hiPixelPairStep_cff import *
from RecoHI.HiTracking.hiDetachedTripletStep_cff import *
from RecoHI.HiTracking.hiMixedTripletStep_cff import *
from RecoHI.HiTracking.MergeTrackCollectionsHI_cff import *

hiIterTracking = cms.Sequence(
    hiSecondPixelTripletStep
    *hiPixelPairStep
    *hiGeneralTracks
    )

# Wei's modification
hiFullIterTracking = cms.Sequence(
    hiSecondPixelTripletStep
    *hiPixelPairStep
    *hiDetachedTripletStep
    *hiMixedTripletStep
    *hiGeneralTracks
    )
