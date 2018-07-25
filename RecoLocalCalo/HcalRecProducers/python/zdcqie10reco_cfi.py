import FWCore.ParameterSet.Config as cms

zdcRecHits = cms.EDProducer("ZDCQIE10Reconstructor",
    digiLabel = cms.InputTag("hcalDigis","ZDC"),
    dropZSmarkedPassed = cms.bool(False),
    tsFromDB = cms.bool(False),
    tsToStart = cms.int32(4),
    satCorrFactor = cms.double(1.5)
)
