import FWCore.ParameterSet.Config as cms

from RecoVertex.V0Producer.generalCascadeCandidates_cfi import *

generalXiCandidates = generalCascadeCandidates.clone(
    v0RecoAlgorithm = cms.InputTag('generalV0Candidates','Lambda'),
    batDauTrkMass = cms.double(0.139570)
)

generalOmegaCandidates = generalCascadeCandidates.clone(
    v0RecoAlgorithm = cms.InputTag('generalV0Candidates','Lambda'),
    batDauTrkMass = cms.double(0.493677)
)

generalDdCandidates = generalCascadeCandidates.clone(
    v0RecoAlgorithm = cms.InputTag('generalV0Candidates','Kshort'),
    batDauTrkMass = cms.double(0.139570)
)

generalDsCandidates = generalCascadeCandidates.clone(
    v0RecoAlgorithm = cms.InputTag('generalV0Candidates','Kshort'),
    batDauTrkMass = cms.double(0.493677)
)
