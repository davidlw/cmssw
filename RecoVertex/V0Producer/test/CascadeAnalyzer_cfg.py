import FWCore.ParameterSet.Config as cms

process = cms.Process("V0Val")
process.load("FWCore.MessageService.MessageLogger_cfi")

### standard includes
process.load('Configuration.StandardSequences.Generator_cff')
process.load('GeneratorInterface.HiGenCommon.VtxSmearedRealisticPPbBoost8TeVCollision_cff')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load("Configuration.StandardSequences.Digi_cff")
process.load("Configuration.StandardSequences.DigiToRaw_cff")
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_P_V43F::All'

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
'/store/hidata/HIRun2013A/PAHighPt/RECO/PromptReco-v1/000/210/634/FA4E6B7E-7366-E211-8DD0-0019B9F581C9.root'
#'root://xrootd1.cmsaf.mit.edu//store/user/vzhukova/HIJING_MB/HIJING_MB_RECO_v3/13a591fee6315e7fb1e99e6ba8e52eaa/reco_hijing_2582_1_tV4.root'
#'root://xrootd1.cmsaf.mit.edu//store/user/vzhukova/HIJING_MB/HIJING_MB_RECO_v3/13a591fee6315e7fb1e99e6ba8e52eaa/reco_hijing_2582_1_tV4.root'
#'root://xrootd1.cmsaf.mit.edu//store/user/vzhukova/HIJING_MB/HIJING_MB_RECO_v3/13a591fee6315e7fb1e99e6ba8e52eaa/reco_hijing_1000_1_Bov.root'
#'/store/user/vzhukova/HYDGET_PERIPH_batch7/HYDGET_PERIPH_RECO_batch7/b7d33bba7673cdb1ee6f4983c0800c79/HYDGET_PERIPH_RECO_10_1_7bq.root'
#'root://xrootd1.cmsaf.mit.edu//store/user/vzhukova/HIJING_GEN-SIM_YUE-SHI_Minbias_2_v1/HIJING_RECO_YUE-SHI_Minbias__2_v1/b7d33bba7673cdb1ee6f4983c0800c79/hijing_reco_fix_4_2_RU0.root'
#'root://xrootd1.cmsaf.mit.edu//store/himc/HiWinter13/Hijing_PPb502_MinimumBias/GEN-SIM-RECO/pa_STARTHI53_V25-v1/30000/64B44ED9-EF77-E211-825A-00266CF97FF4.root'
                )
                            )
process.load("RiceHIG.V0Analysis.v0selector_cff")
process.load("RecoVertex.V0Producer.generalCascadeCandidates_cff")
process.generalDdCandidates.v0RecoAlgorithm = cms.InputTag('selectV0CandidatesNewkshort:Kshort')

process.generalV0CandidatesNew = process.generalV0Candidates.clone (
    tkNhitsCut = cms.int32(0),
    tkChi2Cut = cms.double(7.0),
    dauTransImpactSigCut = cms.double(1.0),
    dauLongImpactSigCut = cms.double(1.0),
    xiVtxSignificance3DCut = cms.double(0.0),
    xiVtxSignificance2DCut = cms.double(0.0),
    vtxSignificance2DCut = cms.double(0.0),
    vtxSignificance3DCut = cms.double(4.0)
)   

process.selectV0CandidatesNewlambda.v0CollName = cms.string("generalV0CandidatesNew")
process.selectV0CandidatesNewkshort.v0CollName = cms.string("generalV0CandidatesNew")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

# Additional output definition
#process.TFileService = cms.Service("TFileService",
#                                   fileName = cms.string('v0validation.root')
#                                   )

process.charm_seq = cms.Sequence(process.generalV0CandidatesNew*process.selectV0CandidatesNewlambda*process.selectV0CandidatesNewkshort*process.generalDdCandidates)

process.p = cms.Path(process.charm_seq)

#####################################################################################
# Edm Output
#####################################################################################
#process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.RECOoutput = cms.OutputModule("PoolOutputModule",
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    outputCommands = process.RECOEventContent.outputCommands,
#    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('trackRecoAndSelection')),
    fileName = cms.untracked.string('test.root')
)
process.RECOoutput.outputCommands.extend(cms.untracked.vstring('drop *_*_*_*'))
process.RECOoutput.outputCommands.extend(cms.untracked.vstring('keep *_*Dd*_*_*'))

process.output = cms.EndPath(process.RECOoutput)

#process.schedule = cms.Schedule(process.p)
