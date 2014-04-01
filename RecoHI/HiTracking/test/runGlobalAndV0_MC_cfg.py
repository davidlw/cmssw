import FWCore.ParameterSet.Config as cms

process = cms.Process('TRACKATTACK')

doLargeD0=True
doRegit=False
rawORreco=True
isEmbedded=False

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

#####################################################################################
# Input source
#####################################################################################

process.source = cms.Source("PoolSource",
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                            fileNames = cms.untracked.vstring(
# '/store/himc/HiFall13/Hydjet1p8_TuneDrum_Quenched_MinBias_2760GeV/GEN-SIM/STARTHI53_V28-v2/00000/00699AE5-5A5E-E311-83B4-008CFA007B98.root'
#'file:/net/hisrv0001/home/davidlw/scratch1/RECO_MC_53x.root'
       'root://xrootd1.cmsaf.mit.edu//store/user/davidlw/Hydjet1p8_TuneDrum_Quenched_MinBias_2760GeV/RECO_v1/0e6ab43111ff2b7e6299139d8f86b661/reco_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_100_1_KCp.root'
    ))

process.Timing = cms.Service("Timing")

# Number of events we want to process, -1 = all events
process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(82))


#####################################################################################
# Load some general stuff
#####################################################################################

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
#process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff')
# Data Global Tag 44x 
#process.GlobalTag.globaltag = 'GR_P_V27A::All'
# Data Global Tag 53x 
#process.GlobalTag.globaltag = 'GR_P_V43F::All'
# MC Global Tag 53x 
process.GlobalTag.globaltag = 'STARTHI53_V17::All'

#Trigger Selection
### Comment out for the timing being assuming running on secondary dataset with trigger bit selected already
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHIMinBias = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
#process.hltHIMinBias.HLTPaths = ['HLT_HIMinBiasHfOrBSC_*'] # for allphysics
process.hltHIMinBias.HLTPaths = ['*'] # for allphysics
process.hltHIMinBias.andOr = cms.bool(True)
process.hltHIMinBias.throw = cms.bool(False)

# load centrality
from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import *
overrideCentrality(process)
process.HeavyIonGlobalParameters = cms.PSet(
	centralityVariable = cms.string("HFtowers"),
	nonDefaultGlauberModel = cms.string("Hydjet_Drum"),
	centralitySrc = cms.InputTag("hiCentrality")
	)

#process.hiCentrality.pixelBarrelOnly = False

process.load("RecoHI.HiCentralityAlgos.CentralityFilter_cfi")
#process.centralityFilter.selectedBins = [4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40]
process.centralityFilter.selectedBins = [29,30,31,32,33,34,35,36,37,38,39,40]

#process.eventFilter = cms.Sequence (process.collisionEventSelection * process.hltHIMinBias * process.centralityFilter)
process.eventFilter = cms.Sequence ( process.centralityFilter )

process.load("RecoVertex.V0Producer.generalV0Candidates_cff")
process.generalV0CandidatesHI = process.generalV0Candidates.clone (
    trackRecoAlgorithm = cms.InputTag('hiGeneralAndRegitTracks'),
    vertexRecoAlgorithm = cms.InputTag('hiSelectedVertex'),
    trackQualities = cms.vstring('loose'),
    tkNhitsCut = cms.int32(0),
    tkChi2Cut = cms.double(10000.0),
    tkDCACut = cms.double(1.0),
    vtxChi2Cut = cms.double(10000.0),
    vtxSignificance3DCut = cms.double(0.0),
    vtxSignificance2DCut = cms.double(0.0),
    dauTransImpactSigCut = cms.double(1.0),
    dauLongImpactSigCut = cms.double(1.0),
    collinearityCut = cms.double(0.0),
    innerHitPosCut = cms.double(-1),
    selectXis = cms.bool(False),
    selectOmegas = cms.bool(False)
)

### validation-specific includes
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
process.load("RiceHIG.V0Analysis.v0selector_cff")
process.load("RiceHIG.V0Analysis.v0validator_cff")

process.v0ValidatorHI.kShortCollection = cms.InputTag('selectV0CandidatesNewkshort:Kshort')
process.v0ValidatorHI.lambdaCollection = cms.InputTag('selectV0CandidatesNewlambda:Lambda')
process.v0ValidatorHI.vertexCollection = cms.InputTag('hiSelectedVertex')
#process.v0ValidatorHI.isMatchByHitsOrChi2 = cms.bool(True)
process.TrackAssociatorByHits.Cut_RecoToSim = cms.double(0.5)

process.selectV0CandidatesNewlambda.v0CollName = cms.string("generalV0CandidatesHI")
process.selectV0CandidatesNewkshort.v0CollName = cms.string("generalV0CandidatesHI")
process.selectV0CandidatesNewlambda.vertexCollName = cms.InputTag("hiSelectedVertex")
process.selectV0CandidatesNewkshort.vertexCollName = cms.InputTag("hiSelectedVertex")
process.selectV0CandidatesNewkshort.cosThetaCut = cms.double(0.99)
process.selectV0CandidatesNewkshort.decayLSigCut = cms.double(3.0)
process.selectV0CandidatesNewkshort.misIDMassCut   = cms.double(0.015)
process.selectV0CandidatesNewkshort.misIDMassCutEE = cms.double(0.015)
process.selectV0CandidatesNewlambda.cosThetaCut = cms.double(0.99)
process.selectV0CandidatesNewlambda.decayLSigCut = cms.double(3.0)
process.selectV0CandidatesNewlambda.misIDMassCut   = cms.double(0.015)
process.selectV0CandidatesNewlambda.misIDMassCutEE = cms.double(0.015)

process.v0validation = cms.Sequence(process.selectV0CandidatesNewlambda*process.selectV0CandidatesNewkshort*process.v0ValidatorHI)

#process.ana = cms.EDAnalyzer("V0algo",
#)

#####################################################################################
# Define tree output
#####################################################################################

process.TFileService = cms.Service("TFileService",
                                  fileName=cms.string("v0AnalyzerNT.root"))

#####################################################################################
# Additional Reconstruction 
#####################################################################################


# redo reco or just tracking

if rawORreco:
    process.rechits = cms.Sequence(process.siPixelRecHits * process.siStripMatchedRecHits)
    process.hiTrackReco = cms.Sequence(process.rechits * process.heavyIonTracking)

    process.trackRecoAndSelection = cms.Path(
        process.eventFilter*
        process.hiTrackReco
        )
    
else:
    process.reco_extra = cms.Path(
        process.eventFilter*
        process.RawToDigi * process.reconstructionHeavyIons)
    

    
# tack on iteative tracking, particle flow and calo-matching

#iteerative tracking
process.load("RecoHI.HiTracking.hiIterTracking_cff")
process.hiDetachedTripletStepSeeds.RegionFactoryPSet.RegionPSet.ptMin = 0.2
process.hiDetachedTripletStepSeeds.RegionFactoryPSet.RegionPSet.originRadius = 1.5
process.hiDetachedTripletStepSeeds.RegionFactoryPSet.RegionPSet.fixedError = 5 
process.hiMixedTripletStepSeedsA.RegionFactoryPSet.RegionPSet.ptMin = 0.2
process.hiMixedTripletStepSeedsA.RegionFactoryPSet.RegionPSet.originRadius = 1.5
process.hiMixedTripletStepSeedsA.RegionFactoryPSet.RegionPSet.fixedError = 5
process.hiMixedTripletStepSeedsB.RegionFactoryPSet.RegionPSet.ptMin = 0.4
process.hiMixedTripletStepSeedsB.RegionFactoryPSet.RegionPSet.originRadius = 1.5
process.hiMixedTripletStepSeedsB.RegionFactoryPSet.RegionPSet.fixedError = 5

if doLargeD0:
    process.heavyIonTracking *= process.hiFullIterTracking
else:
    process.heavyIonTracking *= process.hiIterTracking

# Now do more tracking around the jets
if doRegit:
    process.load("RecoHI.HiTracking.hiRegitTracking_cff")
 
    process.selectedJets = cms.EDFilter("CandViewSelector",
       cut = cms.string  ('pt > 20.0'),
       src = cms.InputTag("iterativeConePu5CaloJets"),
    )
   
    process.hiRegitInitialStepSeeds.RegionFactoryPSet.RegionPSet.JetSrc = cms.InputTag("selectedJets")
    process.hiRegitLowPtTripletStepSeeds.RegionFactoryPSet.RegionPSet.JetSrc = cms.InputTag("selectedJets")
    process.hiRegitPixelPairStepSeeds.RegionFactoryPSet.RegionPSet.JetSrc = cms.InputTag("selectedJets")
    process.hiRegitDetachedTripletStepSeeds.RegionFactoryPSet.RegionPSet.JetSrc = cms.InputTag("selectedJets")
    process.hiRegitMixedTripletStepSeedsA.RegionFactoryPSet.RegionPSet.JetSrc = cms.InputTag("selectedJets")
    process.hiRegitMixedTripletStepSeedsB.RegionFactoryPSet.RegionPSet.JetSrc = cms.InputTag("selectedJets")

#    process.hiRegitInitialStepSeeds.RegionFactoryPSet.RegionPSet.JetSrc = cms.InputTag("iterativeConePu5CaloJets")
#    process.hiRegitLowPtTripletStepSeeds.RegionFactoryPSet.RegionPSet.JetSrc = cms.InputTag("iterativeConePu5CaloJets")
#    process.hiRegitPixelPairStepSeeds.RegionFactoryPSet.RegionPSet.JetSrc = cms.InputTag("iterativeConePu5CaloJets")
#    process.hiRegitDetachedTripletStepSeeds.RegionFactoryPSet.RegionPSet.JetSrc = cms.InputTag("iterativeConePu5CaloJets")
#    process.hiRegitMixedTripletStepSeedsA.RegionFactoryPSet.RegionPSet.JetSrc = cms.InputTag("iterativeConePu5CaloJets")
#    process.hiRegitMixedTripletStepSeedsB.RegionFactoryPSet.RegionPSet.JetSrc = cms.InputTag("iterativeConePu5CaloJets")

    # merged with the global, iterative tracking
    process.load("RecoHI.HiTracking.MergeRegit_cff")

    process.regionalTracking = cms.Path(
        process.eventFilter *
        process.selectedJets *
        process.hiRegitTracking *
        process.hiGeneralAndRegitTracks
        )

if doRegit:
    process.generalV0CandidatesHI.trackRecoAlgorithm = 'hiGeneralAndRegitTracks'
    process.v0ValidatorHI.trackCollection = cms.InputTag('hiGeneralAndRegitTracks')
else:
    process.generalV0CandidatesHI.trackRecoAlgorithm = 'hiGeneralTracks'
    process.v0ValidatorHI.trackCollection = cms.InputTag('hiGeneralTracks')

process.ana_step = cms.Path(   
    process.eventFilter *
    process.generalV0CandidatesHI *
    process.v0validation
)

#####################################################################################
# Edm Output
#####################################################################################
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.RECOoutput = cms.OutputModule("PoolOutputModule",
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    outputCommands = process.RECOEventContent.outputCommands,
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('trackRecoAndSelection')),
    fileName = cms.untracked.string('test.root')
)
process.RECOoutput.outputCommands.extend(cms.untracked.vstring('drop *_*_*_*'))
process.RECOoutput.outputCommands.extend(cms.untracked.vstring('keep *_hiGeneral*_*_*'))
process.RECOoutput.outputCommands.extend(cms.untracked.vstring('keep *_hiSelectedVertex_*_*'))
process.RECOoutput.outputCommands.extend(cms.untracked.vstring('keep *_*V0*_*_*'))

#process.output = cms.EndPath(process.RECOoutput)
