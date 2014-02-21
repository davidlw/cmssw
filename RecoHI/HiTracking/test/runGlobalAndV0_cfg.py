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
 'root://xrootd1.cmsaf.mit.edu//store/user/davidlw/hireco_2011.root'
# '/store/data/HIRun2011/HIHighPt/RECO/hiHighPt-PromptSkim-v1/0000/96CC0844-A315-E111-A461-842B2B6AEE8B.root'
    ))

process.Timing = cms.Service("Timing")

# Number of events we want to process, -1 = all events
process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(21))


#####################################################################################
# Load some general stuff
#####################################################################################

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff')
# Data Global Tag 44x 
process.GlobalTag.globaltag = 'GR_P_V27A::All'

#Trigger Selection
### Comment out for the timing being assuming running on secondary dataset with trigger bit selected already
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHIMinBias = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
#process.hltHIMinBias.HLTPaths = ['HLT_HIMinBiasHfOrBSC_*'] # for allphysics
process.hltHIMinBias.HLTPaths = ['*'] # for allphysics
process.hltHIMinBias.andOr = cms.bool(True)
process.hltHIMinBias.throw = cms.bool(False)

# load centrality
from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)
process.HeavyIonGlobalParameters = cms.PSet(
	centralityVariable = cms.string("HFtowers"),
#	nonDefaultGlauberModel = cms.string("Hydjet_2760GeV"),
	centralitySrc = cms.InputTag("hiCentrality")
	)

#process.hiCentrality.pixelBarrelOnly = False

process.load("RecoHI.HiCentralityAlgos.CentralityFilter_cfi")
#process.centralityFilter.selectedBins = [4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40]
process.centralityFilter.selectedBins = [0]

process.eventFilter = cms.Sequence (process.collisionEventSelection * process.hltHIMinBias * process.centralityFilter)

process.load("RecoVertex.V0Producer.generalV0Candidates_cff")
process.generalV0CandidatesHI = process.generalV0Candidates.clone (
    trackRecoAlgorithm = cms.InputTag('hiGeneralAndRegitTracks'),
    vertexRecoAlgorithm = cms.InputTag('hiSelectedVertex'),
    trackQualities = cms.vstring(''),
    tkNhitsCut = cms.int32(0),
    tkChi2Cut = cms.double(10000.0),
    tkDCACut = cms.double(0.5),
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

process.ana = cms.EDAnalyzer("V0algo",
)

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
process.hiDetachedTripletStepSeeds.RegionFactoryPSet.RegionPSet.ptMin = 0.7 
process.hiDetachedTripletStepSeeds.RegionFactoryPSet.RegionPSet.originRadius = 1.0 
process.hiDetachedTripletStepSeeds.RegionFactoryPSet.RegionPSet.fixedError = 0.5

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
else:
    process.generalV0CandidatesHI.trackRecoAlgorithm = 'hiGeneralTracks'

process.ana_step = cms.Path(   
    process.eventFilter *
    process.generalV0CandidatesHI *
    process.ana
)

#####################################################################################
# Edm Output
#####################################################################################
#process.load('Configuration.EventContent.EventContentHeavyIons_cff')
#process.RECOoutput = cms.OutputModule("PoolOutputModule",
#    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
#    outputCommands = process.RECOEventContent.outputCommands,
#    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('trackRecoAndSelection')),
#    fileName = cms.untracked.string('test.root')
#)
#process.RECOoutput.outputCommands.extend(cms.untracked.vstring('drop *_*_*_*'))
#process.RECOoutput.outputCommands.extend(cms.untracked.vstring('keep *_hiGeneralAndRegitTracks_*_*'))
#process.RECOoutput.outputCommands.extend(cms.untracked.vstring('keep *_hiSelectedVertex_*_*'))
#process.RECOoutput.outputCommands.extend(cms.untracked.vstring('keep *_hiFullGeneralTracks_*_*'))
#process.RECOoutput.outputCommands.extend(cms.untracked.vstring('keep *_offlineBeamSpot_*_*'))

#process.output = cms.EndPath(process.RECOoutput)
