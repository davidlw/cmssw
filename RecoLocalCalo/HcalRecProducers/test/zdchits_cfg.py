process = cms.Process("ZDCTest")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1))

#
#   Command Line Input(Copied from DQM for now)
#
import sys

#
#   Change the filename to process
#
runNumber = sys.argv[2]

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/q/qwang/work/public/ZDC2018/example/HFanalysis_'+runNumber+'.root'
    ),
    labelRawDataLikeMC = cms.untracked.bool(False)
)

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(False)
        )
#
#   For Debugging: Create a Pool Output Module
#
process.output = cms.OutputModule(
        'PoolOutputModule',
        fileName = cms.untracked.string('ZDCRecHits_'+runNumber+'.root')
)

process.load('Configuration.Geometry.GeometryIdeal_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.autoCond import autoCond
#from CondCore.DBCommon.CondDBSetup_cfi import *
from CondCore.CondDB.CondDB_cfi import *

process.GlobalTag.globaltag = autoCond['startup']

process.zdcdigi = cms.EDProducer('QWZDC2018Producer',
		Src = cms.untracked.InputTag('hcalDigis', 'ZDC')
		)

process.load("RecoLocalCalo.HcalRecProducers.zdcqie10reco_cfi")

process.p = cms.Path( process.zdcRecHits )
process.outpath = cms.EndPath(process.output)
