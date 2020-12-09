import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
from FWCore.ParameterSet.VarParsing import VarParsing
import os


options = VarParsing ('analysis')

options.register( 'executionMode',
                  0,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "0 for particle data, 1 for jet data, 2 for jet prediction"
               )

process = cms.Process("AK4jets")



# QG likelihood
process.load("JetNtupleProducerTool.JetAnalyzer.QGLikelihood_cfi")
process.load("RecoJets.JetProducers.QGTagger_cfi")
process.QGTagger.srcJets = cms.InputTag("slimmedJets")
process.QGTagger.jetsLabel = cms.string("QGL_AK4PFchs")

# File service
process.load("CommonTools.UtilAlgos.TFileService_cfi")
process.TFileService.fileName=cms.string("JetNtuple_RunIISummer16_13TeV_MC.root")

# Load up the filelist
filePath=os.environ["CMSSW_BASE"]+"/src/JetNtupleProducerTool/JetAnalyzer/python/"
if(options.executionMode == 0):
    process.source = cms.Source("PoolSource",
	## Process whole data set
	#fileNames = cms.untracked.vstring(*fileList)
	## Process just one file
	fileNames = cms.untracked.vstring(fileList[0])
    )
elif(options.executionMode == 1):
    process.source = cms.Source("PoolSource",
	## Process whole data set
	#fileNames = cms.untracked.vstring(*fileList)
	## Process just one file
	fileNames = cms.untracked.vstring(fileList[1])
    )
else:
    process.source = cms.Source("PoolSource",
	## Process whole data set
	#fileNames = cms.untracked.vstring(*fileList)
	## Process just one file
	fileNames = cms.untracked.vstring(fileList[2])
    )



process.AK4jets = cms.EDAnalyzer("JetAnalyzer",
	## jet, PF and generator level collections ##
	jets = cms.InputTag("slimmedJets"),
	pfCands = cms.InputTag("packedPFCandidates"),
	genJets = cms.InputTag("slimmedGenJets"),
	genEventInfo = cms.InputTag("generator"),
	## good primary vertices ##
	vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
	confGoodVtxNdof = cms.double(4),
	confGoodVtxZ = cms.double(24),
	confGoodVtxRho = cms.double(2),
	## pileup and rhos ##
	pileupInfo = cms.InputTag("slimmedAddPileupInfo"),
	pfRhoAll = cms.InputTag("fixedGridRhoFastjetAll"),
	pfRhoCentral = cms.InputTag("fixedGridRhoFastjetCentral"),
	pfRhoCentralNeutral = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
	pfRhoCentralChargedPileUp = cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp"),
    executionMode = cms.untracked.int(options.executionMode)
)

# Choose how many events to process (-1 = all)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Report execution progress
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.p = cms.Path(process.QGTagger + process.AK4jets)
