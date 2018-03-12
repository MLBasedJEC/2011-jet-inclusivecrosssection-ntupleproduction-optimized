
# Forked from SMPJ Analysis Framework
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/SMPJAnalysisFW
# https://github.com/cms-smpj/SMPJ/tree/v1.0
# (further optimized to improve performance)


## Skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *
import FWCore.Utilities.FileUtils as FileUtils
import sys

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# True : when running in OpenData virtual machine
# False: when runing in lxplus 
runOnVM = False

# Local input
fileList = FileUtils.loadListFromFile('CMS_MonteCarlo2012_Summer12_DR53X_QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6_AODSIM_NoPileUp_START53_V7N-v1_combined_file_index.txt')
#fileList = FileUtils.loadListFromFile('CMS_MonteCarlo2011_Summer11LegDR_QCD_Pt-80to120_TuneZ2_7TeV_pythia6_AODSIM_PU_S13_START53_LV6-v1_00000_file_index.txt')
#fileList = ['file:04C05FC3-35B7-E311-9924-003048679182.root']
process.source.fileNames = cms.untracked.vstring(*fileList)

if runOnVM:
    process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/START53_V27.db')

# Global tag for Summer12DR53X-NoPileUp_START53__V7N-v1
process.GlobalTag.globaltag = cms.string('START53_V27::All')

# Select good vertices
process.goodOfflinePrimaryVertices = cms.EDFilter(
    "VertexSelector",
    filter = cms.bool(False),
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
    )

# -------- The Tracking failure filter ------#
from RecoMET.METFilters.trackingFailureFilter_cfi import trackingFailureFilter
process.trackingFailureFilter = trackingFailureFilter.clone()
process.trackingFailureFilter.VertexSource = cms.InputTag('goodOfflinePrimaryVertices')

# Load jet correction services for all jet algoritms
process.load("JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff")

################### EDAnalyzer ##############################3
algsize = {'ak':[5,7]} #Options: ak
jettype = ['PF'] #Options: PF, PFchs, PFPuppi
corrs = ['L2L3'] #Options: L1Fast, L2L3, L2L3Residual
jetCollections = []
jetCorrections = []
for k, v in algsize.iteritems():
    for s in v:
        for j in jettype:
            jetCollections.append((str(k+str(s)+j+"Jets"),str(k.upper()+str(s)+j),str(k+str(s)+"GenJets")))
            corrStr = ""
            for c in corrs:
                corrStr += c
            jetCorrections.append(str(k+str(s)+j+corrStr))

for ic, collection in enumerate(jetCollections):
    treeProducer = cms.EDAnalyzer('OpenDataTreeProducerOptimized',
        ## numpy output                                                                                                          
        maxRows = cms.untracked.int32(10000000),
        ## jet collections ###########################
        jets            = cms.InputTag(collection[0]),
        ## MET collection ####
        pfmet           = cms.InputTag('pfMET7'),
        ## database entry for the uncertainties ######
        PFPayloadName   = cms.string(collection[1]),
    
        ## set the conditions for good Vtx counting ##
        offlineVertices = cms.InputTag('goodOfflinePrimaryVertices'),
        goodVtxNdof     = cms.double(4), 
        goodVtxZ        = cms.double(24),
        ## rho #######################################
        srcPFRho        = cms.InputTag('kt6PFJets','rho'),
        ## preselection cuts #########################
        maxY            = cms.double(99.0), 
        maxEta          = cms.double(10.0), 
        minPFPt         = cms.double(10.),
        minNPFJets      = cms.int32(1),
        minGenPt        = cms.untracked.double(30),
        minJJMass       = cms.double(-1),
        isMCarlo        = cms.untracked.bool(True),
        genjets         = cms.untracked.InputTag(collection[2]),
        genparticles    = cms.untracked.InputTag('genParticles'),
        useGenInfo      = cms.untracked.bool(True),
        ## trigger ###################################
        printTriggerMenu = cms.untracked.bool(True),
        processName     = cms.string('HLT'),
        triggerNames    = cms.vstring(
                                    'HLT_Jet30', 'HLT_Jet60', 'HLT_Jet80', 'HLT_Jet110', 
                                    'HLT_Jet150','HLT_Jet190','HLT_Jet240','HLT_Jet370',
                                    ),
        triggerResults  = cms.InputTag("TriggerResults","","HLT"),
        triggerEvent    = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
        ## jet energy correction labels ##############
        jetCorr          = cms.string(jetCorrections[ic]),
        # PF Candidates
        pfCandidates     = cms.InputTag("particleFlow","","RECO"),
    )

    ############# hlt filter #########################
    process.hltFilter = cms.EDFilter('HLTHighLevel',
        TriggerResultsTag  = cms.InputTag('TriggerResults','','HLT'),
        HLTPaths           = cms.vstring('HLT_Jet*', 'HLT_DiJetAve*'),
        eventSetupPathsKey = cms.string(''),
        andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
        throw              = cms.bool(False)
    )

    setattr(process,treeProducer.jets.getModuleLabel(),treeProducer)

    # Let it run
    path = cms.Path(
        process.goodOfflinePrimaryVertices*
        process.hltFilter *
        process.trackingFailureFilter *
        treeProducer
    )
    setattr(process, treeProducer.jets.getModuleLabel() + 'Path', path)


# Approximate processing time on VM (Intel Core i5-5300U 2.3GHz laptop):
# 50000 events per 1 hour (both for DATA and MC)

# Change number of events here:
process.maxEvents.input = 10

process.MessageLogger.cerr.FwkReport.reportEvery = 500

# Output file
process.TFileService = cms.Service("TFileService", fileName = cms.string('OpenDataTree_mc.root'))

# To suppress long output at the end of the job
#process.options.wantSummary = False   

del process.outpath
