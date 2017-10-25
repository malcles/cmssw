# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: miniAOD-prod -s PAT --eventcontent MINIAOD --runUnscheduled --data --filein /store/data/Run2017A/DoubleMuon/AOD/PromptReco-v2/000/296/172/00000/4AAD2027-544C-E711-802C-02163E01A515.root --conditions 92X_dataRun2_Prompt_v4 -n 100 --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('PAT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.PAT_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'root://cms-xrd-global.cern.ch//store/data/Run2017A/DoubleMuon/AOD/PromptReco-v2/000/296/172/00000/4AAD2027-544C-E711-802C-02163E01A515.root'
        'root://cms-xrd-global.cern.ch//store/data/Run2017D/DoubleMuon/AOD/PromptReco-v1/000/302/031/00000/02367FB0-4C8F-E711-A7D0-02163E0135EB.root',
#        '/store/data/Run2017D/DoubleMuon/AOD/PromptReco-v1/000/302/031/00000/064FB49F-428F-E711-A18B-02163E011BF8.root',
#        '/store/data/Run2017D/DoubleMuon/AOD/PromptReco-v1/000/302/031/00000/0C2BD442-398F-E711-97DE-02163E019DE8.root',
),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('miniAOD-prod nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

#added JM: beg

from Validation.RecoVertex.myNoMuonTrackProducer_cfi import *
process.myNoMuonTrackProducerNoMu=myNoMuonTrackProducer.clone()
process.myNoMuonTrackProducerNoMu.doRemoveMuons=cms.untracked.bool(True)
process.myNoMuonTrackProducerNoMu.verbose=cms.untracked.bool(True)
process.myNoMuonTrackProducerNoMu.isData=cms.untracked.bool(True)

from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi import *
process.offlinePrimaryVerticesNoMu=offlinePrimaryVertices.clone()
process.offlinePrimaryVerticesNoMu.TrackLabel = cms.InputTag("myNoMuonTrackProducerNoMu")


from Validation.RecoVertex.myNoMuonTrackProducer_cfi import *
process.myNoMuonTrackProducerWithMu=myNoMuonTrackProducer.clone()
process.myNoMuonTrackProducerWithMu.doRemoveMuons=cms.untracked.bool(False)
process.myNoMuonTrackProducerWithMu.verbose=cms.untracked.bool(False)
process.myNoMuonTrackProducerWithMu.isData=cms.untracked.bool(True)

from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi import *
process.offlinePrimaryVerticesWithMu=offlinePrimaryVertices.clone()
process.offlinePrimaryVerticesWithMu.TrackLabel = cms.InputTag("myNoMuonTrackProducerWithMu")

from PhysicsTools.PatAlgos.selectionLayer1.countPatCandidates_cff import *
from PhysicsTools.PatAlgos.selectionLayer1.muonCountFilter_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.muonSelector_cfi import *
 


process.myFilter = cms.EDFilter(
    "myNoLeptonFilter",
    allTrackTag=cms.InputTag('myNoMuonTrackProducerWithMu'),
    noLepTrackTag=cms.InputTag('myNoMuonTrackProducerNoMu')
    )

process.selectionNoLeptonFilter = cms.Path((process.myNoMuonTrackProducerNoMu+process.myNoMuonTrackProducerWithMu)*process.myFilter)
process.vtxRefit=cms.Path(process.offlinePrimaryVerticesNoMu+process.offlinePrimaryVerticesWithMu)

#added JM: end

# Output definition



process.MINIAODoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string('miniAOD-prod_PAT_my.root'),
    outputCommands = process.MINIAODEventContent.outputCommands,
    overrideInputFileSplitLevels = cms.untracked.bool(True),
    SelectEvents = cms.untracked.PSet( #added JM
        SelectEvents = cms.vstring('selectionNoLeptonFilter')#added JM
    )#added JM
)

process.MINIAODoutput.outputCommands.append('keep *_offlinePrimaryVertices*_*_*') #added JM
process.MINIAODoutput.outputCommands.append('keep *_myNoMuonTrackProducer*_*_*') #added JM


# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '92X_dataRun2_Prompt_v4', '')
#92X_dataRun2_Prompt_v8 for Run2017 C

process.load('RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi')#added JM
process.offlinePrimaryVertices.verbose = cms.untracked.bool(False)#added JM

# Path and EndPath definitions
process.Flag_trackingFailureFilter = cms.Path(process.goodVertices+process.trackingFailureFilter)
process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)
process.Flag_trkPOGFilters = cms.Path(process.trkPOGFilters)
process.Flag_HcalStripHaloFilter = cms.Path(process.HcalStripHaloFilter)
process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~process.logErrorTooManyClusters)
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)
process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)
process.Flag_eeBadScFilter = cms.Path(process.eeBadScFilter)
process.Flag_METFilters = cms.Path(process.metFilters)
process.Flag_chargedHadronTrackResolutionFilter = cms.Path(process.chargedHadronTrackResolutionFilter)
process.Flag_globalTightHalo2016Filter = cms.Path(process.globalTightHalo2016Filter)
process.Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(process.CSCTightHaloTrkMuUnvetoFilter)
process.Flag_HBHENoiseIsoFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseIsoFilter)
process.Flag_BadChargedCandidateSummer16Filter = cms.Path(process.BadChargedCandidateSummer16Filter)
process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)
process.Flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)
process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseFilter)
process.Flag_trkPOG_toomanystripclus53X = cms.Path(~process.toomanystripclus53X)
process.Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(process.EcalDeadCellBoundaryEnergyFilter)
process.Flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter)
process.Flag_trkPOG_manystripclus53X = cms.Path(~process.manystripclus53X)
process.Flag_BadPFMuonSummer16Filter = cms.Path(process.BadPFMuonSummer16Filter)
process.Flag_muonBadTrackFilter = cms.Path(process.muonBadTrackFilter)
process.Flag_CSCTightHalo2015Filter = cms.Path(process.CSCTightHalo2015Filter)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.MINIAODoutput_step = cms.EndPath(process.MINIAODoutput)

# Schedule definition
process.schedule = cms.Schedule(process.Flag_HBHENoiseFilter,process.Flag_HBHENoiseIsoFilter,process.Flag_CSCTightHaloFilter,process.Flag_CSCTightHaloTrkMuUnvetoFilter,process.Flag_CSCTightHalo2015Filter,process.Flag_globalTightHalo2016Filter,process.Flag_globalSuperTightHalo2016Filter,process.Flag_HcalStripHaloFilter,process.Flag_hcalLaserEventFilter,process.Flag_EcalDeadCellTriggerPrimitiveFilter,process.Flag_EcalDeadCellBoundaryEnergyFilter,process.Flag_goodVertices,process.Flag_eeBadScFilter,process.Flag_ecalLaserCorrFilter,process.Flag_trkPOGFilters,process.Flag_chargedHadronTrackResolutionFilter,process.Flag_muonBadTrackFilter,process.Flag_BadChargedCandidateFilter,process.Flag_BadPFMuonFilter,process.Flag_BadChargedCandidateSummer16Filter,process.Flag_BadPFMuonSummer16Filter,process.Flag_trkPOG_manystripclus53X,process.Flag_trkPOG_toomanystripclus53X,process.Flag_trkPOG_logErrorTooManyClusters,process.Flag_METFilters,process.selectionNoLeptonFilter,process.vtxRefit,process.endjob_step,process.MINIAODoutput_step)
process.schedule.associate(process.patTask)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllData 

#call to customisation function miniAOD_customizeAllData imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
process = miniAOD_customizeAllData(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
