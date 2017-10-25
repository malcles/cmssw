from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

#config.General.requestName = 'DoubleMuon-Run2017A-PromptReco-NewMini-v2'
#config.General.requestName = 'DoubleMuon-Run2017A-PromptReco-NewMini-v1'

#config.General.requestName = 'DoubleMuon-Run2017B-PromptReco-NewMini-v1'

#config.General.requestName = 'DoubleMuon-Run2017A-PromptReco-NewMini-v3'
#config.General.requestName = 'DoubleMuon-Run2017B-PromptReco-NewMini-v2'
#config.General.requestName = 'DoubleMuon-Run2017C-PromptReco-NewMini-v1'
#config.General.requestName = 'DoubleMuon-Run2017C-PromptReco-NewMini-v2'
#config.General.requestName = 'DoubleMuon-Run2017C-PromptReco-NewMini-v3'
config.General.requestName = 'DoubleMuon-Run2017D-PromptReco-NewMini-v1'


config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'

config.JobType.psetName = 'Validation/RecoVertex/python/miniAOD-prod_PAT_DATA_AllVerticesCollFiltered.py'

#config.Data.inputDataset = '/DoubleMuon/Run2017A-PromptReco-v2/AOD'
#config.Data.inputDataset = '/DoubleMuon/Run2017A-PromptReco-v1/AOD'
config.Data.inputDataset = '/DoubleMuon/Run2017B-PromptReco-v1/AOD'

#config.Data.inputDataset = '/DoubleMuon/Run2017A-PromptReco-v3/AOD'
#config.Data.inputDataset = '/DoubleMuon/Run2017B-PromptReco-v2/AOD'
#config.Data.inputDataset = '/DoubleMuon/Run2017C-PromptReco-v1/AOD'
#config.Data.inputDataset = '/DoubleMuon/Run2017C-PromptReco-v2/AOD'
#config.Data.inputDataset = '/DoubleMuon/Run2017C-PromptReco-v3/AOD'
config.Data.inputDataset = '/DoubleMuon/Run2017D-PromptReco-v1/AOD'

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 20


config.Data.lumiMask = config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/PromptReco/Cert_294927-304120_13TeV_PromptReco_Collisions17_JSON.txt'

config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())

config.Data.publication = True
#config.Data.outputDatasetTag = 'DoubleMuon-Run2017A-PromptReco-NewMini-v2'
#config.Data.outputDatasetTag = 'DoubleMuon-Run2017A-PromptReco-NewMini-v1'
#config.Data.outputDatasetTag = 'DoubleMuon-Run2017B-PromptReco-NewMini-v1'

#config.Data.outputDatasetTag = 'DoubleMuon-Run2017A-PromptReco-NewMini-v3'
#config.Data.outputDatasetTag = 'DoubleMuon-Run2017B-PromptReco-NewMini-v2'
#config.Data.outputDatasetTag = 'DoubleMuon-Run2017C-PromptReco-NewMini-v1'
#config.Data.outputDatasetTag = 'DoubleMuon-Run2017C-PromptReco-NewMini-v2'
#config.Data.outputDatasetTag = 'DoubleMuon-Run2017C-PromptReco-NewMini-v3'
config.Data.outputDatasetTag = 'DoubleMuon-Run2017D-PromptReco-NewMini-v1'



config.Site.storageSite = 'T2_FR_GRIF_IRFU'
