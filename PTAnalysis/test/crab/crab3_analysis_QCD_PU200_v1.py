from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName   = 'QCD_PU200'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
#config.JobType.inputFiles = ['Summer16_23Sep2016V4_MC_L1FastJet_AK4PFchs.txt','Summer16_23Sep2016V4_MC_L2Relative_AK4PFchs.txt','Summer16_23Sep2016V4_MC_L3Absolute_AK4PFchs.txt','Summer16_23Sep2016V4_MC_L1FastJet_AK4PFPuppi.txt','Summer16_23Sep2016V4_MC_L2Relative_AK4PFPuppi.txt','Summer16_23Sep2016V4_MC_L3Absolute_AK4PFPuppi.txt']
# Name of the CMSSW configuration file
config.JobType.psetName    = 'electronIsolationAnalyzer_cfg.py'
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
config.Data.inputDataset = '/QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8/PhaseIISpr18AODMiniAOD-PU200_93X_upgrade2023_realistic_v5-v1/AODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 30
config.Data.totalUnits = -1
config.Data.publication = False
config.Data.outputDatasetTag = 'QCD_PU200'

config.section_("Site")
config.Site.storageSite = 'T2_CN_Beijing'  #2_CH_CERN'