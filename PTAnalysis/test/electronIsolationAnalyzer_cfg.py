import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 10 )

inputFilesAOD = cms.untracked.vstring(
    # AOD test files from /DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM
    '/store/mc/PhaseIISpr18AODMiniAOD/DYToLL-M-50_2J_14TeV-madgraphMLM-pythia8/AODSIM/noPU_93X_upgrade2023_realistic_v5-v1/210000/F49D8A6B-B348-E811-9DC2-A4BF011254E0.root',
    )
inputFilesGenSimReco = cms.untracked.vstring(
    # GEN-SIM-RECO test files from /DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO
    '/store/mc/PhaseIITDRFall17DR/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/0051AD4B-31AF-E711-8324-7845C4FC3683.root',
    )
inputFilesMiniAOD = cms.untracked.vstring(
    # MiniAOD test files from /DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM
    '/store/mc/PhaseIISpr18AODMiniAOD/DYToLL-M-50_2J_14TeV-madgraphMLM-pythia8/MINIAODSIM/noPU_93X_upgrade2023_realistic_v5-v1/10000/A2B8867E-4944-E811-86F0-001E67792650.root',
    )

# Set up input/output depending on the format
# You can list here either AOD or miniAOD files, but not both types mixed
#
useAOD = False

if useAOD == True :
    inputFiles = inputFilesAOD
    outputFile = "electronIsolation.root"
    print("AOD input files are used")
else :
    inputFiles = inputFilesGenSimReco #inputFiles = inputFilesMiniAOD
    outputFile = "electronIsolation_GEN-SIM-RECO.root" #outputFile = "electronIsolationminiAOD.root"
    print("GEN-SIM-RECO input files are used")
process.source = cms.Source ("PoolSource", fileNames = inputFiles )

#
# Set up electron ID (VID framework)
#
'''from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
 DataFormat.AOD or DataFormat.MiniAOD, as appropriate
if useAOD == True :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)
#switchOnVIDElectronIdProducer(process)
# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff']
#my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff']
# add them to the VID producer
for idmod in my_id_modules:
   setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)'''

process.analysis = cms.EDAnalyzer(
    'ElectronIsolationAnalyzer',
    genXYZTag = cms.untracked.InputTag("genParticles", "xyz0", "HLT"),
    genT0Tag = cms.untracked.InputTag("genParticles", "t0", "HLT"),
    VertexTag3D  = cms.InputTag("offlinePrimaryVertices"),
    VertexTag4D  = cms.InputTag("offlinePrimaryVertices4D"),
    PileUpTag    = cms.InputTag("addPileupInfo"),
    barrelElectronsTag   = cms.untracked.InputTag("gedGsfElectrons"),
    #endcapElectronsTag   = cms.untracked.InputTag("ecalDrivenGsfElectronsFromMultiCl"),
    TracksTag    = cms.InputTag("generalTracks"),
    TrackTimeValueMapTag = cms.InputTag("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModel"),
    PFCandidateTag = cms.InputTag("particleFlow", "", "RECO"),
    genPartTag = cms.untracked.InputTag("genParticles", "", "HLT"),
    genVtxTag = cms.untracked.InputTag("g4SimHits", "", "SIM"),
    genJetsTag = cms.untracked.InputTag("ak4GenJets", "", "HLT"),
    timeResolutions = cms.untracked.vdouble(0.030, 0.050, 0.070),
    isoConeDR = cms.untracked.vdouble(0.2, 0.3, 0.4, 0.5),
    saveTracks = cms.untracked.bool(True),
    maxDz = cms.untracked.double(0.1),
    minDr = cms.untracked.double(0.01),
    isAOD=cms.untracked.bool(True),
    #
    # ID decisions (common to all formats)
    #
    '''eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-veto"),
    eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-loose"),
    eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-medium"),
    eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-tight"),'''
    #eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV70")
)

# Output TFile
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outputFile))
# Make sure to add the ID sequence upstream from the user analysis module
process.p = cms.Path(process.analysis)
#process.p = cms.Path(process.egmGsfElectronIDSequence * process.analysis)
