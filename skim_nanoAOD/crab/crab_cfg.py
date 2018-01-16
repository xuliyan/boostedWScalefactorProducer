from WMCore.Configuration import Configuration

# this will use CRAB client API
from CRABAPI.RawCommand import crabCommand

# talk to DBS to get list of files in this dataset
from dbs.apis.dbsClient import DbsApi
dbs = DbsApi('https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader')

dataset0 = '/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/srappocc-TTJetsTuneCUETP8M113TeV-madgraphMLM-pythia8RunIISummer16MiniAODv2-PUMoriond1780XmcRun2-4a4b356339e753e24c281c17941d0081/USER'
fileDictList = dbs.listFiles(dataset=dataset0)

print ("dataset %s has %d files" % (dataset0, len(fileDictList)))

# DBS client returns a list of dictionaries, but we want a list of Logical File Names
lfnList = [ dic['logical_file_name'] for dic in fileDictList ]


config = Configuration()

config.section_("General")
config.General.requestName = 'NanoPost1'
config.General.transferLogs=True
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'crab_script.sh'
config.JobType.inputFiles = ['crab_script.py','haddnano.py'] #hadd nano will not be needed once nano tools are in cmssw
config.JobType.sendPythonFolder	 = True
config.section_("Data")
# config.Data.inputDataset = '/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/srappocc-TTJetsTuneCUETP8M113TeV-madgraphMLM-pythia8RunIISummer16MiniAODv2-PUMoriond1780XmcRun2-4a4b356339e753e24c281c17941d0081/USER'
# config.Data.inputDBS = 'global'
#config.Data.splitting = 'FileBased'

# following 3 lines are the trick to skip DBS data lookup in CRAB Server
config.Data.userInputFiles = lfnList
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1

# config.Data.splitting = 'EventAwareLumiBased'
# config.Data.unitsPerJob = 100
# config.Data.totalUnits = 2000
# config.Data.inputDBS='phys03'
config.Data.outLFNDirBase = '/store/user/thaarres/NanoTEST/'
config.Data.publication = False
config.Data.outputDatasetTag = 'NanoTestPost'
config.section_("Site")
# since there is no data discovery and no data location lookup in CRAB
# you have to say where the input files are
config.Site.whitelist = ['T3_US_FNALLPC','T3_CH_PSI']

config.Site.storageSite = "T3_CH_PSI"

result = crabCommand('submit', config = config)

print (result)

