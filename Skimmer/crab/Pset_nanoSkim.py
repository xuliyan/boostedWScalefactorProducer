
import FWCore.ParameterSet.Config as cms
process = cms.Process('NANO')
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(),

)
process.source.fileNames = [
        '/store/user/srappocc/SingleMuon/SingleMuon_Run2016G-07Aug17-v1/180113_045720/0000/test_data_80X_NANO_10.root'

]
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

process.options = cms.untracked.PSet()

process.output = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string('94XNanoV0-TTbar_SemiLep.root'), fakeNameForCrab =cms.untracked.bool(True))
process.out = cms.EndPath(process.output)
