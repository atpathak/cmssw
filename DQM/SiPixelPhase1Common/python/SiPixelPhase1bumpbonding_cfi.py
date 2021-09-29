import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDHarvester import DQMEDHarvester

SiPixelPhase1bumpbonding = DQMEDHarvester("SiPixelPhase1bumpbonding",
    TopFolderName = cms.string('PixelPhase1/Phase1_MechanicalView'),
    MinHits = cms.int32(30)
)
