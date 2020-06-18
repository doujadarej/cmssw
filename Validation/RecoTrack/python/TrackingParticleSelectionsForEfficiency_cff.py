import FWCore.ParameterSet.Config as cms

generalTpSelectorBlock = cms.PSet(
    lip = cms.double(30.0),
    chargedOnly = cms.bool(True),
    pdgId = cms.vint32(),
    signalOnly = cms.bool(True),
    intimeOnly = cms.bool(False),
    stableOnly = cms.bool(False),
    minRapidity = cms.double(-2.5),
    minHit = cms.int32(0),
    ptMin = cms.double(0.9),
    ptMax = cms.double(1e100),
    maxRapidity = cms.double(2.5),
    tip = cms.double(3.5),
    minPhi = cms.double(-3.2),
    maxPhi = cms.double(3.2),
)

from Configuration.Eras.Modifier_fastSim_cff import fastSim
fastSim.toModify(generalTpSelectorBlock, stableOnly = True)

TpSelectorForEfficiencyVsEtaBlock = generalTpSelectorBlock.clone()
TpSelectorForEfficiencyVsPhiBlock = generalTpSelectorBlock.clone()
TpSelectorForEfficiencyVsPtBlock = generalTpSelectorBlock.clone(ptMin = 0.050)
TpSelectorForEfficiencyVsVTXRBlock = generalTpSelectorBlock.clone(tip=1e5, lip=1e5) #removed cuts on tip and lip for DQM plots vs vertex r position to have all tracks (displaced)
TpSelectorForEfficiencyVsVTXZBlock = generalTpSelectorBlock.clone(tip=1e5, lip=1e5) #removed cuts on tip and lip for DQM plots vs vertex r position to have all tracks (displaced)

def _modifyForPhase1(pset):
    pset.minRapidity = -3
    pset.maxRapidity = 3
from Configuration.Eras.Modifier_phase1Pixel_cff import phase1Pixel
phase1Pixel.toModify(generalTpSelectorBlock, _modifyForPhase1)
phase1Pixel.toModify(TpSelectorForEfficiencyVsEtaBlock, _modifyForPhase1)
from Configuration.Eras.Modifier_phase2_tracker_cff import phase2_tracker
phase2_tracker.toModify(generalTpSelectorBlock, minRapidity=-4.5, maxRapidity=4.5)
phase2_tracker.toModify(TpSelectorForEfficiencyVsEtaBlock, minRapidity=-4.5, maxRapidity=4.5)
