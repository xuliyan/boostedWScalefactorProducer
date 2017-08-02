
// THIS FILE HAS BEEN GENERATED AUTOMATICALLY. DO NOT EDIT DIRECTLY, CHANGES WILL BE LOST UPON NEXT CODE GENERATION.
// Code produced by Id: CodeIt.py 494 2010-07-30 13:41:32Z svn 

// Local include(s):
#include "../include/TauNtupleObject.h"

namespace Ntuple {

  TauNtupleObject::TauNtupleObject( SCycleBaseNTuple* parent )
    : SInputVariables< SCycleBaseNTuple >( parent ) {
      m_connectsucceeded.resize(kEnd);
  }

  void TauNtupleObject::setConnectSucceeded(const unsigned int index, const bool success) {
    if (m_connectsucceeded.size() < index+1)  m_connectsucceeded.resize(index+1);
    m_connectsucceeded.at(index) = success;
  }

  void TauNtupleObject::ConnectVariables( const TString& treeName,
                                         const TString& prefix,
                                         const TString& ntupleType ) throw( SError ) {
   
     TauNtupleObject::ConnectVariables( treeName, Ntuple::TauAll, prefix, ntupleType);
                                         
  }                                         

  void TauNtupleObject::ConnectVariables( const TString& treeName,
                                         UInt_t detail_level,
                                         const TString& prefix,
                                         const TString& ntupleType ) throw( SError ) {
                                         
    // get instance of NtupleObjectNames
    NtupleObjectNames m_objectNames(ntupleType);
    

    // particle vector size
    ConnectVariable( treeName, prefix + m_objectNames.getName("N"), N );
    
    // connect variables that are defined in Particle.cxx
    ConnectVariable( treeName, prefix + m_objectNames.getName("e"), e );
    ConnectVariable( treeName, prefix + m_objectNames.getName("pt"), pt );
    ConnectVariable( treeName, prefix + m_objectNames.getName("eta"), eta );
    ConnectVariable( treeName, prefix + m_objectNames.getName("phi"), phi );
    ConnectVariable( treeName, prefix + m_objectNames.getName("m"), m );
        
    

    // connect object specific variables
if(  ((detail_level & Ntuple::TauAdvancedID) == Ntuple::TauAdvancedID)  ) {
     setConnectSucceeded(5, ConnectVariable( treeName, prefix + m_objectNames.getName("decayModeFindingNewDMs"), decayModeFindingNewDMs)); 
    setConnectSucceeded(6, ConnectVariable( treeName, prefix + m_objectNames.getName("decayModeFinding"), decayModeFinding)); 
    setConnectSucceeded(7, ConnectVariable( treeName, prefix + m_objectNames.getName("byLooseCombinedIsolationDeltaBetaCorr3Hits"), byLooseCombinedIsolationDeltaBetaCorr3Hits)); 
    setConnectSucceeded(8, ConnectVariable( treeName, prefix + m_objectNames.getName("byMediumCombinedIsolationDeltaBetaCorr3Hits"), byMediumCombinedIsolationDeltaBetaCorr3Hits)); 
    setConnectSucceeded(9, ConnectVariable( treeName, prefix + m_objectNames.getName("byTightCombinedIsolationDeltaBetaCorr3Hits"), byTightCombinedIsolationDeltaBetaCorr3Hits)); 
    setConnectSucceeded(10, ConnectVariable( treeName, prefix + m_objectNames.getName("byCombinedIsolationDeltaBetaCorrRaw3Hits"), byCombinedIsolationDeltaBetaCorrRaw3Hits)); 
    setConnectSucceeded(11, ConnectVariable( treeName, prefix + m_objectNames.getName("chargedIsoPtSum"), chargedIsoPtSum)); 
    setConnectSucceeded(12, ConnectVariable( treeName, prefix + m_objectNames.getName("neutralIsoPtSum"), neutralIsoPtSum)); 
    setConnectSucceeded(13, ConnectVariable( treeName, prefix + m_objectNames.getName("puCorrPtSum"), puCorrPtSum)); 
    setConnectSucceeded(14, ConnectVariable( treeName, prefix + m_objectNames.getName("byIsolationMVA3oldDMwLTraw"), byIsolationMVA3oldDMwLTraw)); 
    setConnectSucceeded(15, ConnectVariable( treeName, prefix + m_objectNames.getName("byVLooseIsolationMVA3oldDMwLT"), byVLooseIsolationMVA3oldDMwLT)); 
    setConnectSucceeded(16, ConnectVariable( treeName, prefix + m_objectNames.getName("byLooseIsolationMVA3oldDMwLT"), byLooseIsolationMVA3oldDMwLT)); 
    setConnectSucceeded(17, ConnectVariable( treeName, prefix + m_objectNames.getName("byMediumIsolationMVA3oldDMwLT"), byMediumIsolationMVA3oldDMwLT)); 
    setConnectSucceeded(18, ConnectVariable( treeName, prefix + m_objectNames.getName("byTightIsolationMVA3oldDMwLT"), byTightIsolationMVA3oldDMwLT)); 
    setConnectSucceeded(19, ConnectVariable( treeName, prefix + m_objectNames.getName("byVTightIsolationMVA3oldDMwLT"), byVTightIsolationMVA3oldDMwLT)); 
    setConnectSucceeded(20, ConnectVariable( treeName, prefix + m_objectNames.getName("byIsolationMVA3newDMwLTraw"), byIsolationMVA3newDMwLTraw)); 
    setConnectSucceeded(21, ConnectVariable( treeName, prefix + m_objectNames.getName("byVLooseIsolationMVA3newDMwLT"), byVLooseIsolationMVA3newDMwLT)); 
    setConnectSucceeded(22, ConnectVariable( treeName, prefix + m_objectNames.getName("byLooseIsolationMVA3newDMwLT"), byLooseIsolationMVA3newDMwLT)); 
    setConnectSucceeded(23, ConnectVariable( treeName, prefix + m_objectNames.getName("byMediumIsolationMVA3newDMwLT"), byMediumIsolationMVA3newDMwLT)); 
    setConnectSucceeded(24, ConnectVariable( treeName, prefix + m_objectNames.getName("byTightIsolationMVA3newDMwLT"), byTightIsolationMVA3newDMwLT)); 
    setConnectSucceeded(25, ConnectVariable( treeName, prefix + m_objectNames.getName("byVTightIsolationMVA3newDMwLT"), byVTightIsolationMVA3newDMwLT)); 
    setConnectSucceeded(26, ConnectVariable( treeName, prefix + m_objectNames.getName("againstElectronMVA5raw"), againstElectronMVA5raw)); 
    setConnectSucceeded(27, ConnectVariable( treeName, prefix + m_objectNames.getName("againstElectronMVA5category"), againstElectronMVA5category)); 
    setConnectSucceeded(28, ConnectVariable( treeName, prefix + m_objectNames.getName("againstElectronVLooseMVA5"), againstElectronVLooseMVA5)); 
    setConnectSucceeded(29, ConnectVariable( treeName, prefix + m_objectNames.getName("againstElectronLooseMVA5"), againstElectronLooseMVA5)); 
    setConnectSucceeded(30, ConnectVariable( treeName, prefix + m_objectNames.getName("againstElectronMediumMVA5"), againstElectronMediumMVA5)); 
    setConnectSucceeded(31, ConnectVariable( treeName, prefix + m_objectNames.getName("againstElectronTightMVA5"), againstElectronTightMVA5)); 
    setConnectSucceeded(32, ConnectVariable( treeName, prefix + m_objectNames.getName("againstElectronVTightMVA5"), againstElectronVTightMVA5)); 
    setConnectSucceeded(33, ConnectVariable( treeName, prefix + m_objectNames.getName("againstMuonLoose3"), againstMuonLoose3)); 
    setConnectSucceeded(34, ConnectVariable( treeName, prefix + m_objectNames.getName("againstMuonTight3"), againstMuonTight3)); 
    setConnectSucceeded(35, ConnectVariable( treeName, prefix + m_objectNames.getName("byPileupWeightedIsolationRaw3Hits"), byPileupWeightedIsolationRaw3Hits)); 
    setConnectSucceeded(36, ConnectVariable( treeName, prefix + m_objectNames.getName("byLoosePileupWeightedIsolation3Hits"), byLoosePileupWeightedIsolation3Hits)); 
    setConnectSucceeded(37, ConnectVariable( treeName, prefix + m_objectNames.getName("byMediumPileupWeightedIsolation3Hits"), byMediumPileupWeightedIsolation3Hits)); 
    setConnectSucceeded(38, ConnectVariable( treeName, prefix + m_objectNames.getName("byTightPileupWeightedIsolation3Hits"), byTightPileupWeightedIsolation3Hits)); 
    setConnectSucceeded(39, ConnectVariable( treeName, prefix + m_objectNames.getName("byPhotonPtSumOutsideSignalCone"), byPhotonPtSumOutsideSignalCone)); 
    setConnectSucceeded(40, ConnectVariable( treeName, prefix + m_objectNames.getName("footprintCorrection"), footprintCorrection)); 
} // end of detail level AdvancedID

if(  ((detail_level & Ntuple::TauBasic) == Ntuple::TauBasic)  ) {
     setConnectSucceeded(1, ConnectVariable( treeName, prefix + m_objectNames.getName("pdgId"), pdgId)); 
    setConnectSucceeded(2, ConnectVariable( treeName, prefix + m_objectNames.getName("charge"), charge)); 
    setConnectSucceeded(3, ConnectVariable( treeName, prefix + m_objectNames.getName("d0"), d0)); 
} // end of detail level Basic

if(  ((detail_level & Ntuple::TauID) == Ntuple::TauID)  ) {
     setConnectSucceeded(4, ConnectVariable( treeName, prefix + m_objectNames.getName("TauType"), TauType)); 
}


        
    // save actual detail_level
    detailLevel = detail_level;

    return;

  }

} // namespace Ntuple
