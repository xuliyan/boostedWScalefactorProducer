
// THIS FILE HAS BEEN GENERATED AUTOMATICALLY. DO NOT EDIT DIRECTLY, CHANGES WILL BE LOST UPON NEXT CODE GENERATION.
// Code produced by Id: CodeIt.py 494 2010-07-30 13:41:32Z svn 

#ifndef __UZHTOP_Tau_H__
#define __UZHTOP_Tau_H__

#include <cmath>
#include "Particle.h"
#include <vector>
#include "TauNtupleObject.h"



namespace Ntuple {
  class TauNtupleObject;
}

namespace UZH {
  
  /**
   *  @short Class that maps TauNtupleObjects to Tau Particle class
   *
   *         This class can be used to map the offline Tau information from
   *         TauNtupleObjects to UZH::Tau class. All particles inherit 
   *         from UZH::Particle.
   *
   * @author Code produced by Id: CodeIt.py 494 2010-07-30 13:41:32Z svn 
   *
   */

  class Tau : public Basic 
 , public Particle 
  {
  public:

    /// default c'tor
    Tau();
    /// default d'tor
    ~Tau();
    
    /// c'tor with index
    Tau( const Ntuple::TauNtupleObject* ana, const Int_t idx );

    TLorentzVector* getTLV() const;
    TLorentzVector tlv() const;


    
    // variable definitions



    floatingnumber* m_decayModeFindingNewDMs;
    floatingnumber* m_decayModeFinding;
    floatingnumber* m_byLooseCombinedIsolationDeltaBetaCorr3Hits;
    floatingnumber* m_byMediumCombinedIsolationDeltaBetaCorr3Hits;
    floatingnumber* m_byTightCombinedIsolationDeltaBetaCorr3Hits;
    floatingnumber* m_byCombinedIsolationDeltaBetaCorrRaw3Hits;
    floatingnumber* m_chargedIsoPtSum;
    floatingnumber* m_neutralIsoPtSum;
    floatingnumber* m_puCorrPtSum;
    floatingnumber* m_byIsolationMVA3oldDMwLTraw;
    floatingnumber* m_byVLooseIsolationMVA3oldDMwLT;
    floatingnumber* m_byLooseIsolationMVA3oldDMwLT;
    floatingnumber* m_byMediumIsolationMVA3oldDMwLT;
    floatingnumber* m_byTightIsolationMVA3oldDMwLT;
    floatingnumber* m_byVTightIsolationMVA3oldDMwLT;
    floatingnumber* m_byIsolationMVA3newDMwLTraw;
    floatingnumber* m_byVLooseIsolationMVA3newDMwLT;
    floatingnumber* m_byLooseIsolationMVA3newDMwLT;
    floatingnumber* m_byMediumIsolationMVA3newDMwLT;
    floatingnumber* m_byTightIsolationMVA3newDMwLT;
    floatingnumber* m_byVTightIsolationMVA3newDMwLT;
    floatingnumber* m_againstElectronMVA5raw;
    floatingnumber* m_againstElectronMVA5category;
    floatingnumber* m_againstElectronVLooseMVA5;
    floatingnumber* m_againstElectronLooseMVA5;
    floatingnumber* m_againstElectronMediumMVA5;
    floatingnumber* m_againstElectronTightMVA5;
    floatingnumber* m_againstElectronVTightMVA5;
    floatingnumber* m_againstMuonLoose3;
    floatingnumber* m_againstMuonTight3;
    floatingnumber* m_byPileupWeightedIsolationRaw3Hits;
    floatingnumber* m_byLoosePileupWeightedIsolation3Hits;
    floatingnumber* m_byMediumPileupWeightedIsolation3Hits;
    floatingnumber* m_byTightPileupWeightedIsolation3Hits;
    floatingnumber* m_byPhotonPtSumOutsideSignalCone;
    floatingnumber* m_footprintCorrection;
    int* m_pdgId;
    floatingnumber* m_charge;
    floatingnumber* m_d0;
    int* m_TauType;







    // check level given here must be consistent with ...NtupleObject.cxx, otherwise you'll get a segfault
    floatingnumber decayModeFindingNewDMs() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kdecayModeFindingNewDMs)) std::cout<<"decayModeFindingNewDMs not connected!"<<std::endl;*/ return *(m_decayModeFindingNewDMs); } 
    floatingnumber decayModeFinding() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kdecayModeFinding)) std::cout<<"decayModeFinding not connected!"<<std::endl;*/ return *(m_decayModeFinding); } 
    floatingnumber byLooseCombinedIsolationDeltaBetaCorr3Hits() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kbyLooseCombinedIsolationDeltaBetaCorr3Hits)) std::cout<<"byLooseCombinedIsolationDeltaBetaCorr3Hits not connected!"<<std::endl;*/ return *(m_byLooseCombinedIsolationDeltaBetaCorr3Hits); } 
    floatingnumber byMediumCombinedIsolationDeltaBetaCorr3Hits() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kbyMediumCombinedIsolationDeltaBetaCorr3Hits)) std::cout<<"byMediumCombinedIsolationDeltaBetaCorr3Hits not connected!"<<std::endl;*/ return *(m_byMediumCombinedIsolationDeltaBetaCorr3Hits); } 
    floatingnumber byTightCombinedIsolationDeltaBetaCorr3Hits() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kbyTightCombinedIsolationDeltaBetaCorr3Hits)) std::cout<<"byTightCombinedIsolationDeltaBetaCorr3Hits not connected!"<<std::endl;*/ return *(m_byTightCombinedIsolationDeltaBetaCorr3Hits); } 
    floatingnumber byCombinedIsolationDeltaBetaCorrRaw3Hits() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kbyCombinedIsolationDeltaBetaCorrRaw3Hits)) std::cout<<"byCombinedIsolationDeltaBetaCorrRaw3Hits not connected!"<<std::endl;*/ return *(m_byCombinedIsolationDeltaBetaCorrRaw3Hits); } 
    floatingnumber chargedIsoPtSum() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kchargedIsoPtSum)) std::cout<<"chargedIsoPtSum not connected!"<<std::endl;*/ return *(m_chargedIsoPtSum); } 
    floatingnumber neutralIsoPtSum() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kneutralIsoPtSum)) std::cout<<"neutralIsoPtSum not connected!"<<std::endl;*/ return *(m_neutralIsoPtSum); } 
    floatingnumber puCorrPtSum() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kpuCorrPtSum)) std::cout<<"puCorrPtSum not connected!"<<std::endl;*/ return *(m_puCorrPtSum); } 
    floatingnumber byIsolationMVA3oldDMwLTraw() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kbyIsolationMVA3oldDMwLTraw)) std::cout<<"byIsolationMVA3oldDMwLTraw not connected!"<<std::endl;*/ return *(m_byIsolationMVA3oldDMwLTraw); } 
    floatingnumber byVLooseIsolationMVA3oldDMwLT() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kbyVLooseIsolationMVA3oldDMwLT)) std::cout<<"byVLooseIsolationMVA3oldDMwLT not connected!"<<std::endl;*/ return *(m_byVLooseIsolationMVA3oldDMwLT); } 
    floatingnumber byLooseIsolationMVA3oldDMwLT() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kbyLooseIsolationMVA3oldDMwLT)) std::cout<<"byLooseIsolationMVA3oldDMwLT not connected!"<<std::endl;*/ return *(m_byLooseIsolationMVA3oldDMwLT); } 
    floatingnumber byMediumIsolationMVA3oldDMwLT() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kbyMediumIsolationMVA3oldDMwLT)) std::cout<<"byMediumIsolationMVA3oldDMwLT not connected!"<<std::endl;*/ return *(m_byMediumIsolationMVA3oldDMwLT); } 
    floatingnumber byTightIsolationMVA3oldDMwLT() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kbyTightIsolationMVA3oldDMwLT)) std::cout<<"byTightIsolationMVA3oldDMwLT not connected!"<<std::endl;*/ return *(m_byTightIsolationMVA3oldDMwLT); } 
    floatingnumber byVTightIsolationMVA3oldDMwLT() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kbyVTightIsolationMVA3oldDMwLT)) std::cout<<"byVTightIsolationMVA3oldDMwLT not connected!"<<std::endl;*/ return *(m_byVTightIsolationMVA3oldDMwLT); } 
    floatingnumber byIsolationMVA3newDMwLTraw() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kbyIsolationMVA3newDMwLTraw)) std::cout<<"byIsolationMVA3newDMwLTraw not connected!"<<std::endl;*/ return *(m_byIsolationMVA3newDMwLTraw); } 
    floatingnumber byVLooseIsolationMVA3newDMwLT() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kbyVLooseIsolationMVA3newDMwLT)) std::cout<<"byVLooseIsolationMVA3newDMwLT not connected!"<<std::endl;*/ return *(m_byVLooseIsolationMVA3newDMwLT); } 
    floatingnumber byLooseIsolationMVA3newDMwLT() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kbyLooseIsolationMVA3newDMwLT)) std::cout<<"byLooseIsolationMVA3newDMwLT not connected!"<<std::endl;*/ return *(m_byLooseIsolationMVA3newDMwLT); } 
    floatingnumber byMediumIsolationMVA3newDMwLT() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kbyMediumIsolationMVA3newDMwLT)) std::cout<<"byMediumIsolationMVA3newDMwLT not connected!"<<std::endl;*/ return *(m_byMediumIsolationMVA3newDMwLT); } 
    floatingnumber byTightIsolationMVA3newDMwLT() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kbyTightIsolationMVA3newDMwLT)) std::cout<<"byTightIsolationMVA3newDMwLT not connected!"<<std::endl;*/ return *(m_byTightIsolationMVA3newDMwLT); } 
    floatingnumber byVTightIsolationMVA3newDMwLT() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kbyVTightIsolationMVA3newDMwLT)) std::cout<<"byVTightIsolationMVA3newDMwLT not connected!"<<std::endl;*/ return *(m_byVTightIsolationMVA3newDMwLT); } 
    floatingnumber againstElectronMVA5raw() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kagainstElectronMVA5raw)) std::cout<<"againstElectronMVA5raw not connected!"<<std::endl;*/ return *(m_againstElectronMVA5raw); } 
    floatingnumber againstElectronMVA5category() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kagainstElectronMVA5category)) std::cout<<"againstElectronMVA5category not connected!"<<std::endl;*/ return *(m_againstElectronMVA5category); } 
    floatingnumber againstElectronVLooseMVA5() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kagainstElectronVLooseMVA5)) std::cout<<"againstElectronVLooseMVA5 not connected!"<<std::endl;*/ return *(m_againstElectronVLooseMVA5); } 
    floatingnumber againstElectronLooseMVA5() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kagainstElectronLooseMVA5)) std::cout<<"againstElectronLooseMVA5 not connected!"<<std::endl;*/ return *(m_againstElectronLooseMVA5); } 
    floatingnumber againstElectronMediumMVA5() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kagainstElectronMediumMVA5)) std::cout<<"againstElectronMediumMVA5 not connected!"<<std::endl;*/ return *(m_againstElectronMediumMVA5); } 
    floatingnumber againstElectronTightMVA5() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kagainstElectronTightMVA5)) std::cout<<"againstElectronTightMVA5 not connected!"<<std::endl;*/ return *(m_againstElectronTightMVA5); } 
    floatingnumber againstElectronVTightMVA5() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kagainstElectronVTightMVA5)) std::cout<<"againstElectronVTightMVA5 not connected!"<<std::endl;*/ return *(m_againstElectronVTightMVA5); } 
    floatingnumber againstMuonLoose3() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kagainstMuonLoose3)) std::cout<<"againstMuonLoose3 not connected!"<<std::endl;*/ return *(m_againstMuonLoose3); } 
    floatingnumber againstMuonTight3() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kagainstMuonTight3)) std::cout<<"againstMuonTight3 not connected!"<<std::endl;*/ return *(m_againstMuonTight3); } 
    floatingnumber byPileupWeightedIsolationRaw3Hits() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kbyPileupWeightedIsolationRaw3Hits)) std::cout<<"byPileupWeightedIsolationRaw3Hits not connected!"<<std::endl;*/ return *(m_byPileupWeightedIsolationRaw3Hits); } 
    floatingnumber byLoosePileupWeightedIsolation3Hits() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kbyLoosePileupWeightedIsolation3Hits)) std::cout<<"byLoosePileupWeightedIsolation3Hits not connected!"<<std::endl;*/ return *(m_byLoosePileupWeightedIsolation3Hits); } 
    floatingnumber byMediumPileupWeightedIsolation3Hits() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kbyMediumPileupWeightedIsolation3Hits)) std::cout<<"byMediumPileupWeightedIsolation3Hits not connected!"<<std::endl;*/ return *(m_byMediumPileupWeightedIsolation3Hits); } 
    floatingnumber byTightPileupWeightedIsolation3Hits() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kbyTightPileupWeightedIsolation3Hits)) std::cout<<"byTightPileupWeightedIsolation3Hits not connected!"<<std::endl;*/ return *(m_byTightPileupWeightedIsolation3Hits); } 
    floatingnumber byPhotonPtSumOutsideSignalCone() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kbyPhotonPtSumOutsideSignalCone)) std::cout<<"byPhotonPtSumOutsideSignalCone not connected!"<<std::endl;*/ return *(m_byPhotonPtSumOutsideSignalCone); } 
    floatingnumber footprintCorrection() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kfootprintCorrection)) std::cout<<"footprintCorrection not connected!"<<std::endl;*/ return *(m_footprintCorrection); } 
    int pdgId() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kpdgId)) std::cout<<"pdgId not connected!"<<std::endl;*/ return *(m_pdgId); } 
    floatingnumber charge() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kcharge)) std::cout<<"charge not connected!"<<std::endl;*/ return *(m_charge); } 
    floatingnumber d0() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kd0)) std::cout<<"d0 not connected!"<<std::endl;*/ return *(m_d0); } 
    int TauType() const { /*if(!m_ana->getConnectSucceeded(Ntuple::TauNtupleObject::kTauType)) std::cout<<"TauType not connected!"<<std::endl;*/ return *(m_TauType); } 
    
    void decayModeFindingNewDMs( const floatingnumber& val){ *(m_decayModeFindingNewDMs)=val; } 
    void decayModeFinding( const floatingnumber& val){ *(m_decayModeFinding)=val; } 
    void byLooseCombinedIsolationDeltaBetaCorr3Hits( const floatingnumber& val){ *(m_byLooseCombinedIsolationDeltaBetaCorr3Hits)=val; } 
    void byMediumCombinedIsolationDeltaBetaCorr3Hits( const floatingnumber& val){ *(m_byMediumCombinedIsolationDeltaBetaCorr3Hits)=val; } 
    void byTightCombinedIsolationDeltaBetaCorr3Hits( const floatingnumber& val){ *(m_byTightCombinedIsolationDeltaBetaCorr3Hits)=val; } 
    void byCombinedIsolationDeltaBetaCorrRaw3Hits( const floatingnumber& val){ *(m_byCombinedIsolationDeltaBetaCorrRaw3Hits)=val; } 
    void chargedIsoPtSum( const floatingnumber& val){ *(m_chargedIsoPtSum)=val; } 
    void neutralIsoPtSum( const floatingnumber& val){ *(m_neutralIsoPtSum)=val; } 
    void puCorrPtSum( const floatingnumber& val){ *(m_puCorrPtSum)=val; } 
    void byIsolationMVA3oldDMwLTraw( const floatingnumber& val){ *(m_byIsolationMVA3oldDMwLTraw)=val; } 
    void byVLooseIsolationMVA3oldDMwLT( const floatingnumber& val){ *(m_byVLooseIsolationMVA3oldDMwLT)=val; } 
    void byLooseIsolationMVA3oldDMwLT( const floatingnumber& val){ *(m_byLooseIsolationMVA3oldDMwLT)=val; } 
    void byMediumIsolationMVA3oldDMwLT( const floatingnumber& val){ *(m_byMediumIsolationMVA3oldDMwLT)=val; } 
    void byTightIsolationMVA3oldDMwLT( const floatingnumber& val){ *(m_byTightIsolationMVA3oldDMwLT)=val; } 
    void byVTightIsolationMVA3oldDMwLT( const floatingnumber& val){ *(m_byVTightIsolationMVA3oldDMwLT)=val; } 
    void byIsolationMVA3newDMwLTraw( const floatingnumber& val){ *(m_byIsolationMVA3newDMwLTraw)=val; } 
    void byVLooseIsolationMVA3newDMwLT( const floatingnumber& val){ *(m_byVLooseIsolationMVA3newDMwLT)=val; } 
    void byLooseIsolationMVA3newDMwLT( const floatingnumber& val){ *(m_byLooseIsolationMVA3newDMwLT)=val; } 
    void byMediumIsolationMVA3newDMwLT( const floatingnumber& val){ *(m_byMediumIsolationMVA3newDMwLT)=val; } 
    void byTightIsolationMVA3newDMwLT( const floatingnumber& val){ *(m_byTightIsolationMVA3newDMwLT)=val; } 
    void byVTightIsolationMVA3newDMwLT( const floatingnumber& val){ *(m_byVTightIsolationMVA3newDMwLT)=val; } 
    void againstElectronMVA5raw( const floatingnumber& val){ *(m_againstElectronMVA5raw)=val; } 
    void againstElectronMVA5category( const floatingnumber& val){ *(m_againstElectronMVA5category)=val; } 
    void againstElectronVLooseMVA5( const floatingnumber& val){ *(m_againstElectronVLooseMVA5)=val; } 
    void againstElectronLooseMVA5( const floatingnumber& val){ *(m_againstElectronLooseMVA5)=val; } 
    void againstElectronMediumMVA5( const floatingnumber& val){ *(m_againstElectronMediumMVA5)=val; } 
    void againstElectronTightMVA5( const floatingnumber& val){ *(m_againstElectronTightMVA5)=val; } 
    void againstElectronVTightMVA5( const floatingnumber& val){ *(m_againstElectronVTightMVA5)=val; } 
    void againstMuonLoose3( const floatingnumber& val){ *(m_againstMuonLoose3)=val; } 
    void againstMuonTight3( const floatingnumber& val){ *(m_againstMuonTight3)=val; } 
    void byPileupWeightedIsolationRaw3Hits( const floatingnumber& val){ *(m_byPileupWeightedIsolationRaw3Hits)=val; } 
    void byLoosePileupWeightedIsolation3Hits( const floatingnumber& val){ *(m_byLoosePileupWeightedIsolation3Hits)=val; } 
    void byMediumPileupWeightedIsolation3Hits( const floatingnumber& val){ *(m_byMediumPileupWeightedIsolation3Hits)=val; } 
    void byTightPileupWeightedIsolation3Hits( const floatingnumber& val){ *(m_byTightPileupWeightedIsolation3Hits)=val; } 
    void byPhotonPtSumOutsideSignalCone( const floatingnumber& val){ *(m_byPhotonPtSumOutsideSignalCone)=val; } 
    void footprintCorrection( const floatingnumber& val){ *(m_footprintCorrection)=val; } 
    void pdgId( const int& val){ *(m_pdgId)=val; } 
    void charge( const floatingnumber& val){ *(m_charge)=val; } 
    void d0( const floatingnumber& val){ *(m_d0)=val; } 
    void TauType( const int& val){ *(m_TauType)=val; } 
    

  private:
    const Ntuple::TauNtupleObject* m_ana;
  }; // class Tau

  typedef std::vector< Tau > TauVec;
  typedef std::vector< Tau >::iterator TauVecIt;
  typedef std::vector< Tau >::const_iterator TauVecConstIt;



  /// sort Taus by pT
  bool operator<( const Tau& e1, const Tau& e2 );

  /// function class to sort Tau vector contents by pT
  class sortTauPt {
  public:
    bool operator()( const Tau& e1,
                     const Tau& e2 );
  };



} // end of namespace UZH

/// output stream operator overloaded for Tau objects
std::ostream& operator<<( std::ostream& out,
                          const UZH::Tau& rhs );


#endif //__UZH_Tau_H__
