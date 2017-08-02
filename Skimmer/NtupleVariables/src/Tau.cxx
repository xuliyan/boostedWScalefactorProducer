
// THIS FILE HAS BEEN GENERATED AUTOMATICALLY. DO NOT EDIT DIRECTLY, CHANGES WILL BE LOST UPON NEXT CODE GENERATION.
// Code produced by Id: CodeIt.py 494 2010-07-30 13:41:32Z svn 

#include "../include/Tau.h"

using namespace std;
using namespace UZH;

Tau::Tau() {
  
}


Tau::Tau( const Ntuple::TauNtupleObject* ana, const Int_t idx ) 
: Basic( idx )
 , Particle() 


{
  m_ana=ana;
  // copy variables defined in Particle.h
    m_e = &((*ana->e)[idx]); 
    m_pt = &((*ana->pt)[idx]); 
    m_eta = &((*ana->eta)[idx]); 
    m_phi = &((*ana->phi)[idx]); 
    m_m = &((*ana->m)[idx]); 

  m_lvl    = ana->detailLevel;


  // copy rest of variables
  /*${ {AllNoBools:    printf("acc#name#\n"); if (ana->m_connectsucceeded[#index#]) 
         {printf("?\n"); m_#name# = &((*ana->#name#)[idx]);}
    else {printf("fak\n"); m_#name# = new #type#(); *m_#name# = #default#; } }}
  */
if(  ((ana->detailLevel & Ntuple::TauAdvancedID) == Ntuple::TauAdvancedID)  ) {
     if (ana->m_connectsucceeded[5]) m_decayModeFindingNewDMs = &((*ana->decayModeFindingNewDMs)[idx]); else m_decayModeFindingNewDMs = 0; 
    if (ana->m_connectsucceeded[6]) m_decayModeFinding = &((*ana->decayModeFinding)[idx]); else m_decayModeFinding = 0; 
    if (ana->m_connectsucceeded[7]) m_byLooseCombinedIsolationDeltaBetaCorr3Hits = &((*ana->byLooseCombinedIsolationDeltaBetaCorr3Hits)[idx]); else m_byLooseCombinedIsolationDeltaBetaCorr3Hits = 0; 
    if (ana->m_connectsucceeded[8]) m_byMediumCombinedIsolationDeltaBetaCorr3Hits = &((*ana->byMediumCombinedIsolationDeltaBetaCorr3Hits)[idx]); else m_byMediumCombinedIsolationDeltaBetaCorr3Hits = 0; 
    if (ana->m_connectsucceeded[9]) m_byTightCombinedIsolationDeltaBetaCorr3Hits = &((*ana->byTightCombinedIsolationDeltaBetaCorr3Hits)[idx]); else m_byTightCombinedIsolationDeltaBetaCorr3Hits = 0; 
    if (ana->m_connectsucceeded[10]) m_byCombinedIsolationDeltaBetaCorrRaw3Hits = &((*ana->byCombinedIsolationDeltaBetaCorrRaw3Hits)[idx]); else m_byCombinedIsolationDeltaBetaCorrRaw3Hits = 0; 
    if (ana->m_connectsucceeded[11]) m_chargedIsoPtSum = &((*ana->chargedIsoPtSum)[idx]); else m_chargedIsoPtSum = 0; 
    if (ana->m_connectsucceeded[12]) m_neutralIsoPtSum = &((*ana->neutralIsoPtSum)[idx]); else m_neutralIsoPtSum = 0; 
    if (ana->m_connectsucceeded[13]) m_puCorrPtSum = &((*ana->puCorrPtSum)[idx]); else m_puCorrPtSum = 0; 
    if (ana->m_connectsucceeded[14]) m_byIsolationMVA3oldDMwLTraw = &((*ana->byIsolationMVA3oldDMwLTraw)[idx]); else m_byIsolationMVA3oldDMwLTraw = 0; 
    if (ana->m_connectsucceeded[15]) m_byVLooseIsolationMVA3oldDMwLT = &((*ana->byVLooseIsolationMVA3oldDMwLT)[idx]); else m_byVLooseIsolationMVA3oldDMwLT = 0; 
    if (ana->m_connectsucceeded[16]) m_byLooseIsolationMVA3oldDMwLT = &((*ana->byLooseIsolationMVA3oldDMwLT)[idx]); else m_byLooseIsolationMVA3oldDMwLT = 0; 
    if (ana->m_connectsucceeded[17]) m_byMediumIsolationMVA3oldDMwLT = &((*ana->byMediumIsolationMVA3oldDMwLT)[idx]); else m_byMediumIsolationMVA3oldDMwLT = 0; 
    if (ana->m_connectsucceeded[18]) m_byTightIsolationMVA3oldDMwLT = &((*ana->byTightIsolationMVA3oldDMwLT)[idx]); else m_byTightIsolationMVA3oldDMwLT = 0; 
    if (ana->m_connectsucceeded[19]) m_byVTightIsolationMVA3oldDMwLT = &((*ana->byVTightIsolationMVA3oldDMwLT)[idx]); else m_byVTightIsolationMVA3oldDMwLT = 0; 
    if (ana->m_connectsucceeded[20]) m_byIsolationMVA3newDMwLTraw = &((*ana->byIsolationMVA3newDMwLTraw)[idx]); else m_byIsolationMVA3newDMwLTraw = 0; 
    if (ana->m_connectsucceeded[21]) m_byVLooseIsolationMVA3newDMwLT = &((*ana->byVLooseIsolationMVA3newDMwLT)[idx]); else m_byVLooseIsolationMVA3newDMwLT = 0; 
    if (ana->m_connectsucceeded[22]) m_byLooseIsolationMVA3newDMwLT = &((*ana->byLooseIsolationMVA3newDMwLT)[idx]); else m_byLooseIsolationMVA3newDMwLT = 0; 
    if (ana->m_connectsucceeded[23]) m_byMediumIsolationMVA3newDMwLT = &((*ana->byMediumIsolationMVA3newDMwLT)[idx]); else m_byMediumIsolationMVA3newDMwLT = 0; 
    if (ana->m_connectsucceeded[24]) m_byTightIsolationMVA3newDMwLT = &((*ana->byTightIsolationMVA3newDMwLT)[idx]); else m_byTightIsolationMVA3newDMwLT = 0; 
    if (ana->m_connectsucceeded[25]) m_byVTightIsolationMVA3newDMwLT = &((*ana->byVTightIsolationMVA3newDMwLT)[idx]); else m_byVTightIsolationMVA3newDMwLT = 0; 
    if (ana->m_connectsucceeded[26]) m_againstElectronMVA5raw = &((*ana->againstElectronMVA5raw)[idx]); else m_againstElectronMVA5raw = 0; 
    if (ana->m_connectsucceeded[27]) m_againstElectronMVA5category = &((*ana->againstElectronMVA5category)[idx]); else m_againstElectronMVA5category = 0; 
    if (ana->m_connectsucceeded[28]) m_againstElectronVLooseMVA5 = &((*ana->againstElectronVLooseMVA5)[idx]); else m_againstElectronVLooseMVA5 = 0; 
    if (ana->m_connectsucceeded[29]) m_againstElectronLooseMVA5 = &((*ana->againstElectronLooseMVA5)[idx]); else m_againstElectronLooseMVA5 = 0; 
    if (ana->m_connectsucceeded[30]) m_againstElectronMediumMVA5 = &((*ana->againstElectronMediumMVA5)[idx]); else m_againstElectronMediumMVA5 = 0; 
    if (ana->m_connectsucceeded[31]) m_againstElectronTightMVA5 = &((*ana->againstElectronTightMVA5)[idx]); else m_againstElectronTightMVA5 = 0; 
    if (ana->m_connectsucceeded[32]) m_againstElectronVTightMVA5 = &((*ana->againstElectronVTightMVA5)[idx]); else m_againstElectronVTightMVA5 = 0; 
    if (ana->m_connectsucceeded[33]) m_againstMuonLoose3 = &((*ana->againstMuonLoose3)[idx]); else m_againstMuonLoose3 = 0; 
    if (ana->m_connectsucceeded[34]) m_againstMuonTight3 = &((*ana->againstMuonTight3)[idx]); else m_againstMuonTight3 = 0; 
    if (ana->m_connectsucceeded[35]) m_byPileupWeightedIsolationRaw3Hits = &((*ana->byPileupWeightedIsolationRaw3Hits)[idx]); else m_byPileupWeightedIsolationRaw3Hits = 0; 
    if (ana->m_connectsucceeded[36]) m_byLoosePileupWeightedIsolation3Hits = &((*ana->byLoosePileupWeightedIsolation3Hits)[idx]); else m_byLoosePileupWeightedIsolation3Hits = 0; 
    if (ana->m_connectsucceeded[37]) m_byMediumPileupWeightedIsolation3Hits = &((*ana->byMediumPileupWeightedIsolation3Hits)[idx]); else m_byMediumPileupWeightedIsolation3Hits = 0; 
    if (ana->m_connectsucceeded[38]) m_byTightPileupWeightedIsolation3Hits = &((*ana->byTightPileupWeightedIsolation3Hits)[idx]); else m_byTightPileupWeightedIsolation3Hits = 0; 
    if (ana->m_connectsucceeded[39]) m_byPhotonPtSumOutsideSignalCone = &((*ana->byPhotonPtSumOutsideSignalCone)[idx]); else m_byPhotonPtSumOutsideSignalCone = 0; 
    if (ana->m_connectsucceeded[40]) m_footprintCorrection = &((*ana->footprintCorrection)[idx]); else m_footprintCorrection = 0; 
} // end of detail level AdvancedID

if(  ((ana->detailLevel & Ntuple::TauBasic) == Ntuple::TauBasic)  ) {
     if (ana->m_connectsucceeded[1]) m_pdgId = &((*ana->pdgId)[idx]); else m_pdgId = 0; 
    if (ana->m_connectsucceeded[2]) m_charge = &((*ana->charge)[idx]); else m_charge = 0; 
    if (ana->m_connectsucceeded[3]) m_d0 = &((*ana->d0)[idx]); else m_d0 = 0; 
} // end of detail level Basic

if(  ((ana->detailLevel & Ntuple::TauID) == Ntuple::TauID)  ) {
     if (ana->m_connectsucceeded[4]) m_TauType = &((*ana->TauType)[idx]); else m_TauType = 0; 
}







}


Tau::~Tau() {

}

ostream& operator<<( ostream& out,
                     const Tau& rhs ) {
  
   out << "Tau -" << ( Basic) rhs; 



  ;
if(  ((rhs.getLvl() & Ntuple::TauAdvancedID) == Ntuple::TauAdvancedID)  ) {
   out << " decayModeFindingNewDMs " << rhs.decayModeFindingNewDMs();
  out << " decayModeFinding " << rhs.decayModeFinding();
  out << " byLooseCombinedIsolationDeltaBetaCorr3Hits " << rhs.byLooseCombinedIsolationDeltaBetaCorr3Hits();
  out << " byMediumCombinedIsolationDeltaBetaCorr3Hits " << rhs.byMediumCombinedIsolationDeltaBetaCorr3Hits();
  out << " byTightCombinedIsolationDeltaBetaCorr3Hits " << rhs.byTightCombinedIsolationDeltaBetaCorr3Hits();
  out << " byCombinedIsolationDeltaBetaCorrRaw3Hits " << rhs.byCombinedIsolationDeltaBetaCorrRaw3Hits();
  out << " chargedIsoPtSum " << rhs.chargedIsoPtSum();
  out << " neutralIsoPtSum " << rhs.neutralIsoPtSum();
  out << " puCorrPtSum " << rhs.puCorrPtSum();
  out << " byIsolationMVA3oldDMwLTraw " << rhs.byIsolationMVA3oldDMwLTraw();
  out << " byVLooseIsolationMVA3oldDMwLT " << rhs.byVLooseIsolationMVA3oldDMwLT();
  out << " byLooseIsolationMVA3oldDMwLT " << rhs.byLooseIsolationMVA3oldDMwLT();
  out << " byMediumIsolationMVA3oldDMwLT " << rhs.byMediumIsolationMVA3oldDMwLT();
  out << " byTightIsolationMVA3oldDMwLT " << rhs.byTightIsolationMVA3oldDMwLT();
  out << " byVTightIsolationMVA3oldDMwLT " << rhs.byVTightIsolationMVA3oldDMwLT();
  out << " byIsolationMVA3newDMwLTraw " << rhs.byIsolationMVA3newDMwLTraw();
  out << " byVLooseIsolationMVA3newDMwLT " << rhs.byVLooseIsolationMVA3newDMwLT();
  out << " byLooseIsolationMVA3newDMwLT " << rhs.byLooseIsolationMVA3newDMwLT();
  out << " byMediumIsolationMVA3newDMwLT " << rhs.byMediumIsolationMVA3newDMwLT();
  out << " byTightIsolationMVA3newDMwLT " << rhs.byTightIsolationMVA3newDMwLT();
  out << " byVTightIsolationMVA3newDMwLT " << rhs.byVTightIsolationMVA3newDMwLT();
  out << " againstElectronMVA5raw " << rhs.againstElectronMVA5raw();
  out << " againstElectronMVA5category " << rhs.againstElectronMVA5category();
  out << " againstElectronVLooseMVA5 " << rhs.againstElectronVLooseMVA5();
  out << " againstElectronLooseMVA5 " << rhs.againstElectronLooseMVA5();
  out << " againstElectronMediumMVA5 " << rhs.againstElectronMediumMVA5();
  out << " againstElectronTightMVA5 " << rhs.againstElectronTightMVA5();
  out << " againstElectronVTightMVA5 " << rhs.againstElectronVTightMVA5();
  out << " againstMuonLoose3 " << rhs.againstMuonLoose3();
  out << " againstMuonTight3 " << rhs.againstMuonTight3();
  out << " byPileupWeightedIsolationRaw3Hits " << rhs.byPileupWeightedIsolationRaw3Hits();
  out << " byLoosePileupWeightedIsolation3Hits " << rhs.byLoosePileupWeightedIsolation3Hits();
  out << " byMediumPileupWeightedIsolation3Hits " << rhs.byMediumPileupWeightedIsolation3Hits();
  out << " byTightPileupWeightedIsolation3Hits " << rhs.byTightPileupWeightedIsolation3Hits();
  out << " byPhotonPtSumOutsideSignalCone " << rhs.byPhotonPtSumOutsideSignalCone();
  out << " footprintCorrection " << rhs.footprintCorrection();
;
} // end of detail level AdvancedID

if(  ((rhs.getLvl() & Ntuple::TauBasic) == Ntuple::TauBasic)  ) {
   out << " pdgId " << rhs.pdgId();
  out << " charge " << rhs.charge();
  out << " d0 " << rhs.d0();
;
} // end of detail level Basic

if(  ((rhs.getLvl() & Ntuple::TauID) == Ntuple::TauID)  ) {
   out << " TauType " << rhs.TauType();
;
}


  return out;
}



bool sortTauPt::operator()( const Tau& e1,
                                 const Tau& e2 ) {
  return ( e1.pt() > e2.pt() ) ? true : false;
}

bool operator<<( const Tau& e1, const Tau& e2 ) {
  sortTauPt sort;
  return sort( e1, e2 );
}




TLorentzVector* Tau::getTLV() const {

  TLorentzVector* tlv = new TLorentzVector();
  tlv->SetPtEtaPhiE(*(m_pt), *(m_eta), *(m_phi), *(m_e));
  return tlv;

}


TLorentzVector Tau::tlv() const {

  TLorentzVector tlv;
  tlv.SetPtEtaPhiE(*(m_pt), *(m_eta), *(m_phi), *(m_e));
  return tlv;

}








