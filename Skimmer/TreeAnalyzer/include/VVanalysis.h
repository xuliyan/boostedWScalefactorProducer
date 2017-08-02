// Dear emacs, this is -*- c++ -*-

#ifndef VVanalysis_H
#define VVanalysis_H

// STL include(s):
#include <vector>
#include <string>

// SFrame include(s):
#include "core/include/SCycleBase.h"
#include "plug-ins/include/SSummedVar.h"


// ROOT include(s):
#include <TBits.h>

// External include(s):
#include "../NtupleVariables/include/JetNtupleObject.h"
#include "../NtupleVariables/include/Jet.h"
#include "../NtupleVariables/include/EventInfoNtupleObject.h"
#include "../NtupleVariables/include/ElectronNtupleObject.h"
#include "../NtupleVariables/include/Electron.h"
#include "../NtupleVariables/include/MuonNtupleObject.h"
#include "../NtupleVariables/include/Muon.h"
#include "../NtupleVariables/include/MissingEtNtupleObject.h"
#include "../NtupleVariables/include/MissingEt.h"
#include "../NtupleVariables/include/GenParticleNtupleObject.h"
#include "../NtupleVariables/include/GenParticle.h"
// #include "../GoodRunsLists/include/TGoodRunsList.h"
#include "../PileupReweightingTool/include/PileupReweightingTool.h"
#include "../BTaggingTools/include/BTaggingScaleTool.h"

class TH1D;
class TH2D;
class TRandom3;
class TBits;
namespace UZH {
  class Jet;
  class Electron;
  class Muon;
  class MissingEt;
  class GenParticle;
}

/**
 *   @short Analyse flat UZH ntuple
 *
 *
 *  @author Thea Arrestad
 * @version $Revision: 1 $
 */
class VVanalysis : public SCycleBase {

public:
   /// Default constructor
   VVanalysis();
   /// Default destructor
   ~VVanalysis();

   /// Function called at the beginning/end of the cycle
   virtual void BeginCycle() throw( SError );
   virtual void EndCycle() throw( SError );
   
   /// Function called at the beginning of/finished processing a new input data
   virtual void BeginInputData( const SInputData& ) throw( SError );
   virtual void EndInputData  ( const SInputData& ) throw( SError );
   
   ///  Function called before/after processing each InputData block
   virtual void BeginMasterInputData( const SInputData& ) throw( SError );
   virtual void EndMasterInputData  ( const SInputData& ) throw( SError );
   


   /// Function called after opening each new input file
   virtual void BeginInputFile( const SInputData& ) throw( SError );

   /// Function called for every event
   virtual void ExecuteEvent( const SInputData&, Double_t ) throw( SError );

private:
   
  // Input variable objects:
  Ntuple::JetNtupleObject         m_jetAK4;       ///< jet container
  Ntuple::JetNtupleObject         m_jetAK8;       ///< jet container
  Ntuple::JetNtupleObject         m_jetAK8Puppi;  ///< jet container
  Ntuple::EventInfoNtupleObject   m_eventInfo;    ///< event info container
  Ntuple::ElectronNtupleObject    m_electron;     ///< electron container
  Ntuple::MuonNtupleObject        m_muon;         ///< muon container
  Ntuple::MissingEtNtupleObject   m_missingEt;    ///< missing E_T container
  Ntuple::GenParticleNtupleObject m_genParticle;  ///< gen particle container
  
  
  // Further objects
  // Root::TGoodRunsList m_grl;
  PileupReweightingTool m_pileupReweightingTool;
  BTaggingScaleTool m_bTaggingScaleTool;
  
  // Some counters:
  //
  SSummedVar< Int_t > m_allEvents; //!
  SSummedVar< Int_t > m_passedLoose; //!
  SSummedVar< Int_t > m_passedPuppi; //!
  SSummedVar< Int_t > m_passedMjj; //!
  SSummedVar< Int_t > m_passedDEta; //!
  SSummedVar< Int_t > m_passedEvents; //!
  SSummedVar< std::vector< Int_t > > m_test; //!
  
  Ntuple::JetNtupleObject         m_genjetAK8;    ///< jet container
  // The output variables
  
  float     m_o_mjj                       ; 
  float     m_o_genmjj                    ; 
  float     m_o_mpuppisoftdrop_jet1       ; 
  float     m_o_mpuppisoftdrop_jet2       ; 
  float     m_o_mgensoftdrop_jet1         ; 
  float     m_o_mgensoftdrop_jet2         ; 
  float     m_o_tau1_jet1                 ; 
  float     m_o_tau1_jet2                 ; 
  float     m_o_tau2_jet1                 ; 
  float     m_o_tau2_jet2                 ; 
  int       Flag_goodVertices             ; 
  int       Flag_globalTightHalo2016Filter; 
  int       Flag_HBHENoiseFilter          ; 
  int       Flag_HBHENoiseIsoFilter       ; 
  int       Flag_eeBadScFilter            ; 
  int       Flag_badChargedHadronFilter   ; 
  int       Flag_badMuonFilter            ;    
  int       Flag_ECALDeadCell             ;                                     
  float     HLTJet360_TrimMass30          ; 
  float     HLTHT700_TrimMass50           ; 
  float     HLTHT650_MJJ950DEtaJJ1p5      ; 
  float     HLTHT650_MJJ900DEtaJJ1p5      ; 
  float     HLTHT800                      ;                                        
  int       nLeptonOverlap                ; 
  int       jj_mergedVTruth               ; 
  float     b_weight                      ;
  float     b_weightGen                   ;
  float     b_weightPU                    ;
  float     b_weightBtag                  ;                 
  
  
 
  
  //naming
  std::string m_jetAK8Name;
  std::string m_genjetAK8Name;
  std::string m_jetAK8PuppiName;
  std::string m_genParticleName;
  
  
  
  // Names of the input/output trees:
  std::string m_recoTreeName;       ///< name of tree with reconstructed objects in ntuple
  std::string m_outputTreeName;    ///< name of output tree
  
  enum ValHistsType { GENERAL, ELECTRON, MUON, JETS };
  void FillValidationHists( ValHistsType, const TString& status );
  double getEventWeight( );
  void clearBranches( );
  
  std::string m_jsonName;
  bool        m_isData;
  
  
  

   // Macro adding the functions for dictionary generation
   ClassDef( VVanalysis, 0 );

}; // class VVanalysis

#endif // VVanalysis_H

