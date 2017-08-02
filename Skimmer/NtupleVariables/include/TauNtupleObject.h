// Dear emacs, this is -*- c++ -*-
// $Id: TauNtupleObject.h 37457 2010-07-05 12:04:33Z mann $

// THIS FILE HAS BEEN GENERATED AUTOMATICALLY. DO NOT EDIT DIRECTLY, CHANGES WILL BE LOST UPON NEXT CODE GENERATION.
// Code produced by Id: CodeIt.py 494 2010-07-30 13:41:32Z svn 


#ifndef SFRAME_NtupleVARIABLES_TauNtupleObject_H
#define SFRAME_NtupleVARIABLES_TauNtupleObject_H

// Local include(s):
#include "NtupleObjectNames.h"

// STL include(s):
#include <vector>
#include <string>

// ROOT include(s):
#include <TString.h>

// SFrame include(s):
#include "core/include/SError.h"
#include "core/include/SCycleBaseNTuple.h"
#include "plug-ins/include/SInputVariables.h"

namespace Ntuple {

  /**
  *  @short Class that can read the variables produced by TauNtupleObject
  *
  *         This class can be used to read the offline muon information from
  *         an ntuple produced by the SingleTopDPDMaker code.
  *
  * @author Code produced by Id: CodeIt.py 494 2010-07-30 13:41:32Z svn   
  *
  */
  
  enum TauDetails {
    TauBasic = 1,
    TauID = 2,
    TauAdvancedID = 4,
    TauAll = 7,

  };
  
  // forward declaration of NtupleObjectNames
  class NtupleObjectNames;
  class TauNtupleObject : public SInputVariables< SCycleBaseNTuple > {

    public:
    /// Constructor specifying the parent of the object
    TauNtupleObject( SCycleBaseNTuple* parent );

    /// remember if connect succeeded
    void setConnectSucceeded(const unsigned int index, const bool success);
    bool getConnectSucceeded(const unsigned int index) const {return m_connectsucceeded.at(index);}  

    /// Connect the variables to the input branches
    void ConnectVariables( const TString& treeName,
                           UInt_t detail_level = 0,
                           const TString& prefix = "Tau_",
                           const TString& ntupleType = "NtupleMakerNtuple" ) throw( SError );

    void ConnectVariables( const TString& treeName,
                           const TString& prefix = "Tau_",
                           const TString& ntupleType = "NtupleMakerNtuple" ) throw( SError );

    int getDetailLevel() const {return detailLevel;}   


    // particle vector size
    Int_t                   N;
   
    enum ConnectionIndex { 
     kdecayModeFindingNewDMs=5, 
     kdecayModeFinding=6, 
     kbyLooseCombinedIsolationDeltaBetaCorr3Hits=7, 
     kbyMediumCombinedIsolationDeltaBetaCorr3Hits=8, 
     kbyTightCombinedIsolationDeltaBetaCorr3Hits=9, 
     kbyCombinedIsolationDeltaBetaCorrRaw3Hits=10, 
     kchargedIsoPtSum=11, 
     kneutralIsoPtSum=12, 
     kpuCorrPtSum=13, 
     kbyIsolationMVA3oldDMwLTraw=14, 
     kbyVLooseIsolationMVA3oldDMwLT=15, 
     kbyLooseIsolationMVA3oldDMwLT=16, 
     kbyMediumIsolationMVA3oldDMwLT=17, 
     kbyTightIsolationMVA3oldDMwLT=18, 
     kbyVTightIsolationMVA3oldDMwLT=19, 
     kbyIsolationMVA3newDMwLTraw=20, 
     kbyVLooseIsolationMVA3newDMwLT=21, 
     kbyLooseIsolationMVA3newDMwLT=22, 
     kbyMediumIsolationMVA3newDMwLT=23, 
     kbyTightIsolationMVA3newDMwLT=24, 
     kbyVTightIsolationMVA3newDMwLT=25, 
     kagainstElectronMVA5raw=26, 
     kagainstElectronMVA5category=27, 
     kagainstElectronVLooseMVA5=28, 
     kagainstElectronLooseMVA5=29, 
     kagainstElectronMediumMVA5=30, 
     kagainstElectronTightMVA5=31, 
     kagainstElectronVTightMVA5=32, 
     kagainstMuonLoose3=33, 
     kagainstMuonTight3=34, 
     kbyPileupWeightedIsolationRaw3Hits=35, 
     kbyLoosePileupWeightedIsolation3Hits=36, 
     kbyMediumPileupWeightedIsolation3Hits=37, 
     kbyTightPileupWeightedIsolation3Hits=38, 
     kbyPhotonPtSumOutsideSignalCone=39, 
     kfootprintCorrection=40, 
     kpdgId=1, 
     kcharge=2, 
     kd0=3, 
     kTauType=4, 
 
      kEnd 
    }; 


    // vectors of properties defined in Particle.h
    std::vector< floatingnumber >  *e;
    std::vector< floatingnumber >  *pt;
    std::vector< floatingnumber >  *eta;
    std::vector< floatingnumber >  *phi;
    std::vector< floatingnumber >  *m;
    

    
    // vectors of object specific variables
    std::vector< floatingnumber >  *decayModeFindingNewDMs;
    std::vector< floatingnumber >  *decayModeFinding;
    std::vector< floatingnumber >  *byLooseCombinedIsolationDeltaBetaCorr3Hits;
    std::vector< floatingnumber >  *byMediumCombinedIsolationDeltaBetaCorr3Hits;
    std::vector< floatingnumber >  *byTightCombinedIsolationDeltaBetaCorr3Hits;
    std::vector< floatingnumber >  *byCombinedIsolationDeltaBetaCorrRaw3Hits;
    std::vector< floatingnumber >  *chargedIsoPtSum;
    std::vector< floatingnumber >  *neutralIsoPtSum;
    std::vector< floatingnumber >  *puCorrPtSum;
    std::vector< floatingnumber >  *byIsolationMVA3oldDMwLTraw;
    std::vector< floatingnumber >  *byVLooseIsolationMVA3oldDMwLT;
    std::vector< floatingnumber >  *byLooseIsolationMVA3oldDMwLT;
    std::vector< floatingnumber >  *byMediumIsolationMVA3oldDMwLT;
    std::vector< floatingnumber >  *byTightIsolationMVA3oldDMwLT;
    std::vector< floatingnumber >  *byVTightIsolationMVA3oldDMwLT;
    std::vector< floatingnumber >  *byIsolationMVA3newDMwLTraw;
    std::vector< floatingnumber >  *byVLooseIsolationMVA3newDMwLT;
    std::vector< floatingnumber >  *byLooseIsolationMVA3newDMwLT;
    std::vector< floatingnumber >  *byMediumIsolationMVA3newDMwLT;
    std::vector< floatingnumber >  *byTightIsolationMVA3newDMwLT;
    std::vector< floatingnumber >  *byVTightIsolationMVA3newDMwLT;
    std::vector< floatingnumber >  *againstElectronMVA5raw;
    std::vector< floatingnumber >  *againstElectronMVA5category;
    std::vector< floatingnumber >  *againstElectronVLooseMVA5;
    std::vector< floatingnumber >  *againstElectronLooseMVA5;
    std::vector< floatingnumber >  *againstElectronMediumMVA5;
    std::vector< floatingnumber >  *againstElectronTightMVA5;
    std::vector< floatingnumber >  *againstElectronVTightMVA5;
    std::vector< floatingnumber >  *againstMuonLoose3;
    std::vector< floatingnumber >  *againstMuonTight3;
    std::vector< floatingnumber >  *byPileupWeightedIsolationRaw3Hits;
    std::vector< floatingnumber >  *byLoosePileupWeightedIsolation3Hits;
    std::vector< floatingnumber >  *byMediumPileupWeightedIsolation3Hits;
    std::vector< floatingnumber >  *byTightPileupWeightedIsolation3Hits;
    std::vector< floatingnumber >  *byPhotonPtSumOutsideSignalCone;
    std::vector< floatingnumber >  *footprintCorrection;
    std::vector< int >  *pdgId;
    std::vector< floatingnumber >  *charge;
    std::vector< floatingnumber >  *d0;
    std::vector< int >  *TauType;


    std::vector<int> m_connectsucceeded;

    // save actual detail level and group
    Int_t detailLevel;
    std::string detailGroup;

  }; // class TauNtupleObject

} // namespace Ntuple

#endif // SFRAME_NtupleVARIABLES_TauNtupleObject_H
