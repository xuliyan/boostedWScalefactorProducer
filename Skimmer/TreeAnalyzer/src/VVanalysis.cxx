
// Local include(s):
#include "../include/VVanalysis.h"
#include "../include/VVAnalysisTools.h"

// External include(s):
// #include "../GoodRunsLists/include/TGoodRunsListReader.h"


#include <TMath.h>
#include "TGraph.h"

ClassImp( VVanalysis );

VVanalysis::VVanalysis()
  : SCycleBase()
   , m_jetAK4               ( this )
   , m_jetAK8               ( this )
   , m_jetAK8Puppi          ( this )
   , m_eventInfo            ( this )
   , m_electron             ( this )
   , m_muon                 ( this )
   , m_missingEt            ( this )
   , m_genParticle          ( this )
   , m_pileupReweightingTool( this )
   , m_bTaggingScaleTool    ( this )
   , m_allEvents( "allEvents", this )
   , m_passedLoose( "passedLoose", this )
   , m_passedPuppi( "passedPuppi", this )
   , m_passedMjj( "passedMjj", this )
   , m_passedDEta( "passedDeta", this )
   , m_passedEvents( "passedEvents", this )
   , m_test( "test", this ) 
   , m_genjetAK8( this )
 {
   m_logger << INFO << "Initializing" << SLogger::endmsg;
   SetLogName( GetName() );
   
   DeclareProperty( "JetAK8Name"      ,         m_jetAK8Name              = "jetAK8"        );
   DeclareProperty( "GenJetAK8Name"   ,         m_genjetAK8Name           = "genJetAK8"     );
   DeclareProperty( "JetAK8PuppiName" ,         m_jetAK8PuppiName         = "jetAK8_puppi"  );
   DeclareProperty( "GenParticleName" ,         m_genParticleName         = "genParticle"   );
   
   DeclareProperty( "RecoTreeName"    ,         m_recoTreeName            = "physics" );
   DeclareProperty( "PUPPIJEC"        ,         m_PUPPIJEC                = "weights/puppiCorr.root" );
   // DeclareProperty( "JSONName"        ,         m_jsonName                = std::string (std::getenv("SFRAME_DIR")) + "/../GoodRunsLists/JSON/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver_v2.txt" );
   DeclareProperty( "IsData"          ,         m_isData                  = false );
}

VVanalysis::~VVanalysis() {
  
  m_logger << INFO << "Done!" << SLogger::endmsg;

}

void VVanalysis::BeginCycle() throw( SError ) { //Block used to perform an initial configuration of the cycle (ead some local file likegood data ranges)
  
  // // Load the GRL:
  // if (m_isData) {
  //   m_logger << INFO << "Loading GoodRunsList from file: " << m_jsonName << SLogger::endmsg;
  //   Root::TGoodRunsListReader reader( TString( m_jsonName.c_str() ) );
  //   if( ! reader.Interpret() ) {
  //     m_logger << FATAL << "Couldn't read in the GRL!" << SLogger::endmsg;
  //     throw SError( ( "Couldn't read in file: " + m_jsonName ).c_str(), SError::SkipCycle );
  //   }
  //   m_grl = reader.GetMergedGoodRunsList();
  //   m_grl.Summary();
  //   m_grl.SetName( "MyGoodRunsList" );
  //
  //   // Add it as a configuration object, so that the worker nodes will receive it:
  //   AddConfigObject( &m_grl );
  // }

   return;

}

void VVanalysis::EndCycle() throw( SError ) { //Any finalisation steps should be done here. (Closure of some helper files opened by the user code for instance.)

   return;

}

void VVanalysis::BeginInputFile( const SInputData& ) throw( SError ) { //For each new input file the user has to connect his input variables. (More on this later.) This has to be performed in this function.
  
  m_logger << INFO << "Connecting input variables" << SLogger::endmsg;
  if (m_isData) {
    m_jetAK4.ConnectVariables(       m_recoTreeName.c_str(), Ntuple::JetBasic|Ntuple::JetAnalysis|Ntuple::JetTruth, (m_jetAK4Name + "_").c_str() );
    m_electron.ConnectVariables(     m_recoTreeName.c_str(), Ntuple::ElectronBasic|Ntuple::ElectronID, (m_electronName + "_").c_str() );
    m_muon.ConnectVariables(         m_recoTreeName.c_str(), Ntuple::MuonBasic|Ntuple::MuonID|Ntuple::MuonIsolation, (m_muonName + "_").c_str() );
    m_missingEt.ConnectVariables(    m_recoTreeName.c_str(), Ntuple::MissingEtBasic, (m_missingEtName + "_").c_str() );
    m_jetAK8Puppi .ConnectVariables(  m_recoTreeName.c_str(), Ntuple::JetSubstructure, (m_jetAK8PuppiName + "_").c_str() );
    m_jetAK8      .ConnectVariables(  m_recoTreeName.c_str(), Ntuple::JetBasic|Ntuple::JetAnalysis|Ntuple::JetSubstructure|Ntuple::JetSoftdropSubjets, (m_jetAK8Name + "_").c_str() );
    m_eventInfo   .ConnectVariables(  m_recoTreeName.c_str(), Ntuple::EventInfoBasic|Ntuple::EventInfoTrigger|Ntuple::EventInfoMETFilters, "" );
  }
  else {
    m_jetAK4.ConnectVariables(       m_recoTreeName.c_str(), Ntuple::JetBasic|Ntuple::JetAnalysis|Ntuple::JetTruth, (m_jetAK4Name + "_").c_str() );
    m_electron.ConnectVariables(     m_recoTreeName.c_str(), Ntuple::ElectronBasic|Ntuple::ElectronID, (m_electronName + "_").c_str() );
    m_muon.ConnectVariables(         m_recoTreeName.c_str(), Ntuple::MuonBasic|Ntuple::MuonID|Ntuple::MuonIsolation, (m_muonName + "_").c_str() );
    m_missingEt.ConnectVariables(    m_recoTreeName.c_str(), Ntuple::MissingEtBasic, (m_missingEtName + "_").c_str() );
    m_jetAK8Puppi .ConnectVariables(  m_recoTreeName.c_str(), Ntuple::JetSubstructure, (m_jetAK8PuppiName + "_").c_str() );
    m_jetAK8      .ConnectVariables(  m_recoTreeName.c_str(), Ntuple::JetBasic|Ntuple::JetAnalysis|Ntuple::JetSubstructure|Ntuple::JetTruth|Ntuple::JetSoftdropSubjets|Ntuple::JetSoftdropSubjetsTruth, (m_jetAK8Name + "_").c_str() );
    m_eventInfo   .ConnectVariables(  m_recoTreeName.c_str(), Ntuple::EventInfoBasic|Ntuple::EventInfoTrigger|Ntuple::EventInfoMETFilters|Ntuple::EventInfoTruth, "" );
    m_genParticle .ConnectVariables(  m_recoTreeName.c_str(), Ntuple::GenParticleBasic, (m_genParticleName + "_").c_str() );
    m_genjetAK8   .ConnectVariables(  m_recoTreeName.c_str(), Ntuple::GenJet, (m_genjetAK8Name + "_").c_str() );
  }
  // Unused collections
  // m_jetAK4.ConnectVariables(       m_recoTreeName.c_str(), Ntuple::JetBasic|Ntuple::JetAnalysis|Ntuple::JetTruth, (m_jetAK4Name + "_").c_str() );
  // m_electron.ConnectVariables(     m_recoTreeName.c_str(), Ntuple::ElectronBasic|Ntuple::ElectronID, (m_electronName + "_").c_str() );
  // m_muon.ConnectVariables(         m_recoTreeName.c_str(), Ntuple::MuonBasic|Ntuple::MuonID|Ntuple::MuonIsolation, (m_muonName + "_").c_str() );
  // m_missingEt.ConnectVariables(    m_recoTreeName.c_str(), Ntuple::MissingEtBasic, (m_missingEtName + "_").c_str() );
  
  m_logger << INFO << "Connecting input variables completed" << SLogger::endmsg;
  
  m_logger << INFO << "Initializing inputs for puppi softdop mass corrections" << SLogger::endmsg; 
  TFile* file = TFile::Open( m_PUPPIJEC.c_str(),"READ");
  if(file->IsZombie()) throw SError(SError::StopExecution);
  m_puppisd_corr.push_back( dynamic_cast<TF1*>(file->Get("puppiJECcorr_gen")));
  m_puppisd_corr.push_back( dynamic_cast<TF1*>(file->Get("puppiJECcorr_reco_0eta1v3")));
  m_puppisd_corr.push_back( dynamic_cast<TF1*>(file->Get("puppiJECcorr_reco_1v3eta2v5")));
  
  return;

}


void VVanalysis::BeginInputData( const SInputData& ) throw( SError ) { //called once before processing each of the input data types. If you need to initialise output objects (histograms, etc.) before the event-by-event execution, you should do that here. Also the declaration of the output variables has to be done here.
  
  // VARIABLES TO BE IMPLEMENTED FOR V-TAG CODE!!!
  //DeclareVariable( puweight_     , "puweight"   );
  //DeclareVariable( lumiweight_   , "lumiweight" );
  //DeclareVariable( genweight_    , "genweight"  );
  //DeclareVariable( btagweight_ 	 , "btagweight" );
  //DeclareVariable( ptweight_ 	   , "ptweight"   );
  //DeclareVariable( weight	   	   , "weight"	    );
  //DeclareVariable( jet_csv            , "Whadr_csv"    );
  //DeclareVariable( jet_pt             , "Whadr_pt"     );
  //DeclareVariable( jet_eta            , "Whadr_eta"    );
  //DeclareVariable( jet_phi            , "Whadr_phi"    );
  //DeclareVariable( realW              , "Whadr_isW"    );
  //DeclareVariable( jet_puppi_softdrop , "Whadr_puppi_softdrop"   );
  //DeclareVariable( jet_puppi_tau1     , "Whadr_puppi_tau1"       );
  //DeclareVariable( jet_puppi_tau2     , "Whadr_puppi_tau2"       );
  //DeclareVariable( jet_puppi_tau3     , "Whadr_puppi_tau3"       );
 
  //DeclareVariable( l_pt               , "lept_pt"      );
  //DeclareVariable( l_eta              , "lept_eta"     );
  //DeclareVariable( l_phi              , "lept_phi"     );
  //DeclareVariable( l_iso              , "lept_iso"     );
  //DeclareVariable( Wlept_pt           , "Wlept_pt"     ); 
  //DeclareVariable( WMass              , "Wlept_mass"   );
  //DeclareVariable( pfMET              , "pfMET"        ); 
  //DeclareVariable( pfMETPhi           , "pfMETPhi"     );
  //DeclareVariable( nak4jets           , "nak4jets"     );
  
  
  
  
  
  
  DeclareVariable( m_o_mjj                        , "mjj");
  DeclareVariable( m_o_genmjj                     , "genmjj");
  DeclareVariable( m_o_mpuppisoftdrop_jet1        , "jet1_softDrop_mass");
  DeclareVariable( m_o_mpuppisoftdrop_jet2        , "jet2_softDrop_mass");
  DeclareVariable( m_o_mgensoftdrop_jet1          , "jet1_gen_softDrop_mass");
  DeclareVariable( m_o_mgensoftdrop_jet2          , "jet2_gen_softDrop_mass");
  DeclareVariable( m_o_tau1_jet1                  , "jet1_tau1");  
  DeclareVariable( m_o_tau1_jet2                  , "jet2_tau1");
  DeclareVariable( m_o_tau2_jet1                  , "jet1_tau2");  
  DeclareVariable( m_o_tau2_jet2                  , "jet2_tau2");
  DeclareVariable( Flag_goodVertices              , "Flag_goodVertices");
  DeclareVariable( Flag_globalTightHalo2016Filter , "Flag_globalTightHalo2016Filter");
  DeclareVariable( Flag_HBHENoiseFilter           , "Flag_HBHENoiseFilter");
  DeclareVariable( Flag_HBHENoiseIsoFilter        , "Flag_HBHENoiseIsoFilter");
  DeclareVariable( Flag_eeBadScFilter             , "Flag_eeBadScFilter");
  DeclareVariable( Flag_badChargedHadronFilter    , "Flag_badChargedHadronFilter");
  DeclareVariable( Flag_badMuonFilter             , "Flag_badMuonFilter");
  DeclareVariable( Flag_ECALDeadCell              , "Flag_ECALDeadCell");
  DeclareVariable( HLTJet360_TrimMass30           , "HLTJet360_TrimMass30");
  DeclareVariable( HLTHT700_TrimMass50            , "HLTHT700_TrimMass50");
  DeclareVariable( HLTHT650_MJJ950DEtaJJ1p5       , "HLTHT650_MJJ950DEtaJJ1p5");
  DeclareVariable( HLTHT650_MJJ900DEtaJJ1p5       , "HLTHT650_MJJ900DEtaJJ1p5");
  DeclareVariable( HLTHT800                       , "HLTHT800");
  DeclareVariable( nLeptonOverlap                 , "nLeptonOverlap");
  DeclareVariable( jj_mergedVTruth_jet1           , "jj_mergedVTruth_jet1");
  DeclareVariable( jj_mergedVTruth_jet2           , "jj_mergedVTruth_jet2");

  // // Declare the output histograms:
  // Book( TH1F( "Mjj_hist"     , "Mjj", 100, 0.0, 3000. ) );
  
 
  m_test->resize( 6, 0 ); // Reserve two entries in vector for counting entries:
              
   return;
}

void VVanalysis::EndInputData( const SInputData& ) throw( SError ) {
  
  // static const Int_t n = 5;
 //  Float_t x_array[ n ] = { 0.0, 1.0, 2.0, 3.0, 4.0 };
 //  Float_t y_array[ n ] = { 0.0, 2.0, 4.0, 6.0, 8.0 };
 //  TGraph mygraph( n, x_array, y_array );
 //  mygraph.SetName( "MyGraph" );
 //  WriteObj( mygraph, "graph_dir" );

   return;

}



void VVanalysis::BeginMasterInputData( const SInputData& ) throw( SError ){

   return;

}

void VVanalysis::EndMasterInputData( const SInputData& ) throw( SError ){ //this is a good place to print some summaries, do some final calculations on the created histograms (for instance fitting them), etc
  
  m_logger << INFO << "Number of all processed events      :  "<< *m_allEvents   << "   " << ( m_test->size() > 0 ? ( *m_test )[ 0 ] : 0 )<< SLogger::endmsg;
  m_logger << INFO << "Number of events passing pt, eta    :  "<< *m_passedLoose << "   " << ( m_test->size() > 1 ? ( *m_test )[ 1 ] : 0 )<< SLogger::endmsg;
  m_logger << INFO << "Found two PUPPI jets matched to AK8 :  "<< *m_passedPuppi << "   " << ( m_test->size() > 2 ? ( *m_test )[ 2 ] : 0 )<< SLogger::endmsg;
  m_logger << INFO << "Number of events passing Mjj        :  "<< *m_passedMjj   << "   " << ( m_test->size() > 3 ? ( *m_test )[ 3 ] : 0 )<< SLogger::endmsg;
  m_logger << INFO << "Number of events passing dETa       :  "<< *m_passedDEta  << "   " << ( m_test->size() > 4 ? ( *m_test )[ 4 ] : 0 )<< SLogger::endmsg;
  m_logger << INFO << "Number of events passing selection  :  "<< *m_passedEvents<< "   " << ( m_test->size() > 5 ? ( *m_test )[ 5 ] : 0 )<< SLogger::endmsg;
  
  
  
   return;

}
  
  


void VVanalysis::ExecuteEvent( const SInputData&, Double_t weight) throw( SError ) { //This is the main analysis function that is called for each event. It receives the weight of the event, as it is calculated by the framework from the luminosities and generator cuts defined in the XML configuration.

  clearBranches();
  
  ++m_allEvents;
  ( *m_test )[ 0 ]++;
  
  
  // W-top scalefactor pseudocode/old code to be rewritten
  // setWeight();
  // if( !applyJSON() ) throw SError( SError::SkipEvent );
  // if( !passedTrigger() ) throw SError( SError::SkipEvent );
  // if( !passedFilters( infile ) ) throw SError( SError::SkipEvent );
  // if( passedTTbarSelections( infile) ){ <--- SEE METHOD BELOW!!!
    
   //GEN MATCHING
    // int isW_def1 = 0;
    // int isW_def2 = 0;
    //
    // if( runOnMC_ ){
    //   TLorentzVector genP;
    //   std::vector<TLorentzVector> daus;
    //   TLorentzVector hadW;
    //
    //   for( unsigned int p = 0; p < (data_.genParticle_pdgId)->size(); ++p ){
    //     bool isHad = false;
    //     if( (*data_.genParticle_pt).at(p) < 0.1) continue;
    //     if( fabs((*data_.genParticle_pdgId).at(p)) > 9) continue;
    //     for( unsigned int d = 0; d < (*data_.genParticle_mother)[p].size(); ++d ) {
    //       if( fabs((*data_.genParticle_mother)[p][d]) != 24 ) continue;
    //       isHad = true;
    //     }
    //     if (!isHad) continue;
    //     genP.SetPtEtaPhiE( (*data_.genParticle_pt).at(p), (*data_.genParticle_eta).at(p), (*data_.genParticle_phi).at(p), (*data_.genParticle_e).at(p) );
    //     daus.push_back(genP);
    //   }
    //   if( daus.size() > 1){
    //
    //     for( unsigned int p = 0; p < (data_.genParticle_pdgId)->size(); ++p ){
    //       bool isHad = false;
    //       if( (*data_.genParticle_pt).at(p) < 0.1) continue;
    //       if( fabs((*data_.genParticle_pdgId).at(p)) != 24) continue;
    //       // if( not infile.Contains( "herwig" ) && (*data_.genParticle_status).at(p) != 4 && (*data_.genParticle_status).at(p) != 22 && (*data_.genParticle_status).at(p) != 23) continue;
    //       // if( infile.Contains( "herwig" ) && (*data_.genParticle_status).at(p) != 11 ) continue;
    //       if ((*data_.genParticle_dau)[p].size() < 2) continue;
    //       for( unsigned int d = 0; d < (*data_.genParticle_dau)[p].size(); ++d ) {
    //         if( fabs((*data_.genParticle_dau)[p][d]) > 9) continue;
    //         isHad = true;
    //       }
    //       if (isHad){
    //         hadW.SetPtEtaPhiE( (*data_.genParticle_pt).at(p), (*data_.genParticle_eta).at(p), (*data_.genParticle_phi).at(p), (*data_.genParticle_e).at(p) );
    //       }
    //     }
    //
    //     int isMatch = 0;
    //     for( unsigned int p = 0; p < daus.size(); ++p ){
    //       float dR = daus.at(p).DeltaR(Vcand.at(0).p4);
    //       if (dR < 0.8 ) isMatch +=1;
    //     }
    //     if (isMatch > 1) isW_def1 = 1;
    //     if( (Vcand.at(0).p4).DeltaR(hadW) < 0.4  && daus.at(0).DeltaR(daus.at(1)) < 0.8 ) isW_def2 = 1;
    //   }
    // }
    //
    // realW_def1 = isW_def1;
    // realW_def2 = isW_def2;
    // weight = weight_;
    // jet_csv = Vcand.at(0).csv;
    // jet_pt = Vcand.at(0).p4.Pt();
    // jet_eta = Vcand.at(0).p4.Eta();
    // jet_phi = Vcand.at(0).p4.Phi();
    // l_pt  = leptonCand_[0].p4.Pt();
    // l_eta = leptonCand_[0].p4.Eta();
    // l_phi = leptonCand_[0].p4.Phi();
    // l_iso = leptonCand_[0].iso;
    // Wlept_pt = WCand_[0].p4.Pt();
    // pfMET = METCand_[0].p4.Pt();
    // pfMETPhi = METCand_[0].p4.Phi();
    // nak4jets = AK4jetCand_.size();
    // jet_puppi_softdrop = Vcand.at(0).puppi_softdropMass;
    // jet_puppi_tau1 = Vcand.at(0).puppi_tau1;
    // jet_puppi_tau2 = Vcand.at(0).puppi_tau2;
    // jet_puppi_tau3 = Vcand.at(0).puppi_tau3;

    
    
  
   
   
   if( m_jetAK8.N < 2 ) throw SError( SError::SkipEvent );
   
  //-------------Select two fat jets-------------//
  std::vector<UZH::Jet> goodFatJets;
  std::vector<UZH::Jet> goodGenJets;
  std::vector<UZH::GenParticle> GenQuarks = FindGeneratedQuarks(m_genParticle, m_isData);

  for ( int i = 0; i < (m_jetAK8.N); ++i ) {
   
    UZH::Jet myjet( &m_jetAK8, i );
    if (! (myjet.pt() > 200       )) continue;
    if (! (fabs(myjet.eta()) < 2.5)) continue;
    if (! (myjet.IDTight()        )) continue;
  
    //Match to puppi jet
    float dRmin = 99.;
    for ( int ii = 0; ii < abs((m_jetAK8Puppi.pt->size())); ++ii ) {
      UZH::Jet mypuppijet( &m_jetAK8Puppi, ii );
      float dR = myjet.DeltaR(mypuppijet);
      if ( dR > dRmin ) continue;
      dRmin = dR;
      myjet.puppi_softdropmass= ApplyPuppiSoftdropMassCorrections(mypuppijet,m_puppisd_corr,m_isData);//mypuppijet.softdrop_massCorr();
      myjet.puppi_tau1        = mypuppijet.tau1();
      myjet.puppi_tau2        = mypuppijet.tau2();
    }
    
    //Match to gen jet
    if(!m_isData){
      dRmin = 99.;
      int jetIdx = 99;
      for ( int j = 0; j < (m_genjetAK8.N); ++j ) {
        UZH::Jet genJet( &m_genjetAK8, j );
        if ( genJet.pt() < 1. ) continue;
        float dR = genJet.DeltaR(myjet);
        if ( dR > dRmin ) continue;
        dRmin = dR;
        jetIdx = j;
      }
      UZH::Jet selectedJet( &m_genjetAK8, jetIdx );
      goodGenJets.push_back(selectedJet);
    }
    goodFatJets.push_back(myjet);
  }
  //-------------Select two fat jets-------------//
  if( goodFatJets.size() < 2 ) throw SError( SError::SkipEvent );
  ++m_passedLoose;
  ( *m_test )[ 1 ]++;

     
   //-------------Make sure we have a PUPPI match-------------// 
  if( goodFatJets[0].puppi_softdropmass == -99 || goodFatJets[1].puppi_softdropmass == -99
    || goodFatJets[0].puppi_tau1 == -99         || goodFatJets[1].puppi_tau1 == -99
    || goodFatJets[0].puppi_tau2 == -99         || goodFatJets[1].puppi_tau2 == -99
     ) throw SError( SError::SkipEvent );
  ++m_passedPuppi;
  ( *m_test )[ 2 ]++;

  // Loose Mjj cut to slim down samples
  TLorentzVector dijet = goodFatJets[0].tlv() + goodFatJets[1].tlv();
  if( dijet.M() < 900. ) throw SError( SError::SkipEvent );
  ++m_passedMjj;
  ( *m_test )[ 3 ]++;

  // dEta cut
  if( abs( (goodFatJets[0].tlv()).Eta() - (goodFatJets[1].tlv()).Eta() )  > 1.3 ) throw SError( SError::SkipEvent );
  ++m_passedDEta;
  ( *m_test )[ 4 ]++;
  
  ++m_passedEvents;
  ( *m_test )[ 5 ]++;

  // Fill tree 
  m_o_mjj                  = dijet.M();
  m_o_mpuppisoftdrop_jet1  = goodFatJets[0].puppi_softdropmass;
  m_o_mpuppisoftdrop_jet2  = goodFatJets[1].puppi_softdropmass;
  m_o_tau1_jet1            = goodFatJets[0].puppi_tau1;
  m_o_tau1_jet2            = goodFatJets[1].puppi_tau2;
  m_o_tau2_jet1            = goodFatJets[0].puppi_tau1;
  m_o_tau2_jet2            = goodFatJets[1].puppi_tau2;


  bool isfired = false;
  for( std::map<std::string,bool>::iterator it = (m_eventInfo.trigDecision)->begin(); it != (m_eventInfo.trigDecision)->end(); ++it){
    if( (it->first).find("AK8PFJet360_TrimMass30")          != std::string::npos) HLTJet360_TrimMass30      = it->second;
    if( (it->first).find("AK8PFHT700_TrimR0p1PT0p03Mass50") != std::string::npos) HLTHT700_TrimMass50       = it->second;
    if( (it->first).find("PFHT650_WideJetMJJ950DEtaJJ1p5")  != std::string::npos) HLTHT650_MJJ950DEtaJJ1p5  = it->second;
    if( (it->first).find("PFHT650_WideJetMJJ900DEtaJJ1p5")  != std::string::npos) HLTHT650_MJJ900DEtaJJ1p5  = it->second;
    if( (it->first).find("PFHT800_v")                       != std::string::npos) HLTHT800                  = it->second;
  }

  Flag_goodVertices               = m_eventInfo.PV_filter;
  Flag_globalTightHalo2016Filter  = m_eventInfo.passFilter_CSCHalo;
  Flag_HBHENoiseFilter            = m_eventInfo.passFilter_HBHELoose;
  Flag_HBHENoiseIsoFilter         = m_eventInfo.passFilter_HBHEIso;
  Flag_eeBadScFilter              = m_eventInfo.passFilter_EEBadSc;
  Flag_badChargedHadronFilter     = m_eventInfo.passFilter_chargedHadronTrackResolution;
  Flag_badMuonFilter              = m_eventInfo.passFilter_muonBadTrack;
  Flag_ECALDeadCell               = m_eventInfo.passFilter_ECALDeadCell;

  if(!m_isData){
    m_o_genmjj                      = (goodGenJets[0].tlv() + goodGenJets[1].tlv()).M();
    m_o_mgensoftdrop_jet1           = goodGenJets[0].softdropmass();
    m_o_mgensoftdrop_jet2           = goodGenJets[1].softdropmass();
  }
  
  // nLeptonOverlap                  =           ;
  jj_mergedVTruth_jet1                =  isMergedVJet(goodFatJets[0].tlv(),GenQuarks) ;
  jj_mergedVTruth_jet2                =  isMergedVJet(goodFatJets[1].tlv(),GenQuarks) ;
  //

   
  
  // }
  // TString StatusString = "Before_Cuts_";
  // FillValidationHists( GENERAL, StatusString );
  // FillValidationHists( ELECTRON, StatusString );
  // FillValidationHists( MUON, StatusString );
  // FillValidationHists( JETS, StatusString );

  return;

}

 

void VVanalysis::FillValidationHists( ValHistsType ht, const TString& status ) {  

   if( ht == GENERAL ){
      TString suffix = status + "General_";
   } 
   else if( ht == ELECTRON ){
      TString suffix = status + "Electron_";
   } 
   else if( ht == MUON ){
      TString suffix = status + "Muon_";
   } 
   else if( ht == JETS ){
      TString suffix = status + "Jets_";
      Book( TH1F( suffix + "mjj", "mjj", 300, 0.0, 3000.0 ) )->Fill( m_o_mjj );
   } 
   else {
      SError error( SError::StopExecution );
      error << GetName() << ": Validation histogram type"  << ht << " not supported";
      throw error;
   }

   return;
}

double VVanalysis::getEventWeight() {
  
  double weight = 1.;
  for( unsigned int v = 0; v < (m_eventInfo.actualIntPerXing)->size(); ++v ){
    if ( (*m_eventInfo.bunchCrossing)[v] == 0 ) {
      b_weightPU = m_pileupReweightingTool.getPileUpweight( (*m_eventInfo.actualIntPerXing)[v] );
      m_logger << VERBOSE << "Weight: " << b_weightPU << " for true: " << (*m_eventInfo.actualIntPerXing)[v] << SLogger::endmsg;
      break;
    }
  }
  b_weightGen = (m_eventInfo.genEventWeight < 0) ? -1 : 1; 
  weight *= b_weightPU*b_weightGen;
  
  return weight;
  
}


void VVanalysis::clearBranches() {
  
  b_weight = 1.;
  b_weightGen = 1.;
  b_weightPU = 1.;
  b_weightBtag = 1.;
  
  m_o_mjj                       = -99;
  m_o_genmjj                    = -99;
  m_o_mpuppisoftdrop_jet1       = -99;
  m_o_mpuppisoftdrop_jet2       = -99;
  m_o_mgensoftdrop_jet1         = -99;
  m_o_mgensoftdrop_jet2         = -99;
  m_o_tau1_jet1                 = -99;
  m_o_tau1_jet2                 = -99;
  m_o_tau2_jet1                 = -99;
  m_o_tau2_jet2                 = -99;
  Flag_goodVertices             = -99;
  Flag_globalTightHalo2016Filter= -99;
  Flag_HBHENoiseFilter          = -99;
  Flag_HBHENoiseIsoFilter       = -99;
  Flag_eeBadScFilter            = -99;
  Flag_badChargedHadronFilter   = -99;
  Flag_badMuonFilter            = -99;
  Flag_ECALDeadCell             = -99;
  HLTJet360_TrimMass30          = -99;
  HLTHT700_TrimMass50           = -99;
  HLTHT650_MJJ950DEtaJJ1p5      = -99;
  HLTHT650_MJJ900DEtaJJ1p5      = -99;
  HLTHT800                      = -99;
  nLeptonOverlap                = -99;
  jj_mergedVTruth_jet1          = -99;
  jj_mergedVTruth_jet2          = -99;

}

// bool ExoDiBosonAnalysis::passedTTbarSelections( TString infile ){
//
//   bool foundLept    = false;
//   bool passVeto     = false;
//   bool foundMET     = false;
//   bool foundW       = false;
//   bool foundAK8     = false;
//   bool passAK8LepDR = false;
//   bool passAK8PrunedMass = false;
//   bool found1b      = false;
//
//   foundAMu  = 0;
//   foundAEle = 0;
//
//   // Find lepton (muon/electron depending on channel)
//   if( findLeptonCandidate() ) foundLept = true;
//   if( !foundLept ) return false;
//   nPassedFoundLept_++;
//
//   Hist( "leptonPT"    )->Fill( leptonCand_.at(0).p4.Pt()  , weight_ );
//   Hist( "leptonPhi"   )->Fill( leptonCand_.at(0).p4.Phi(), weight_ );
//   Hist( "leptonEta"   )->Fill( leptonCand_.at(0).p4.Eta(), weight_ );
//
//
//   // Veto additional loose leptons (both muons/electrons )
//   if ( (foundAMu+foundAEle) > 0) return false;
//   else
//     passVeto = true;
//
//   nPassedVetoLep_++;
//
//   // Find MET
//   if( findMETCandidate() ) {
//     Hist( "MET"         )->Fill(  METCand_.at(0).sumEt, weight_ );
//      if( METCand_[0].p4.Pt() > METCut_ ) foundMET = true;
//    }
//   if( !foundMET ) return false;
//   nPassedFoundMET_++;
//
//   // reconstruct W
//   if( foundLept && foundMET && findWCandidate() ){
//     Hist( "Wcand_pt")->Fill(  WCand_[0].p4.Pt() );
//      if( WCand_[0].p4.Pt() > 200.0 ) foundW = true;
//   }
//
//   if( !foundW ) return false;
//   nPassedFoundW_++;
//   // Find AK8 W-jet candidate
//   if(findJetCandidates( infile )){
//     foundAK8 = true;
//     // if( leptonCand_[0].p4.DeltaR(Vcand.at(0).p4) > 1.0 ) passAK8LepDR = true;
//     if( Vcand.at(0).p4.Pt() >= JetPtCutTight_  ) passAK8LepDR = true;
//     float groomedMass = Vcand.at(0).prunedMass;
//     if(usePuppiSD_) groomedMass = Vcand.at(0).puppi_softdropMass;
//     if( groomedMass >= mWLow_ && groomedMass <= mZHigh_ ) passAK8PrunedMass = true;
//   }
//
//   if( !foundAK8 ) return false;
//   nPassedFoundJet_++;
//
//   if( !passAK8LepDR ) return false;
//   nPassedLepJetDR_++;
//
//   if( !passAK8PrunedMass ) return false;
//   nPassedJetPrunedMass_++;
//
//
//   // Find b-tagged AK4 jet
//   findAK4Jets();
//   bool najetCut = false;
//   int najet = 0;
//   for( unsigned int j = 0; j < abs(AK4jetCand_.size()); j++ ){
//      if( AK4jetCand_[j].csv > 0.800 ) najet++;
//   }
//
//   if( najet  > 0 ) found1b = true;
//   if( !found1b ) return false;
//   nPassed1bTag_++;
//
//   if(foundLept && passVeto && foundMET && foundW && foundAK8 && passAK8LepDR && passAK8PrunedMass && found1b) return true;
//
//   else
//     return false;
// }
//
// bool ExoDiBosonAnalysis::findLeptonCandidate( void ){
//
//   if( Channel_.find("el") != std::string::npos ) return findElectronCandidate();
//   else if( Channel_.find("mu") != std::string::npos ) return findMuonCandidate();
//   else return false;
//
// }
//
// //==============================================================================================
// bool ExoDiBosonAnalysis::findMuonCandidate( void ){
//
//   foundAEle = 0;
//   foundAMu = 0;
//
//   TLorentzVector TLV;
//   bool foundLept = false;
//   double ptMu = -999;
//   int muIndex = 999;
//   double scale = 0.;
//   double newpt = 0.;
//   double scale_ = 0.;
//   for( int l = 0; l < data_.mu_N; ++l ){
//     // if( !isLowMass && (*data_.mu_isHighPtMuon)[l] != 1 ) continue;
//     // if(  isLowMass && (*data_.mu_isTightMuon)[l] != 1 ) continue;
//     if( (*data_.mu_isTightMuon)[l] != 1 ) continue;
//     //if( (*data_.mu_isPFMuon)[l] != 1 ) continue;
//     scale = getMuonScale((*data_.mu_pt)[l]);
//     newpt = (*data_.mu_pt)[l]+scale;
//     if( (*data_.mu_trackIso)[l]/newpt >= 0.1 ) continue;
//     if( newpt <= leptPtCut_ ) continue;
//     if( fabs((*data_.mu_eta)[l]) >= leptEtaCut_ ) continue;
//     if( fabs((*data_.mu_eta)[l]) > 1.2 && newpt > 500 ) continue;
//     foundLept = true;
//     if( newpt > ptMu ){
//       ptMu = newpt;
//       TLV.SetPtEtaPhiE( newpt, (*data_.mu_eta)[l], (*data_.mu_phi)[l], (*data_.mu_e)[l]+scale );
//       muIndex = l;
//       scale_ = scale;
//     }
//   }
//
//    if( foundLept ){
//      LeptonCandidate muCand(TLV);
//      muCand.iso = (*data_.mu_trackIso)[muIndex]/ptMu;
//      // muCand.scale = scale_;
//  //     muCand.isGlobalMuon = (*data_.mu_isGlobalMuon)[muIndex];
//  //     muCand.globalHits = (*data_.mu_globalHits)[muIndex];
//  //     muCand.matchedStations = (*data_.mu_matchedStations)[muIndex];
//  //     muCand.bestTrack_ptErrRel = (*data_.mu_bestTrack_ptErr)[muIndex]/(*data_.mu_bestTrack_pt)[muIndex];
//  //     muCand.d0 = (*data_.mu_d0)[muIndex];
//  //     //muCand.dz = (*data_.mu_dz)[muIndex];
//  //     muCand.pixelHits = (*data_.mu_pixelHits)[muIndex];
//  //     muCand.trackerHits = (*data_.mu_trackerHits)[muIndex];
//      leptonCand_.push_back( muCand );
//    }
//
//    scale = 0.;
//    for( int l = 0; l < data_.mu_N; ++l ){
//       if( muIndex == 999 || l == muIndex ) continue;
//       // if( !isLowMass && (*data_.mu_isHighPtMuon)[l] != 1 ) continue;
//       // if(  isLowMass && (*data_.mu_isLooseMuon)[l] != 1 ) continue;
//       if( (*data_.mu_isLooseMuon)[l] != 1 ) continue;
//       scale = getMuonScale((*data_.mu_pt)[l]);
//       newpt = (*data_.mu_pt)[l]+scale;
//       if( (*data_.mu_trackIso)[l]/newpt >= 0.1 ) continue;
//       if( newpt <= AleptPtCut_ ) continue;
//       if( fabs((*data_.mu_eta)[l]) >= AleptEtaCut_ ) continue;
//       foundAMu++;
//    }
//
//    for( int l = 0; l < data_.el_N; ++l ){
//       // if( !isLowMass && (*data_.el_isHEEP)[l] != 1 ) continue;
//       // if( isLowMass && (*data_.el_isVetoElectron)[l] != 1 ) continue;
//       if( (*data_.el_isVetoElectron)[l] != 1 ) continue;
//       if( (*data_.el_pt)[l] <= 30. ) continue;
//       if( fabs((*data_.el_eta)[l]) >= 1.4442 && fabs((*data_.el_eta)[l]) <= 1.566 ) continue;
//       if( fabs((*data_.el_eta)[l]) >= 2.5 ) continue;
//       foundAEle++;
//    }
//
//    return foundLept;
//
//    //if( (foundAMu+foundAEle)<1 && foundLept ) return true;
//    //else return false;
//
// }
// //==============================================================================================
// bool ExoDiBosonAnalysis::findElectronCandidate( void ){
//
//    foundAMu = 0;
//    foundAEle = 0;
//
//    TLorentzVector TLV;
//    bool foundEle = false;
//    double ptEle = -999;
//    int eleIndex = 999;
//    for( int l = 0; l < data_.el_N; ++l ){
//       if( (*data_.el_pt)[l] <= leptPtCut_ ) continue;
//       if( fabs((*data_.el_eta)[l]) >= 1.4442 && fabs((*data_.el_eta)[l]) <= 1.566 ) continue;
//       if( fabs((*data_.el_eta)[l]) >= leptEtaCut_ ) continue;
//       if( (*data_.el_isHEEP)[l] != 1 ) continue;
//       // if(  isLowMass && (*data_.el_isTightElectron)[l] != 1 ) continue;
//
//       foundEle = true;
//       if( (*data_.el_pt)[l] > ptEle ){
//          ptEle = (*data_.el_pt)[l];
//          TLV.SetPtEtaPhiE( (*data_.el_pt)[l], (*data_.el_eta)[l], (*data_.el_phi)[l], (*data_.el_e)[l] );
//          eleIndex = l;
//       }
//    }
//
//    if( foundEle )
//     leptonCand_.push_back( LeptonCandidate( TLV ) );
//
//    for( int l = 0; l < data_.el_N; ++l ){
//       if( eleIndex == 999 || l == eleIndex ) continue;
//       // if( !isLowMass && (*data_.el_isHEEP)[l] != 1 ) continue;
//       if( (*data_.el_isVetoElectron)[l] != 1 ) continue;
//       if( (*data_.el_pt)[l] <= 30. ) continue;
//       if( fabs((*data_.el_eta)[l]) >= 1.4442 && fabs((*data_.el_eta)[l]) <= 1.566 ) continue;
//       if( fabs((*data_.el_eta)[l]) >= 2.5 ) continue;
//       foundAEle++;
//    }
//    for( int l = 0; l < data_.mu_N; ++l ){
//      // if( !isLowMass && (*data_.mu_isHighPtMuon)[l] != 1 ) continue;
//      if(  (*data_.mu_isLooseMuon)[l] != 1 ) continue;
//       if( (*data_.mu_trackIso)[l]/(*data_.mu_pt)[l] >= 0.1 ) continue;
//       if( (*data_.mu_pt)[l] <= 20. ) continue;
//       if( fabs((*data_.mu_eta)[l]) >= 2.4 ) continue;
//       foundAMu++;
//    }
//    return foundEle;
//  }
//
// //==============================================================================================
// void ExoDiBosonAnalysis::findAK4Jets( void ){
//   TLorentzVector jet;
//   for( int j = 0; j < data_.njetsAK4; ++j ){
//      jet.SetPtEtaPhiE( (*data_.jetAK4_pt)[j], (*data_.jetAK4_eta)[j], (*data_.jetAK4_phi)[j], (*data_.jetAK4_e)[j] );
//      if( jet.DeltaR(Vcand[0].p4) < 0.8 ) continue;
//      if( (*data_.jetAK4_IDLoose)[j] != 1 ) continue;
//      if( (*data_.jetAK4_pt)[j] <= JetPtCutLoose_ ) continue;
//      if( fabs((*data_.jetAK4_eta)[j]) >= JetEtaCut_ ) continue;
//      if( jet.DeltaR(leptonCand_[0].p4) < 0.3 ) continue;
//      JetCandidate Ajet(jet);
//      Ajet.csv = (*data_.jetAK4_csv)[j];
//      // Ajet.flavor = 0;
//      // if( runOnMC_ ) Ajet.flavor = abs((*data_.jetAK4_flavor)[j]);
//      AK4jetCand_.push_back(Ajet);
//   }
// }
//
// //==============================================================================================
// bool ExoDiBosonAnalysis::findMETCandidate( void ){
//
//   TLorentzVector TLV;
//   bool foundMET = false;
//   if( (*data_.MET_et)[0] > 0. ){
//
//      TLV.SetPtEtaPhiE( (*data_.MET_et)[0], 0., (*data_.MET_phi)[0], 0. );
//
//      METCandidate metC(TLV);
//      metC.sumEt = (*data_.MET_sumEt)[0];
//      METCand_.push_back( METCandidate( TLV ) );
//      foundMET = true;
//
//   }
//
//   return foundMET;
//
// }
//
// //==============================================================================================
// bool ExoDiBosonAnalysis::findWCandidate( void ){
//
//   TLorentzVector TLV;
//   bool foundW = false;
//   if( (METCand_[0].p4+leptonCand_[0].p4).Pt() > 0. ){
//      foundW = true;
//      TLV = METCand_[0].p4+leptonCand_[0].p4;
//      WCand_.push_back( TLV );
//   }
//
//   return foundW;
// }
// //==============================================================================================
// bool ExoDiBosonAnalysis::findHadronicWCandidate( void ){
//
//   TLorentzVector TLV,TLV1, TLV2;
//   bool foundW = false;
//   Double_t deltaWmassMin = 999.;
//
//   for(int i=0; i < abs(AK4jetCand_.size()); i++ ){
//     if (AK4jetCand_.at(i).p4.Pt() == Leptonicb_.at(0).p4.Pt() || AK4jetCand_.at(i).p4.Pt() <= 0. ) continue;
//     for(int j=0; j != i && j < abs(AK4jetCand_.size()); j++ ){
//       Double_t deltaWmass = fabs(pdgMw - ((AK4jetCand_.at(i).p4+ AK4jetCand_.at(j).p4).M()));
//       if( deltaWmass < deltaWmassMin ){
//         deltaWmassMin = deltaWmass;
//         foundW = true;
//         TLV = AK4jetCand_.at(i).p4+AK4jetCand_.at(j).p4;
//         TLV1 = AK4jetCand_.at(i).p4;
//         TLV2 = AK4jetCand_.at(j).p4;
//       }
//     }
//   }
//   HadronicJ_.push_back( TLV1 );
//   HadronicJ_.push_back( TLV2 );
//   Hist( "DeltaRJJ")->Fill ( HadronicJ_.at(1).p4.DeltaR(HadronicJ_.at(0).p4) , weight_  );
//
//   WCand_.push_back( TLV );
//   return foundW;
//
// }

