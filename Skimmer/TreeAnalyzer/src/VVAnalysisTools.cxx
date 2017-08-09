#include "../include/VVAnalysisTools.h"
#include "../include/VVanalysis.h"


bool isMergedVJet(TLorentzVector goodFatJet, std::vector<UZH::GenParticle> VdecayProducts) {
    int associatedQuarks=0;
    for(unsigned int i=0;i< VdecayProducts.size();i++)
    {
      TLorentzVector q = VdecayProducts[i].tlv();
      if (goodFatJet.DeltaR(q) < 0.8) associatedQuarks +=1;

    }
    if (associatedQuarks ==2) return 1;
    else return 0;   
}


std::vector<UZH::GenParticle> FindGeneratedQuarks(Ntuple::GenParticleNtupleObject m_genParticle,bool m_isData)
{
  std::vector<UZH::GenParticle> GenQuarks;
  if(m_isData) return GenQuarks;
  for(int i=0;i< m_genParticle.N;i++)
  {
    bool containsV=0;
    if(TMath::Abs(m_genParticle.pdgId->at(i)) > 6 or TMath::Abs(m_genParticle.pdgId->at(i)) <1) continue;
       
    for(int ii=0;ii< m_genParticle.nMoth->at(i);ii++)
    {
           
      if(TMath::Abs(m_genParticle.mother->at(i).at(ii)) == 24 or TMath::Abs(m_genParticle.mother->at(i).at(ii)) == 23) containsV = 1;  

    }
    if (!containsV) continue;
    UZH::GenParticle MyQuark( &m_genParticle, i );
    GenQuarks.push_back(MyQuark);
  }
  return GenQuarks;    
}


float ApplyPuppiSoftdropMassCorrections(UZH::Jet puppiJet,std::vector<TF1*> m_puppisd_corr, bool m_isData){
 float genCorr =1;
 if(!m_isData) genCorr = m_puppisd_corr[0]->Eval(puppiJet.pt());
 float recoCorr = 1;
 if( fabs(puppiJet.eta()) <= 1.3) recoCorr = m_puppisd_corr[1]->Eval(puppiJet.pt());
 else if (fabs(puppiJet.eta()) > 1.3) recoCorr = m_puppisd_corr[2]->Eval(puppiJet.pt());
    
 return puppiJet.softdrop_massCorr()*genCorr*recoCorr;
}
