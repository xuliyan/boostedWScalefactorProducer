import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from xsec import getXsec

class selectionProducer(Module):
    def __init__(self, jetSelection):
        self.jetSel = jetSelection
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
		self.out = wrappedOutputTree
		self.out.branch("dr_LepJet",  "F")
		self.out.branch("dphi_LepJet",  "F")
		self.out.branch("dphi_MetJet",  "F")
		self.out.branch("dphi_WJet"  ,  "F")
		self.out.branch("FatJet_isW",  "F")
		self.out.branch("FatJet_softDrop_mass",  "F")
		self.out.branch("FatJet_tau21",  "F")
		self.out.branch("FatJet_tau21_ddt",  "F")
		self.out.branch("FatJet_tau21_ddt_retune",  "F")
		self.out.branch("W_type",  "F")
		self.out.branch("W_pt",  "F")
		self.out.branch("MET",  "F")
		self.out.branch("xsec",  "F")
		xsec = getXsec(inputFile.GetName())
		print inputFile.GetName()
		print xsec
		self.out.fillBranch("xsec",xsec)
		
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
        isMC = event.run == 1
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        jets = list(Collection(event, "Jet"))
        fatjets = list(Collection(event, "FatJet"))
        met = Object(event, "MET")
        triggerMu = Object(event, "HLT_Mu50")
        triggerEl = Object(event, "HLT_Ele115_CaloIdVT_GsfTrkIdT")
        	
			
		# Fat jet selection
        wFatJets =  [x for x in fatjets if x.pt>200 and abs(x.eta)<2.5]
        if len(wFatJets) < 1 : return False
        wFatJets.sort(key=lambda x:x.pt,reverse=True)
        
		# Medium b-tag
        passedjets = [x for x in jets if x.jetId>0 and x.pt>35 and abs(x.eta)<2.5 and x.btagCMVA> 0.89]
        if len(passedjets) < 1 : return False

		#lepton selection
        wElectrons = [x for x in electrons if x.cutBased_HEEP and x.pt > 35 ]	 #loose pt cut for veto 
        wMuons = [x for x in muons if x.pt > 20 and x.highPtId >= 1 ]   			 #loose pt cut for veto
        wMuons.sort(key=lambda x:x.pt,reverse=True)
        wElectrons.sort(key=lambda x:x.pt,reverse=True)
        
        # for j in filter(self.jetSel,fatjets):
        Vtype = -1
        vLeptons = [] # decay products of V
        if len(wElectrons) + len(wMuons) == 1:
            if len(wMuons) == 1:
                if wMuons[0].pt < 53: return False
                Vtype = 0
                vLeptons = [wMuons[0]]
                if triggerMu == 0: return False
            if len(wElectrons) == 1:
                if wElectrons[0].pt < 115: return False
                Vtype=1
                vLeptons = [wElectrons[0]]
                if triggerEl == 0: return False
        else: return False
        
        vLepton_4vec = ROOT.TLorentzVector()
        vLepton_4vec.SetPtEtaPhiM(vLeptons[0].pt,vLeptons[0].eta,vLeptons[0].phi,vLeptons[0].mass)
        
        jet_4v = ROOT.TLorentzVector()
        jet_4v.SetPtEtaPhiM(wFatJets[0].pt,wFatJets[0].eta,wFatJets[0].phi,wFatJets[0].mass)
        dR_jetlep = jet_4v.DeltaR(vLepton_4vec)
        
        if dR_jetlep < 1.0: return False
		
		
        isW = 0
        if isMC:
			gens = Collection(event, "GenPart")
			daus =  [x for x in gens if x.pt>1 and 0<abs(x.pdgId)<9]
			moms =  [x for x in gens if x.pt>10 and abs(x.pdgId)==24]
        
			realVs = []
			if len(moms)>0 and len(daus)>0:
				for dau in daus:
					for mom in moms:
						try:
							if mom == moms[dau.genPartIdxMother]: realVs.append(mom)	
						except:
							continue
			for V in realVs:
				gen_4v = ROOT.TLorentzVector()
				gen_4v.SetPtEtaPhiM(V.pt,V.eta,V.phi,80.)
				dR = jet_4v.DeltaR(gen_4v)
				if dR < 0.8: isW = 1
        	
        
        met_4v  = ROOT.TLorentzVector()
        met_4v.SetPtEtaPhiE(met.sumEt, 0., met.phi, met.sumEt)
        if met.sumEt < 40: return False
        ## add branches for some basic V kinematics
        V = ROOT.TLorentzVector()
        for vLepton in vLeptons:
            vLeptons_4vec = ROOT.TLorentzVector()
            vLeptons_4vec.SetPtEtaPhiM(vLepton.pt,vLepton.eta,vLepton.phi,vLepton.mass)
            V = V + vLeptons_4vec
        self.out.fillBranch("dr_LepJet"  ,dR_jetlep)
        self.out.fillBranch("dphi_LepJet",jet_4v.DeltaPhi(vLepton_4vec))
        self.out.fillBranch("dphi_MetJet",jet_4v.DeltaPhi(met_4v))
        self.out.fillBranch("dphi_WJet"  ,jet_4v.DeltaPhi(V))
        self.out.fillBranch("W_type",Vtype)
        self.out.fillBranch("W_pt", V.Perp()+met.sumEt )
        self.out.fillBranch("MET", met.sumEt )
        self.out.fillBranch("FatJet_isW", isW)
        self.out.fillBranch("FatJet_softDrop_mass",  wFatJets[0].msoftdrop)
        self.out.fillBranch("FatJet_tau21", wFatJets[0].tau2/wFatJets[0].tau1)
        self.out.fillBranch("FatJet_tau21_ddt", wFatJets[0].tau2/wFatJets[0].tau1+0.063*ROOT.TMath.Log(wFatJets[0].msoftdrop**2/wFatJets[0].pt))
        self.out.fillBranch("FatJet_tau21_ddt_retune", wFatJets[0].tau2/wFatJets[0].tau1+0.082*ROOT.TMath.Log(wFatJets[0].msoftdrop**2/wFatJets[0].pt))
        
        return True
        

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

selectionModule = lambda : selectionProducer(jetSelection= lambda j : j.pt > 200) 
 
