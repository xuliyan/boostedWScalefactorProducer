import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import *
from xsec import getXsec

import random
import array

class TTbar_SemiLep(Module):
    def __init__(self ):
        self.writeHistFile = True
        self.verbose = False
    def beginJob(self, histFile, histDirName):
        Module.beginJob(self, histFile, histDirName)

        self.isttbar = False
        if 'TTJets_' in histFile :
           self.isttbar = True 



        ### Set bins for Pt dependent scale factor calculation    

        self.TopcandPtBins = [[200,300], [300,400], [400,500],[500, -1]]
        # e.g. h_wpt_Tptbin0 is a 1D histogram for W candidate subjets (most massive SD subjet) within Top candidates of pt 200-300 GeV  

        self.WcandPtBins = [[200,300], [300-500], [500, -1]]
        # e.g. h_wpt_ptbin0 1D histogram is for W candidate subjets (most massive SD subjet) with pt 200-300 GeV   

        self.minMupt = 53.
        self.maxMuEta = 2.4
        self.maxRelIso = 0.1
        self.minMuMETPt = 40.

        ### Figure our tree branch for HighPtMuon ???
        #is High Pt


        #remove  AK8 jet within 1.0 of lepton
        self.mindRLepJet = 1.0 
        #veto:
        # High pT muon ID
        #pT > 20 GeV, eta < 2.4??
        #relIso < 0.1


        self.minElpt = 120.
        self.minElMETPt = 80.
        #self.goodElEta = if eta < 1.44, 1.56 < eta < 2.5
        # HEEP v7 + iso
        #veto
        # HEEP + iso pt > 35 remove ecal crack region eta < 1.44, 1.56 < eta < 2.5
        #

        self.minLepWPt = 150.

        self.minJetPt = 200.
        self.maxJetEta = 2.5

        self.minBDisc = 0.8484
        ### Medium https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco

        #>= 1 CSVmedium akt4 jet
        self.minAK4Pt = 30.

        #Angular selection (not used by Thea now):

        #dR( lepton, leading AK8 jet) > pi/2
        #dPhi(leading AK8 jet, MET) > 2
        
        #dPhi (leading AK8 jet, leptonic W) >2
        #self.minDPhiWJet = 2.  


                  
    def endJob(self):
        Module.endJob(self)
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        # This is a loose selection to select events for T Tbar semileptonic events where 
        # Type 1 and Type 2 events are included in the selection:

        # In SF code the final cuts need to be made to choose either type 1 or type 2 selection:
        # e.g. for type 1 the W leptonic Pt cut should be tightened to 200 GeV and dPhi cuts applied
        # e.g. for type 2 the AK8 Pt cut should be tightened to 400 GeV and dPhi cuts applied

        # Type 1 - martially merged Hadronic Top Quark (W is AK8, b is AK4) 
        #(AK8 Pt > 200 GeV)

        # Type 2 - fully merged Top (Top is AK8, W is most massive SD subjet, b is less massive subjet, require 1 subjet b-tag) 
        #(AK8 Pt > 400 GeV): 


        # selection aligned with previous SF measurement standard selection
        # https://www.evernote.com/shard/s282/sh/7e5d6baa-d100-4025-8bf8-a61bf1adfbc1/f7e86fde2c2a165e
        

        # 1 AK8 Pt > 200 GeV, |eta| < 2.5 , dR(Ak8, lep) > 1.0
        # 1 AK4 Pt > 30 GeV, |eta| < 2.5
        # 1 lepton , mu pt > 53 GeV or el pt > 120 GeV
        # MET Pt > 40(mu) or 80(el) GeV
        #Leptonic W - lepton + MET has Pt > 150 GeV # did not apply this since we are missing MET eta

        self.out.branch("xsec",  "F")
        self.out.branch("genmatchedAK8Subjet",  "F")
        self.out.branch("AK8Subjet0isMoreMassive",  "F")
        self.out.branch("genmatchedAK8",  "F")

        self.xs = getXsec(inputFile.GetName())
        print inputFile.GetName()
        print self.xs
        

        pass

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def getSubjets(self, p4, subjets, dRmax=0.8):
        ret = []
        for subjet in subjets :
            if p4.DeltaR(subjet.p4()) < dRmax and len(ret) < 2 :
                ret.append(subjet.p4())
        return ret

    def printP4( self, c ):
        if hasattr( c, "p4"):
            s = ' %6.2f %5.2f %5.2f %6.2f ' % ( c.p4().Perp(), c.p4().Eta(), c.p4().Phi(), c.p4().M() )
        else :
            s = ' %6.2f %5.2f %5.2f %6.2f ' % ( c.Perp(), c.Eta(), c.Phi(), c.M() )
        return s
    def printCollection(self,coll):
        for ic,c in enumerate(coll):
            s = self.printP4( c )
            print ' %3d : %s' % ( ic, s )
            
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        weight = 1.0

        isMC = event.run == 1
        if self.verbose:
            print '------------------------ ', event.event

        if isMC:

            ### Look at generator level particles
            ### find events where :
            ### a W decays to quarks (Type 1 - partially merged)
            ###    OR
            ### a Top decays to W + b (Type 2 - fully merged top quark)
            gens = Collection(event, "GenPart")
            Wdaus =  [x for x in gens if x.pt>1 and 0<abs(x.pdgId)<9]
            Wmoms =  [x for x in gens if x.pt>10 and abs(x.pdgId)==24]

            TWdaus =  [x for x in gens if x.pt>1 and  0<abs(x.pdgId)<4]
            Tdaus =  [x for x in gens if x.pt>1 and (abs(x.pdgId)==5  or  abs(x.pdgId)==24 )]
            Tmoms =  [x for x in gens if x.pt>10 and abs(x.pdgId)==6]
                    
            realVs = []

            realTs = []
            realWs = []
            realqs = []
            self.matchedJ = 0
            self.matchedSJ = 0

            if len(Wmoms)>0 and len(Wdaus)>0:
                for dau in Wdaus:
                    for mom in Wmoms:
                        try:
                            if mom == Wmoms[dau.genPartIdxMother]: realVs.append(mom)    
                        except:
                            continue    

            if len(Tmoms)>0 and len(Tdaus)>0:
                for gdau in TWdaus :
                    for dau in Tdaus:
                        for mom in Tmoms:
                            try:
                                if mom == Tmoms[dau.genPartIdxMother] and dau == Tdaus[gdau.genPartIdxMother]: 
                                    realTs.append(mom)
                                    realWs.append(dau)
                                    realqs.append(gdau)    
                            except:
                                continue  


            ###### Get gen Top candidate #######
            genleptons = Collection(event, "GenDressedLepton")

            if len(genleptons) < 1 :
                return False
            if abs(genleptons[0].pdgId) != 13 and abs(genleptons[0].pdgId)!= 12 :
                return False
            ########################

            ### Gen Selection
            ########################

            ### We want the AK4 nearest the lepton
            ### This is the b candidate
            ### check b disc in tagging script


            genAK4s = Collection(event, "GenJet")    
            genAK4jets = [ x for x in genAK4s if x.p4().Perp() > self.minAK4Pt * 0.8 and abs(x.p4().Eta()) < self.maxJetEta ]
            

            METgen_pt = event.GenMET_pt

           
            ### Gen Electron
            ### Pick the high Pt (self.minElpt * 0.8), low eta, HEEPv7 Electrons    
            if abs(genleptons[0].pdgId) == 12 :
                if genleptons[0].p4().Perp() < self.minElpt * 0.8 :
                    return False
                if abs(genleptons[0].p4().Eta()) < 1.44 or (abs(genleptons[0].p4().Eta()) > 1.56 and  abs(genleptons[0].p4().Eta()) < 2.5  ) :
                    return False 
                if METgen_pt < self.minElMETPt :
                    return False
            ### Gen Muons  
            ### Pick high pt, low eta muons       
            if abs(genleptons[0].pdgId) == 13 :
                if genleptons[0].p4().Perp() < self.minMupt * 0.8 :
                    return False
                if abs(genleptons[0].p4().Eta())  > self.maxMuEta :
                    return False 
                if METgen_pt < self.minMuMETPt :
                    return False
            if self.verbose :
                print '----'
                print 'Gen leptons:'
                self.printCollection( genleptons )


            # we want the ak4 closest to the lepton
            #mindRObs = 5.0
            #genbHad = ROOT.TLorentzVector()
            #for ibcand, bcand in enumerate(genAK4jets) :
            #    tempdR = bcand.p4().DeltaR(genleptons[0].p4())
            #    if  tempdR < mindRObs :
            #        mindRObs = tempdR
            #        genbHad.SetPtEtaPhiM( bcand.p4().Perp(), bcand.p4().Eta() , bcand.p4().Phi() , bcand.p4().M()  )

            
            METgen_phi = event.GenMET_phi

            METgen = ROOT.TLorentzVector()

            METgen.SetPtEtaPhiM(METgen_pt,0.0,METgen_phi, 0.0  )
            WLep = genleptons[0].p4() + METgen
            if WLep.Perp() < self.minLepWPt * 0.9 :
                return False

            if self.verbose:
                print '-----'
                print 'Gen W Leptonic:'
                print self.printP4( WLep )

    

            ###### Get list of gen jets #######
            # List of gen jets:
            allgenjets = list(Collection(event, "GenJetAK8"))
            if self.verbose:
                print '-----'
                print 'all genjets:'
                self.printCollection( allgenjets )
            genjets = [ x for x in allgenjets if x.p4().Perp() > self.minJetPt * 0.8 and abs( x.p4().Eta()) < self.maxJetEta]

            # List of gen subjets (no direct link from Genjet):
            #gensubjets = list(Collection(event, "SubGenJetAK8"))
            # Dictionary to hold ungroomed-->groomed for gen
            #genjetsGroomed = {}
            # Get the groomed gen jets
            #maxSubjetMass = 1.
            '''
            WHad = ROOT.TLorentzVector()
            
            for igen,gen in enumerate(genjets):
                hasBtagSJ = None
                gensubjetsMatched = self.getSubjets( p4=gen.p4(),subjets=gensubjets, dRmax=0.8)
                for isub,sub in enumerate(gensubjetsMatched) : 
                    #if sub.btagCSVV2 > self.minBDisc :
                    #    hasBtagSJ = True
                    ### Note : will need to check for a b-tagged subjet in post selection
                    if sub.M() > maxSubjetMass : # and hasBtagSJ :
                        maxSubjetMass = sub.M() 
                        WHad.SetPtEtaPhiM(sub.Perp(),sub.Eta(),sub.Phi(),sub.M())
                        

                genjetsGroomed[gen] = sum( gensubjetsMatched, ROOT.TLorentzVector() ) if len(gensubjetsMatched) > 0   else None
                
            if self.verbose:
                print '----'
                print 'opposite-LepW genjets:'
                for genjet in genjets:
                    sdmassgen = genjetsGroomed[genjet].M() if genjet in genjetsGroomed else -1.0
                    print '         : %s %6.2f' % ( self.printP4(genjet), sdmassgen )            
            
            ''' 
            
        ###### Get reco Top/W candidate #######
        # List of reco muons
        allmuons = Collection(event, "Muon")
        allelectrons = Collection(event, "Electron")
        # Select reco muons:
        muons = [ x for x in allmuons if x.tightId and x.pfRelIso03_all and x.p4().Perp() > self.minMupt and abs(x.p4().Eta())  < self.maxMuEta ]
        # Select reco muons:
        electrons = [ x for x in allelectrons if  x.cutBased_HEEP and x.p4().Perp() > self.minElpt  and (abs(x.p4().Eta()) < 1.44 or (abs(x.p4().Eta()) > 1.56 and  abs(x.p4().Eta()) < 2.5  ))]
        
        print "N leptons =?? " ,len(muons) + len(electrons)
        if ( len(muons) + len(electrons) ) < 1 :
            return False

        ### Choose the leading lepton in the event
        lepton = ROOT.TLorentzVector()
        isMu = None #False

        # Keep only events with exactly 1 lepton passing all above pt , eta and cut based ID cuts

        if  ( len(muons) ) == 1 and  ( len(electrons) ) < 1  :
            lepton = muons[0].p4()
            isMu = True

        if  ( len(muons) ) < 1 and  ( len(electrons) ) == 1  :
            lepton = electrons[0].p4()
            isMu = False
      
        if  ( len(muons) ) > 0 and  ( len(electrons) ) > 0  :
            # Ignore events with muons and electrons
            return False


        MET_pt = event.PuppiMET_pt     
        
        if isMu and MET_pt < self.minMuMETPt :
            return False
        if not isMu and MET_pt < self.minElMETPt :
            return False


        MET = ROOT.TLorentzVector()

        MET.SetPtEtaPhiM(MET_pt, 0.0, event.PuppiMET_phi , event.PuppiMET_sumEt)

        WcandLep = lepton + MET
        if WcandLep.Perp() < self.minLepWPt :
            return False

        allrecoAK4jets = list(Collection(event, "Jet")) # are these AK4s ? 
        recojetsAK4 = [ x for x in allrecoAK4jets if x.p4().Perp() > self.minAK4Pt and abs(x.p4().Eta()) < self.maxJetEta]
        if len(recojetsAK4) < 1:
            return False
        '''
        mindRObs = 5.0
        bHadreco = ROOT.TLorentzVector()
        for ibcand, bcand in enumerate(recojetsAK4 ) :
            tempdR = bcand.p4().DeltaR(genleptons[0].p4())
            if  tempdR < mindRObs :
                mindRObs = tempdR
                bHadreco.SetPtEtaPhiM( bcand.p4().Perp(), bcand.p4().Eta() , bcand.p4().Phi() , bcand.p4().M()  )
        '''


        #Topcandreco =  WcandLep +  bHadreco

        #if self.verbose:
        #    print '-----'
        #    print ' reco Top Leptonic:', self.printP4( Topcandreco)
        
        ###### Get list of reco jets #######
        # List of reco jets:
        allrecojets = list(Collection(event, "FatJet"))
        if self.verbose:
            print '----'
            print 'all recojets:'
            self.printCollection( allrecojets )
        recojets = [ x for x in allrecojets if x.p4().Perp() > self.minJetPt and  abs(x.p4().Eta()) < self.maxJetEta ]
        if len(recojets) < 1 : return False
        recojets.sort(key=lambda x:x.pt,reverse=True)

        jet_4v = ROOT.TLorentzVector()
        jet_4v.SetPtEtaPhiM(recojets[0].pt,recojets[0].eta,recojets[0].phi,recojets[0].mass)
        dR_jetlep = jet_4v.DeltaR(lepton )
        
        if dR_jetlep < self.mindRLepJet : return False

        self.isW = 0
        self.SJ0isW = -1
        if isMC == False:
            genjets = [None] * len(recojets)

        else :
                
            for V in realVs:
                gen_4v = ROOT.TLorentzVector()
                gen_4v.SetPtEtaPhiM(V.pt,V.eta,V.phi,80.)
                dR = jet_4v.DeltaR(gen_4v)
                if dR < 0.8: self.isW = 1
    
        # List of reco subjets:
        recosubjets = list(Collection(event,"SubJet"))
        # Dictionary to hold ungroomed-->groomed for reco
        recojetsGroomed = {}        
        # Get the groomed reco jets
        maxrecoSJmass = 1.
        WHadreco = ROOT.TLorentzVector()
        for ireco,reco in enumerate(recojets):
            if reco.subJetIdx2 >= len(recosubjets) or reco.subJetIdx1 >= len(recosubjets) :
                if self.verbose: print "Reco subjet indices not in Subjet list, Skipping"
                continue
            if reco.subJetIdx1 >= 0 and reco.subJetIdx2 >= 0 :
                recojetsGroomed[reco] = recosubjets[reco.subJetIdx1].p4() + recosubjets[reco.subJetIdx2].p4()
                if recosubjets[reco.subJetIdx1].p4().M() > maxrecoSJmass and recosubjets[reco.subJetIdx1].p4().M() >  recosubjets[reco.subJetIdx2].p4().M() :
                    maxrecoSJmass = recosubjets[reco.subJetIdx1].p4().M() 
                    WHadreco = recosubjets[reco.subJetIdx1].p4()
                    if recosubjets[reco.subJetIdx1].btagCSVV2 >  self.minBDisc  or recosubjets[reco.subJetIdx2].btagCSVV2 >  self.minBDisc :
                        self.SJ0isW = 1
                if recosubjets[reco.subJetIdx2].p4().M() > maxrecoSJmass and recosubjets[reco.subJetIdx1].p4().M() < recosubjets[reco.subJetIdx2].p4().M() :
                    maxrecoSJmass = recosubjets[reco.subJetIdx1].p4().M() 
                    WHadreco = recosubjets[reco.subJetIdx2].p4()
                    if recosubjets[reco.subJetIdx1].btagCSVV2 >  self.minBDisc  or recosubjets[reco.subJetIdx2].btagCSVV2 >  self.minBDisc :
                        self.SJ0isW = 0
                if isMC :
                    for q in realqs:
                        gen_4v = ROOT.TLorentzVector()
                        gen_4v.SetPtEtaPhiM(q.pt,q.eta,q.phi,q.mass)
                        dR = WHadreco.DeltaR(gen_4v)
                        if dR < 0.6 and self.isttbar : self.matchedSJ = 1  
            elif reco.subJetIdx1 >= 0 :
                recojetsGroomed[reco] = recosubjets[reco.subJetIdx1].p4()
                maxrecoSJmass = recosubjets[reco.subJetIdx1].p4().M() 
                WHadreco = recosubjets[reco.subJetIdx1].p4()     

            else :
                recojetsGroomed[reco] = None
                WHadreco = None

        if self.verbose:
            print '----'
            print ' recojets opposite the lepton  (Top/W candidates) :'
            for recojet in recojets:
                sdmassreco = recojetsGroomed[recojet].M() if recojet in recojetsGroomed and recojetsGroomed[recojet] != None else -1.0
                print '         : %s %6.2f' % ( self.printP4( recojet),  sdmassreco )            
        if isMC :
            self.out.fillBranch("genmatchedAK8Subjet", self.matchedSJ)
            self.out.fillBranch("genmatchedAK8",  self.isW)
        self.out.fillBranch("AK8Subjet0isMoreMassive", self.SJ0isW )
        self.out.fillBranch("xsec",self.xs)
        


        return True
# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

ttbar_semilep = lambda : TTbar_SemiLep() 
