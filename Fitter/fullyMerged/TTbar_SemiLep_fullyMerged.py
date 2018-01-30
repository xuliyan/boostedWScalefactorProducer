import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import *


import random
import array

class TTbar_SemiLep_fullyMerged(Module):
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
        # e.g. h_WcandSubjetpt_Tptbin0 is a 1D histogram for W candidate subjets (most massive SD subjet) within Top candidates of pt 200-300 GeV  

        self.WcandPtBins = [[200,300], [300-500], [500, -1]]
        # e.g. h_WcandSubjetpt_ptbin0 1D histogram is for W candidate subjets (most massive SD subjet) with pt 200-300 GeV   

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

        self.minJetPt = 350.
        self.maxJetEta = 2.4
        self.minTopmass = 105.
        self.maxTopmass = 250.
        self.maxtau32Top = 0.7



        self.minBDisc = 0.8484
        ### Medium https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco

        #>= 1 CSVmedium akt4 jet
        self.minAK4Pt = 50.

        ### 2-D cut ###
        ### dR OR PtRel ###
        self.mindRlepAK4 = 0.4
        self.minPtRel_lepAK4 = 30.

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

       
        
        
        self.out.branch("WHadreco_pt"   ,  "F")
        self.out.branch("WHadreco_eta"  ,  "F")
        self.out.branch("WHadreco_phi"  ,  "F")
        self.out.branch("WHadreco_mass" ,  "F")
        self.out.branch("WHadreco_tau21" ,  "F")

        
        self.addObject( ROOT.TH1D('h_lep0pt',          'h_lep0pt',        40, 0, 200 ) )
        self.addObject( ROOT.TH1D('h_lep0eta',         'h_lep0eta',      48, -3, 3 ) )
        self.addObject( ROOT.TH1D('h_lep0phi',         'h_lep0phi',      100, -5, 5 ) )

        self.addObject( ROOT.TH1D('h_hadToppt',          'h_hadToppt',        100, 0, 500 ) )
        self.addObject( ROOT.TH1D('h_hadTopeta',         'h_hadTopeta',      48, -3, 3 ) )
        self.addObject( ROOT.TH1D('h_hadTopphi',         'h_hadTopphi',      100, -5, 5 ) )
        self.addObject( ROOT.TH1D('h_hadTopmass',        'h_hadTopmass',      60, 140, 200 ) )
        '''
        self.addObject( ROOT.TH1D('h_lepToppt',          'h_lepToppt',        100, 0, 500 ) )
        self.addObject( ROOT.TH1D('h_lepTopeta',         'h_lepTopeta',      48, -3, 3 ) )
        self.addObject( ROOT.TH1D('h_lepTopphi',         'h_lepTopphi',      100, -5, 5 ) )
        self.addObject( ROOT.TH1D('h_lepTopmass',        'h_lepTopmass',      60, 140, 200 ) )
        '''
        self.addObject( ROOT.TH1D('h_WcandSubjetpt',          'h_WcandSubjetpt',        100, 0, 500 ) )
        self.addObject( ROOT.TH1D('h_WcandSubjeteta',         'h_WcandSubjeteta',      48, -3, 3 ) )
        self.addObject( ROOT.TH1D('h_WcandSubjetphi',         'h_WcandSubjetphi',      100, -5, 5 ) )
        self.addObject( ROOT.TH1D('h_WcandSubjetmass',        'h_WcandSubjetmass',      100, 50, 150 ) )
        
        self.addObject( ROOT.TH1D('h_WcandSubjetpt_ptbin0',          'h_WcandSubjetpt_ptbin0',        100, 0, 500 ) )
       
        self.addObject( ROOT.TH1D('h_WcandSubjeteta_ptbin0',         'h_WcandSubjeteta_ptbin0',      48, -3, 3 ) )
        self.addObject( ROOT.TH1D('h_WcandSubjetphi_ptbin0',         'h_WcandSubjetphi_ptbin0',      100, -5, 5 ) )
        self.addObject( ROOT.TH1D('h_WcandSubjetmass_ptbin0',        'h_WcandSubjetmass_ptbin0',      100, 50, 150 ) )

        self.addObject( ROOT.TH1D('h_WcandSubjetpt_ptbin1',          'h_WcandSubjetpt_ptbin1',        100, 0, 500 ) )
        self.addObject( ROOT.TH1D('h_WcandSubjeteta_ptbin1',         'h_WcandSubjeteta_ptbin1',      48, -3, 3 ) )
        self.addObject( ROOT.TH1D('h_WcandSubjetphi_ptbin1',         'h_WcandSubjetphi_ptbin1',      100, -5, 5 ) )
        self.addObject( ROOT.TH1D('h_WcandSubjetmass_ptbin1',        'h_WcandSubjetmass_ptbin1',      100, 50, 150 ) )

        self.addObject( ROOT.TH1D('h_WcandSubjetpt_ptbin2',          'h_WcandSubjetpt_ptbin2',        100, 0, 500 ) )
        self.addObject( ROOT.TH1D('h_WcandSubjeteta_ptbin2',         'h_WcandSubjeteta_ptbin2',      48, -3, 3 ) )
        self.addObject( ROOT.TH1D('h_WcandSubjetphi_ptbin2',         'h_WcandSubjetphi_ptbin2',      100, -5, 5 ) )
        self.addObject( ROOT.TH1D('h_WcandSubjetmass_ptbin2',        'h_WcandSubjetmass_ptbin2',      100, 50, 150 ) )
   

        self.WcandSubjetpt = [self.h_WcandSubjetpt_ptbin0, self.h_WcandSubjetpt_ptbin1, self.h_WcandSubjetpt_ptbin2 ]
        self.WcandSubjeteta = [self.h_WcandSubjeteta_ptbin0, self.h_WcandSubjeteta_ptbin1, self.h_WcandSubjeteta_ptbin2 ]
        self.WcandSubjetphi = [self.h_WcandSubjetphi_ptbin0, self.h_WcandSubjetphi_ptbin1, self.h_WcandSubjetphi_ptbin2 ]
        self.WcandSubjetmass = [self.h_WcandSubjetmass_ptbin0, self.h_WcandSubjetmass_ptbin1, self.h_WcandSubjetmass_ptbin2 ]

        '''
        self.addObject( ROOT.TH1D('h_WcandSubjetpt_Tptbin0',          'h_WcandSubjetpt_Tptbin0',        100, 0, 500 ) )
        self.addObject( ROOT.TH1D('h_WcandSubjeteta_Tptbin0',         'h_WcandSubjeteta_Tptbin0',      48, -3, 3 ) )
        self.addObject( ROOT.TH1D('h_WcandSubjetphi_Tptbin0',         'h_WcandSubjetphi_Tptbin0',      100, -5, 5 ) )
        self.addObject( ROOT.TH1D('h_WcandSubjetmass_Tptbin0',        'h_WcandSubjetmass_Tptbin0',      100, 50, 150 ) )

        self.addObject( ROOT.TH1D('h_WcandSubjetpt_Tptbin1',          'h_WcandSubjetpt_Tptbin1',        100, 0, 500 ) )
        self.addObject( ROOT.TH1D('h_WcandSubjeteta_Tptbin1',         'h_WcandSubjeteta_Tptbin1',      48, -3, 3 ) )
        self.addObject( ROOT.TH1D('h_WcandSubjetphi_Tptbin1',         'h_WcandSubjetphi_Tptbin1',      100, -5, 5 ) )
        self.addObject( ROOT.TH1D('h_WcandSubjetmass_Tptbin1',        'h_WcandSubjetmass_Tptbin1',      100, 50, 150 ) )

        self.addObject( ROOT.TH1D('h_WcandSubjetpt_Tptbin2',          'h_WcandSubjetpt_Tptbin2',        100, 0, 500 ) )
        self.addObject( ROOT.TH1D('h_WcandSubjeteta_Tptbin2',         'h_WcandSubjeteta_Tptbin2',      48, -3, 3 ) )
        self.addObject( ROOT.TH1D('h_WcandSubjetphi_Tptbin2',         'h_WcandSubjetphi_Tptbin2',      100, -5, 5 ) )
        self.addObject( ROOT.TH1D('h_WcandSubjetmass_Tptbin2',        'h_WcandSubjetmass_Tptbin2',      100, 50, 150 ) )

        self.addObject( ROOT.TH1D('h_genjetpt',          'h_genjetpt',   100, 0, 500 ) )
        self.addObject( ROOT.TH1D('h_genjeteta',         'h_genjeteta',      48, -3, 3 ) )
        self.addObject( ROOT.TH1D('h_genjetphi',         'h_genjetphi',      100, -5, 5 ) )
        self.addObject( ROOT.TH1D('h_genjetmass',        'h_genjetmass',      300, 0, 300 ) )
        '''

        self.addObject( ROOT.TH1D('h_Wleppt',          'h_Wleppt',        100, 0, 500 ) )
        self.addObject( ROOT.TH1D('h_Wlepeta',         'h_Wlepeta',      48, -3, 3 ) )
        self.addObject( ROOT.TH1D('h_Wlepphi',         'h_Wlepphi',      100, -5, 5 ) )
        self.addObject( ROOT.TH1D('h_Wlepmass',        'h_Wlepmass',      100, 50, 150 ) )



        if self.isttbar :
            self.addObject( ROOT.TH1D('h_matchedAK8Subjetpt',          'h_matchedAK8Subjetpt',      100, 0, 500 ) )
            self.addObject( ROOT.TH1D('h_matchedAK8Subjeteta',         'h_matchedAK8Subjeteta',      48, -3, 3 ) )
            self.addObject( ROOT.TH1D('h_matchedAK8Subjetphi',         'h_matchedAK8Subjetphi',      100, -5, 5 ) )
            self.addObject( ROOT.TH1D('h_matchedAK8Subjetmass',        'h_matchedAK8Subjetmass',      300, 0, 300 ) )

            self.addObject( ROOT.TH1D('h_unmatchedAK8Subjetpt',          'h_unmatchedAK8Subjetpt',      100, 0, 500 ) )
            self.addObject( ROOT.TH1D('h_unmatchedAK8Subjeteta',         'h_unmatchedAK8Subjeteta',      48, -3, 3 ) )
            self.addObject( ROOT.TH1D('h_unmatchedAK8Subjetphi',         'h_unmatchedAK8Subjetphi',      100, -5, 5 ) )
            self.addObject( ROOT.TH1D('h_unmatchedAK8Subjetmass',        'h_unmatchedAK8jetmass',      300, 0, 300 ) )



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

        ###### Get reco Top/W candidate #######
        # List of reco muons
        allmuons = Collection(event, "Muon")
        allelectrons = Collection(event, "Electron")
        # Select reco muons:
        muons = [ x for x in allmuons if x.tightId and x.pfRelIso03_all and x.p4().Perp() > self.minMupt and abs(x.p4().Eta())  < self.maxMuEta ]
        # Select reco muons:
        electrons = [ x for x in allelectrons if  x.cutBased_HEEP and x.p4().Perp() > self.minElpt  and (abs(x.p4().Eta()) < 1.44 or (abs(x.p4().Eta()) > 1.56 and  abs(x.p4().Eta()) < 2.5  ))]

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
        if len(recojetsAK4) < 1:  return False
        mindRObs = 5.0
        bHadreco = ROOT.TLorentzVector()
        for ibcand, bcand in enumerate(recojetsAK4 ) :
            tempdR = bcand.p4().DeltaR(lepton)
            ptrel = (bcand.p4() - lepton).Perp()
            #print ptrel
            if  tempdR < mindRObs and (tempdR >  self.mindRlepAK4 or abs(ptrel) > self.minPtRel_lepAK4 ):
                mindRObs = tempdR
                bHadreco.SetPtEtaPhiM( bcand.p4().Perp(), bcand.p4().Eta() , bcand.p4().Phi() , bcand.p4().M())  
        if self.verbose and bHadreco.Perp() > self.minAK4Pt :
            print '-----'
            print ' reco b candidate AK4:', self.printP4( bHadreco )
        
        ###### Get list of reco jets #######
        # List of reco jets:
        allrecojets = list(Collection(event, "FatJet"))
        if self.verbose:
            print '----'
            print 'all recojets:'
            self.printCollection( allrecojets )

        recojets = [ x for x in allrecojets if x.p4().Perp() > self.minJetPt and  abs(x.p4().Eta()) < self.maxJetEta  and x.p4().M() > self.minTopmass  and  x.p4().M() < self.maxTopmass]
        if len(recojets) < 1 : return False
        recojets.sort(key=lambda x:x.pt,reverse=True)


        
        if isMC == False:
            genjets = [None] * len(recojets)

    
        # List of reco subjets:
        recosubjets = list(Collection(event,"SubJet"))
        # Dictionary to hold reco--> gen matching
        #recoToGen = matchObjectCollection( recojets, genjets, dRmax=0.05 )
        # Dictionary to hold ungroomed-->groomed for reco
        recojetsGroomed = {}        
        # Get the groomed reco jets
        maxrecoSJmass = 1.
        WHadreco = ROOT.TLorentzVector()
        WHadrecoTau21 = -1.
        TopHadreco = ROOT.TLorentzVector()
        TopHadrecoTau32 = -1.
        self.SJ0isW = -1



        for ireco,reco in enumerate(recojets):
            ## Check that this jet is top tagged
            ## Top mass window was already required
            ## Check Nsubjettiness
            if reco.tau3 > 0.00001 :
                TopHadrecoTau32 = reco.tau3/reco.tau2
            if  TopHadrecoTau32 > self.maxtau32Top : continue
            ## Check that leptons are well seperated from the Fat Jet
            if reco.p4().DeltaR(lepton)  < self.mindRLepJet : continue

            if reco.subJetIdx2 >= len(recosubjets) or reco.subJetIdx1 >= len(recosubjets) :
                if self.verbose: print "Reco subjet indices not in Subjet list, Skipping"
                continue
            if reco.subJetIdx1 >= 0 and reco.subJetIdx2 >= 0 :
              
                recojetsGroomed[reco] = recosubjets[reco.subJetIdx1].p4() + recosubjets[reco.subJetIdx2].p4()
                if recosubjets[reco.subJetIdx1].p4().M() > maxrecoSJmass and recosubjets[reco.subJetIdx1].p4().M() >  recosubjets[reco.subJetIdx2].p4().M() :
                    ### Check that one of the subjets is Btagged
                    if recosubjets[reco.subJetIdx1].btagCSVV2 >  self.minBDisc  or recosubjets[reco.subJetIdx2].btagCSVV2 >  self.minBDisc :
                        if recosubjets[reco.subJetIdx1].tau1 > 0.0001 :
                            WHadrecoTau21 = recosubjets[reco.subJetIdx1].tau2 / recosubjets[reco.subJetIdx1].tau1
                            self.SJ0isW = 1
                            WHadreco = recosubjets[reco.subJetIdx1].p4()
                            TopHadreco = reco.p4()
                            break

                if recosubjets[reco.subJetIdx2].p4().M() > maxrecoSJmass and recosubjets[reco.subJetIdx1].p4().M() < recosubjets[reco.subJetIdx2].p4().M() :
                    ### Check that one of the subjets is Btagged
                    if recosubjets[reco.subJetIdx1].btagCSVV2 >  self.minBDisc  or recosubjets[reco.subJetIdx2].btagCSVV2 >  self.minBDisc :
                        if recosubjets[reco.subJetIdx2].tau1 > 0.0001 :
                            WHadrecoTau21 = recosubjets[reco.subJetIdx2].tau2 / recosubjets[reco.subJetIdx2].tau1
                            self.SJ0isW = 0
                            WHadreco = recosubjets[reco.subJetIdx2].p4()
                            TopHadreco = reco.p4()
                            break
            else :
                
                recojetsGroomed[reco] = None
                WHadreco = None
                



        if self.verbose:
            print '----'
            print 'opposite-Z recojets:'
            for recojet in recojets:
                sdmassreco = recojetsGroomed[recojet].M() if recojet in recojetsGroomed and recojetsGroomed[recojet] != None else -1.0
                print '         : %s %6.2f' % ( self.printP4( recojet),  sdmassreco )            
        if self.SJ0isW >= 0 and WHadreco != None and WHadreco.Perp() > 200. :
            self.out.fillBranch("WHadreco_pt", WHadreco.Perp())
            self.out.fillBranch("WHadreco_eta", WHadreco.Eta())
            self.out.fillBranch("WHadreco_phi", WHadreco.Phi())
            self.out.fillBranch("WHadreco_mass", WHadreco.M())
            print " W subjet tau21  %2.2f "%(WHadrecoTau21)
            self.out.fillBranch("WHadreco_tau21", WHadrecoTau21)
        #self.out.fillBranch("genmatchedAK8Subjet", self.matchedSJ)  
        #self.out.fillBranch("AK8Subjet0isMoreMassive", self.SJ0isW )
        if WHadreco.Perp() > self.WcandPtBins[0][0] and event.genmatchedAK8Subjet >= 0 :
            self.h_WcandSubjetpt.Fill(WHadreco.Perp())
            self.h_WcandSubjeteta.Fill(WHadreco.Eta())
            self.h_WcandSubjetphi.Fill(WHadreco.Phi())
            self.h_WcandSubjetmass.Fill(WHadreco.M())

            self.h_lep0pt.Fill(lepton.Perp())
            self.h_lep0eta.Fill(lepton.Eta())
            self.h_lep0phi.Fill(lepton.Phi())

            self.h_Wleppt.Fill(WcandLep.Perp())
            self.h_Wlepeta.Fill(WcandLep.Eta())
            self.h_Wlepphi.Fill(WcandLep.Phi())
            self.h_Wlepmass.Fill(WcandLep.M())

            self.h_hadToppt.Fill(TopHadreco.Perp())
            self.h_hadTopeta.Fill(TopHadreco.Eta())
            self.h_hadTopphi.Fill(TopHadreco.Phi())
            self.h_hadTopmass.Fill(TopHadreco.M())

            if event.genmatchedAK8Subjet > 0 :
                self.h_matchedAK8Subjetpt.Fill(WHadreco.Perp())
                self.h_matchedAK8Subjeteta.Fill(WHadreco.Eta())
                self.h_matchedAK8Subjetphi.Fill(WHadreco.Phi())
                self.h_matchedAK8Subjetmass.Fill(WHadreco.M())
            if event.genmatchedAK8Subjet ==0 :
                self.h_unmatchedAK8Subjetpt.Fill(WHadreco.Perp())
                self.h_unmatchedAK8Subjeteta.Fill(WHadreco.Eta())
                self.h_unmatchedAK8Subjetphi.Fill(WHadreco.Phi())
                self.h_unmatchedAK8Subjetmass.Fill(WHadreco.M())

        for ib, binhist in enumerate(self.WcandPtBins) :

            if WHadreco.Perp() > binhist[0] and WHadreco.Perp() < binhist[1] : 
                self.WcandSubjetpt[ib].Fill(WHadreco.Perp())
                self.WcandSubjeteta[ib].Fill(WHadreco.Eta())
                self.WcandSubjetphi[ib].Fill(WHadreco.Phi())
                self.WcandSubjetmass[ib].Fill(WHadreco.M())

            else : 
                self.WcandSubjetpt[ib].Fill(-1.)
                self.WcandSubjeteta[ib].Fill(-1.)
                self.WcandSubjetphi[ib].Fill(-1.)
                self.WcandSubjetmass[ib].Fill(-1.)        

        return True
# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

ttbar_semilep = lambda : TTbar_SemiLep_fullyMerged( ) 
