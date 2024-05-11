// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Projections/TauFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class WmWm_lvlv : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(WmWm_lvlv);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      _docut = 0;
      if (getOption("DOCUT") =="YES") _docut = 1;
      if (getOption("DOCUT") =="NO") _docut = 0;
      std::cout << "received docut " << _docut <<"\n";

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < _jcuts["eta_tagjets"]);

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs, "jets");
      
      // FinalState of direct photons and bare muons and electrons in the event - ignore taus but if want to include use TauFinder
      DirectFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);
      DirectFinalState photons(Cuts::abspid == PID::PHOTON);
      // Dress the bare direct leptons with direct photons within dR < 0.1,
      // and apply some fiducial cuts on the dressed leptons depending on param passed
      Cut lepton_cuts;
      if (_docut==1){lepton_cuts= Cuts::abseta < _jcuts["eta_lepton"] && Cuts::pT > dbl(_jcuts["pt_lepton"])*GeV;} 
      else{lepton_cuts= Cuts::abseta < 10.0 && Cuts::pT > 0.001*GeV;}
       
      DressedLeptons dressed_leps(photons, bare_leps, 0.1, lepton_cuts);
      declare(dressed_leps, "leptons_stable");

      declare(MissingMomentum(), "METFinder");

      // Book histograms
      // plots common with others
      std::ifstream jet_hist_file(txt_dir + "/jet_hists.json");      
      json jet_hist = json::parse(jet_hist_file);
      for (json::iterator it = jet_hist.begin(); it != jet_hist.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
      }
      //lepton plots
      std::ifstream lep_hist_file(txt_dir + "/lepton_hists.json");      
      json y_hist = json::parse(lep_hist_file);
      for (json::iterator it = y_hist.begin(); it != y_hist.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
      }
      // plots that are not in other ana
      std::ifstream ana_hist_file(txt_dir + "/WmWm_lvlv_hists.json");      
      json ana_hist = json::parse(ana_hist_file);
      for (json::iterator it = ana_hist.begin(); it != ana_hist.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Retrieve dressed leptons, sorted by pT
      Particles leptons_stable = apply<FinalState>(event, "leptons_stable").particlesByPt();
      int nlep_stable = leptons_stable.size();
      if (nlep_stable!=_jcuts["n_lepton_stable"])  vetoEvent; // meaning both are e,mu and not tau

      const Particle& lep1 = leptons_stable[0];
      const Particle& lep2 = leptons_stable[1]; 
      if (lep1.charge() * lep2.charge() < 0) vetoEvent; // want same charge leptons

      const FourMomentum fourvec_ll = lep1.mom() + lep2.mom(); 
      const double m_ll = fourvec_ll.mass()/GeV;
      if (_docut==1 && m_ll<_jcuts["m_ll"]) vetoEvent; 

      const double lep1_pid = lep1.pid();
      const double lep2_pid = lep2.pid();
      double m_z = 91.18; 
      if (_docut==1){
        if (lep1_pid==PID::ELECTRON && lep2_pid==PID::ELECTRON && fabs(m_ll-m_z) < _jcuts["abs_diff_m_z"]) vetoEvent; 
      }
      
      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 30*GeV);
      Jets btagging_jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::absrap < 2.5 && Cuts::pT > 20*GeV);
      // Remove all jets within dR < 0.2 of a dressed lepton
      idiscardIfAny(jets, electrons, [](const Jet& jet, const Particle& lep) {
        return deltaR(lep, jet) < 0.4 && (lep.pT() / jet.pT() > 0.5);
      });

      idiscardIfAnyDeltaRLess(btagging_jets, leptons_stable, 0.2);



      int njets = jets.size();
      if (njets < _jcuts["n_jets"])  vetoEvent;  

      int nbtags = count(btagging_jets, hasBTag());
      if (_docut==1 && nbtags>_jcuts["n_b_jets"]) vetoEvent;

      bool foundVBSJetPair = false; // look in opposite hemispheres and pair should have highest mjj
      double max_mjj = 0;
      int tag1_jet_index = -1 ,tag2_jet_index = -1;
      for (int i = 0; i < njets; i++) {
        const Jet& i_jet = jets[i];
        for (int j = 0; j < njets; j++) {
          if (i!=j){
            const Jet& j_jet = jets[j];
            const double mjj = (i_jet.mom() + j_jet.mom()).mass()/GeV;
            const double eta_prod = i_jet.eta()*j_jet.eta();
            if  (eta_prod < 0.0 && mjj>max_mjj){
              max_mjj = mjj;
              foundVBSJetPair = true;
              tag1_jet_index = i;
              tag2_jet_index = j;
              }
            }
        }
      }
      if (tag2_jet_index < tag1_jet_index) swap(tag1_jet_index, tag2_jet_index); // organize tag jets by pt  
      if (!foundVBSJetPair)  vetoEvent;

      const FourMomentum tag1_jet = jets[tag1_jet_index].mom();
      const FourMomentum tag2_jet = jets[tag2_jet_index].mom();
      if (_docut==1 && (tag1_jet.pT()<_jcuts["pt_tagjet1"] || tag2_jet.pT()<_jcuts["pt_tagjet2"])) vetoEvent; 

      const double m_tagjets = (tag1_jet + tag2_jet).mass()/GeV;
      if (_docut==1 && m_tagjets<_jcuts["m_tagjets"]) vetoEvent;

      const double dy_tagjets = fabs(tag1_jet.rap() - tag2_jet.rap());
      if (_docut==1 && dy_tagjets<_jcuts["dy_tagjets"]) vetoEvent;

      const MissingMomentum& METfinder = apply<MissingMomentum>(event, "METFinder");
      const double scalar_MET = METfinder.missingPt()/GeV;
      if (_docut==1 && scalar_MET<_jcuts["pt_MET"]) vetoEvent;

      const FourMomentum fourvec_MET = METfinder.missingMomentum();
      const double m_T = (fourvec_MET + fourvec_ll).mass()/GeV;
    
      //jet plots
      _h["n_jet"]->fill(njets);
      _h["n_bjet"]->fill(nbtags);
      _h["pt_tagjet1"]->fill(tag1_jet.pt());
      _h["pt_tagjet2"]->fill(tag2_jet.pt());
      // above this worked before
      _h["m_tagjets"]->fill(m_tagjets);
      _h["dy_tagjets"]->fill(dy_tagjets);
      _h2["m_dy_tagjets"]->fill(m_tagjets,dy_tagjets);
      _h["eta_tagjet1"]->fill(tag1_jet.eta());
      _h["eta_tagjet2"]->fill(tag2_jet.eta());
      _h["eta_tagjets"]->fill(tag1_jet.eta()); _h["eta_tagjets"]->fill(tag2_jet.eta()); // fill both to the same hists
      _h["deta_tagjets"]->fill(deltaEta(tag1_jet,tag2_jet));
      _h["phi_tagjet1"]->fill(tag1_jet.phi());
      _h["phi_tagjet2"]->fill(tag2_jet.phi());
      _h["dphi_tagjets"]->fill(deltaPhi(tag1_jet,tag2_jet));
      //lepton plots
      _h2["leptons_pids"]->fill(lep1_pid,lep2_pid); // to check that only have e and mu and of same charge
      _h["n_lepton_stable"]->fill(nlep_stable);
      _h["lepton_pt"]->fill(lep1.pT()); _h["lepton_pt"]->fill(lep2.pT()); // fill both to the same hists
      _h["lepton_eta"]->fill(lep1.eta());
      _h["lepton_eta"]->fill(lep2.eta());
      _h["m_ll"]->fill(m_ll);
      if (lep1_pid == PID::ELECTRON && lep2_pid == PID::ELECTRON){_h["m_ee"]->fill(m_ll);}
      else if (lep1_pid == PID::MUON && lep2_pid == PID::MUON){_h["m_mumu"]->fill(m_ll);}
      else{_h["m_emu"]->fill(m_ll);}
      // other
      _h["MET"]->fill(scalar_MET);
      _h["m_T"]->fill(m_T);
      _c["found_VBS_pair"]->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double veto_survive_sumW = dbl(*_c["found_VBS_pair"]);
      double veto_survive_frac = veto_survive_sumW / sumW();
      std::cout << "survived veto, will norm to this: " << veto_survive_frac << "\n";
      double norm_to = veto_survive_frac*crossSection()/picobarn; // norm to generated cross-section in pb (after cuts)
      
      std::vector<std::string> hist_names_1d = {"n_jet","n_bjet","pt_tagjet1","pt_tagjet2","m_tagjets",
      "dy_tagjets","eta_tagjet1","eta_tagjet2","eta_tagjets", "deta_tagjets", 
      "phi_tagjet1","phi_tagjet2","dphi_tagjets",
      "n_lepton_stable","lepton_pt","lepton_eta",
      "m_ll","m_ee","m_mumu","m_emu","MET","m_T"};       
      for (auto&& i_name : hist_names_1d){ normalize(_h[i_name], norm_to);}
      // also norm few 2d
      normalize(_h2["m_dy_tagjets"], norm_to);
      normalize(_h2["leptons_pids"],norm_to);

    }

    /// @}


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, Histo2DPtr> _h2;
    // map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    int _docut;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(WmWm_lvlv);

}