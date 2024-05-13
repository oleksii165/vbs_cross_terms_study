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
#include "Rivet/Projections/InvisibleFinalState.hh"
#include "Rivet/Tools/Cutflow.hh"
#include "Rivet/Tools/RivetHepMC.hh"
#include <fstream>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace Rivet {


  /// @brief Add a short analysis description here
  class WZ_lllv : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(WZ_lllv);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      std::string out_dir = getOption("OUTDIR");
      
      _docut = 0; // most cuts on number of particles are always applied to avoid segfault
      if (out_dir.find("DOCUT_YES") != string::npos) _docut = 1;
      // if (out_dir.find("DOCUT_NO") != string::npos) _docut = 0;
      std::cout << "++++++received outidir" << out_dir << "meaning _docut is " << _docut << "\n";

      std::string jsonfilestr =  "WZ_lllv_cuts.json";
      std::cout << "++++++assume .json for this WZ_lllv" << "is " << jsonfilestr << "\n";
      std::ifstream json_file(jsonfilestr);
      
      _jcuts = json::parse(json_file);
      std::cout << "++++++ to check 1 var got photon pt min" << _jcuts["m_tagjets"] << "\n";
      _electron_eta_cut = (Cuts::absetaIn(_jcuts["eta_lepton_electron"][0][0], _jcuts["eta_lepton_electron"][0][1])) || 
                                (Cuts::absetaIn(_jcuts["eta_lepton_electron"][1][0], _jcuts["eta_lepton_electron"][1][1]));
      _muon_eta_cut = Cuts::absetaIn(0.0, _jcuts["eta_lepton_muon"]);
      _lepton_stage1_pt_cut = Cuts::pT > dbl(_jcuts["pt_lepton"])*GeV;      

      // The basic final-state projection:
      // all final-state particles within
      const FinalState fs;

      // FinalState of direct photons and bare muons and electrons in the event - ignore taus but if want to include use TauFinder
      DirectFinalState bare_e(Cuts::abspid == PID::ELECTRON);
      DirectFinalState bare_mu(Cuts::abspid == PID::MUON);
      DirectFinalState photons_for_dressing(Cuts::abspid == PID::PHOTON);
      // Dress the bare direct leptons with direct photons within dR < 0.1,
      // and apply some fiducial cuts on the dressed leptons depending on param passed
      DressedLeptons dressed_e(photons_for_dressing, bare_e, 0.1);
      DressedLeptons dressed_mu(photons_for_dressing, bare_mu, 0.1);
      // declare(dressed_leps, "leptons_stable");
      declare(dressed_e, "e_stable");
      declare(dressed_mu, "mu_stable");

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering, also veto electrons(+muons but this is redundant) there
      VetoedFinalState hadrons(FinalState(Cuts::absetaIn(0.0, _jcuts["eta_tagjets"])));
      hadrons.addVetoOnThisFinalState(dressed_e);
      hadrons.addVetoOnThisFinalState(dressed_mu);
      declare(hadrons, "hadrons");
      FastJets jetsfs(hadrons, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetsfs, "jets");

      declare(MissingMomentum(), "METFinder");

      // Book things and save names to normalize later
      
      // plots common with others
      std::ifstream jet_hist_file("jet_hists.json");
      json jet_hist = json::parse(jet_hist_file);
      for (json::iterator it = jet_hist.begin(); it != jet_hist.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
      }
      std::ifstream lep_hist_file("lepton_hists.json");
      json lep_hist = json::parse(lep_hist_file);
      for (json::iterator it = lep_hist.begin(); it != lep_hist.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
      }
      // plots that are not in other ana
      std::ifstream ana_hist_file("WZ_lllv_hists.json");
      json ana_hist = json::parse(ana_hist_file);
      for (json::iterator it = ana_hist.begin(); it != ana_hist.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
      }
      
      //counter for efficiency
      book(_c["pos_w_initial"],"pos_w_initial");
      book(_c["pos_w_final"],"pos_w_final");
      book(_c["neg_w_initial"],"neg_w_initial");
      book(_c["neg_w_final"],"neg_w_final");
      
      // Cut-flows
      _cutflows.addCutflow("WZ_lllv_selections", {"three_leptons", "one_SFOC_or_more", "mZ_mismatch",
                                                  "pt_W_lepton", "m_W_T", "dR__leptons", "n_jets","n_b_jets", 
                                                  "pt_tagjets", "m_tagjets", "tagj_opposite_hemisph"});

    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // save weights before cuts
      double ev_nominal_weight =  event.weights()[0];
      if (ev_nominal_weight>=0){_c["pos_w_initial"]->fill();} // dont need anything in bracket as this will be weight on weight
      else {_c["neg_w_initial"]->fill();}

      _cutflows.fillinit();

      // Retrieve dressed leptons, sorted by pT
      Particles e_stable;
      Particles mu_stable;
      if (_docut==1){
        e_stable = apply<FinalState>(event, "e_stable").particlesByPt(_electron_eta_cut && _lepton_stage1_pt_cut);
        mu_stable = apply<FinalState>(event, "mu_stable").particlesByPt(_muon_eta_cut && _lepton_stage1_pt_cut);
      }
      else{
        e_stable = apply<FinalState>(event, "e_stable").particlesByPt();
        mu_stable = apply<FinalState>(event, "mu_stable").particlesByPt();
      }
      // will ahve either e or mu pair so don't need to sort combined e+mu array
      Particles leptons = e_stable + mu_stable; 

      int nlep = leptons.size();
      if (nlep!=_jcuts["n_lepton_stable"])  vetoEvent; // meaning both are e,mu and not tau
      _cutflows.fillnext();


      // check SFOC and m_ll and save possible pairs
      std::vector<std::vector<int>> pairs_ind;
      for (int i = 0; i < nlep; i++) {
          const Particle& i_lep = leptons[i];
          for (int j = 0; j < nlep; j++) {
              // to avoid comparing to itself and double couting like (i,j),(j,i) as two separate pairs
              if (i==j || j>i) continue;
              const Particle& j_lep = leptons[j];
              int i_sum_pids = i_lep.pid() + j_lep.pid(); // to have SFOC will be 0
              double i_m_ll = (i_lep.mom() + j_lep.mom()).mass()/GeV;
              if (i_sum_pids==0 && i_m_ll>_jcuts["m_ll_all_pairs"]) pairs_ind.push_back({i,j});
          }
      }
      if (pairs_ind.size()<1) vetoEvent;
      _cutflows.fillnext();

      // order found pairs by how close they are m_z
      std::sort(pairs_ind.begin(), pairs_ind.end(), [this,leptons](std::vector<int> &ind_pair_1, std::vector<int> &ind_pair_2) {
          double dist_pair_1 = pair_m_dist_m_z(leptons[ind_pair_1[0]], leptons[ind_pair_1[1]]);
          double dist_pair_2 = pair_m_dist_m_z(leptons[ind_pair_2[0]], leptons[ind_pair_2[1]]);
          return dist_pair_1 < dist_pair_2; // closest to m_z will be first
      });

      const FourMomentum& Z_lep_1 = leptons[pairs_ind[0][0]].mom();
      const FourMomentum& Z_lep_2 = leptons[pairs_ind[0][1]].mom();
      if (_docut==1 && pair_m_dist_m_z(Z_lep_1, Z_lep_2)>_jcuts["abs_diff_m_z"]) vetoEvent;
      _cutflows.fillnext();

      int ind_W_lep = 500;
      for (int i = 0; i < nlep; i++) {
              if (i==pairs_ind[0][0] || i==pairs_ind[0][1]) continue;
              ind_W_lep = i;
        }
      if (ind_W_lep==500) vetoEvent;

      const FourMomentum& W_lep = leptons[ind_W_lep].mom();
      if (_docut==1 && W_lep.pT() < _jcuts["pt_W_lepton"]) vetoEvent;
      _cutflows.fillnext();

      const MissingMomentum& METfinder = apply<MissingMomentum>(event, "METFinder");
      const double pt_MET = METfinder.missingPt()/GeV;
      const FourMomentum fourvec_MET = METfinder.missingMomentum();
      const double d_phi_MET_lep = deltaPhi(fourvec_MET.phi(), W_lep.phi());
      const double m_W_T = sqrt( 2*W_lep.pT()*pt_MET*(1 - cos(d_phi_MET_lep)) );
      if (_docut==1 && m_W_T<_jcuts["m_W_T"]) vetoEvent;
      _cutflows.fillnext();

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 20*GeV); 
      Jets btagging_jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::absrap < 2.5 && Cuts::pT > 20*GeV);
      // Remove all jets within certain dR of a dressed lepton
      idiscardIfAnyDeltaRLess(jets, leptons, _jcuts["dR_lepton_jet"]);
      idiscardIfAnyDeltaRLess(btagging_jets, leptons, _jcuts["dR_lepton_jet"]);
      

      double dR_Z_lep_1_lep_2 = deltaR(Z_lep_1, Z_lep_2);  
      double dR_Z_lep_1_lep_W = deltaR(Z_lep_1, W_lep);  
      double dR_Z_lep_2_lep_W = deltaR(Z_lep_2, W_lep);  
      if (_docut==1 && dR_Z_lep_1_lep_2 < _jcuts["dR_lepton1_lepton2_Z"]) vetoEvent; 
      if (_docut==1 && dR_Z_lep_1_lep_W < _jcuts["dR_Z_leptons_W_lepton"]) vetoEvent; 
      if (_docut==1 && dR_Z_lep_2_lep_W < _jcuts["dR_Z_leptons_W_lepton"]) vetoEvent; 
      _cutflows.fillnext();

      int n_jets = jets.size();
      if (n_jets < _jcuts["n_jets"]) vetoEvent;
      _cutflows.fillnext();

      int n_b_jets = count(btagging_jets, hasBTag());
      if (_docut==1 && n_b_jets>_jcuts["n_b_jets"]) vetoEvent;
      _cutflows.fillnext();

      const FourMomentum tag1_jet = jets[0].mom();
      const FourMomentum tag2_jet = jets[1].mom();
      if (_docut==1 && (tag1_jet.pT()<dbl(_jcuts["pt_tagjet1"])*GeV || tag2_jet.pT()<dbl(_jcuts["pt_tagjet2"])*GeV)) vetoEvent; 
      _cutflows.fillnext();

      const double m_tagjets = (tag1_jet + tag2_jet).mass()/GeV;
      if (_docut==1 && m_tagjets<dbl(_jcuts["m_tagjets"])*GeV) vetoEvent;
      _cutflows.fillnext();

      const double dy_tagjets =  fabs(deltaRap(tag1_jet, tag2_jet));
      const double prod_y_tagjets =  tag1_jet.rap()*tag2_jet.rap();
      if (_docut==1 && prod_y_tagjets > 0) vetoEvent;
      _cutflows.fillnext();

      const double m_WZ_T_term1 =  pow(Z_lep_1.pT() + Z_lep_2.pT() + W_lep.pT() + pt_MET, 2);
      const double m_WZ_T_term2 =  pow(Z_lep_1.px() + Z_lep_2.px() + W_lep.px() + fourvec_MET.px(), 2);
      const double m_WZ_T_term3 =  pow(Z_lep_1.py() + Z_lep_2.py() + W_lep.py() + fourvec_MET.py(), 2);
      const double m_WZ_T = sqrt(m_WZ_T_term1 - m_WZ_T_term2 - m_WZ_T_term3);
      
      //jet common 
      _h["n_jets"]->fill(n_jets);
      _h["pt_tagjet1"]->fill(tag1_jet.pt());
      _h["pt_tagjet2"]->fill(tag2_jet.pt());
      _h["eta_tagjets"]->fill(tag1_jet.eta()); _h["eta_tagjets"]->fill(tag2_jet.eta());
      _h["phi_tagjets"]->fill(tag1_jet.phi()); _h["phi_tagjets"]->fill(tag2_jet.phi());
      _h["m_tagjets"]->fill(m_tagjets);
      _h["dy_tagjets"]->fill(dy_tagjets);
      _h["dphi_tagjets"]->fill(deltaPhi(tag1_jet,tag2_jet));
      //lepton plots
      _h["n_lepton_stable"]->fill(nlep);
      //
      _h["pt_lepton"]->fill(Z_lep_1.pT());
      _h["pt_lepton"]->fill(Z_lep_2.pT()); 
      _h["pt_lepton"]->fill(W_lep.pT()); 
      //
      _h["eta_lepton"]->fill(Z_lep_1.eta()); 
      _h["eta_lepton"]->fill(Z_lep_2.eta());
      _h["eta_lepton"]->fill(W_lep.eta());
      // analysis-specific
      _h["pt_MET"]->fill(pt_MET);
      _h["m_W_T"]->fill(m_W_T);
      _h["m_WZ_T"]->fill(m_WZ_T);
      _h["n_b_jets"]->fill(n_b_jets);
      //
      _h["dR_leptons"]->fill(dR_Z_lep_1_lep_2);
      _h["dR_leptons"]->fill(dR_Z_lep_1_lep_W);
      _h["dR_leptons"]->fill(dR_Z_lep_2_lep_W);

      // save weights after cuts
      if (ev_nominal_weight>=0){_c["pos_w_final"]->fill();}
      else {_c["neg_w_final"]->fill();}

    }

    double pair_m_dist_m_z(const FourMomentum& lep_1, const FourMomentum& lep_2){
        double i_m_ll = (lep_1 + lep_2).mass()/GeV;
        double i_m_ll_dist_m_z = fabs(i_m_ll-91.18);
        return i_m_ll_dist_m_z;
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      std::cout << _cutflows.str();

      double pos_w_sum_initial = dbl(*_c["pos_w_initial"]); // from which also number of entries can be obtained
      double neg_w_sum_initial = dbl(*_c["neg_w_initial"]);
      double pos_w_sum_final = dbl(*_c["pos_w_final"]);
      double neg_w_sum_final = dbl(*_c["neg_w_final"]);
      MSG_INFO("\n pos weights initial final ratio " << pos_w_sum_initial <<" " << pos_w_sum_final <<" "<< pos_w_sum_final/pos_w_sum_initial << "\n" );
      MSG_INFO("\n neg weights initial final ratio " << neg_w_sum_initial <<" " << neg_w_sum_final <<" "<< neg_w_sum_final/neg_w_sum_initial << "\n" );

      // normalize all to 1 since in case of mostly negative weights not clear what it will do
      for (auto & i_name : _hist_names){ 
        std::cout << "normalizeing hist " << i_name <<" to 1; " ;
        normalize(_h[i_name], 1.0);
      }
    }

    /// @}

    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, Histo2DPtr> _h2;
    map<string, CounterPtr> _c;
    int _docut;
    Cut _electron_eta_cut;
    Cut _muon_eta_cut;
    Cut _lepton_stage1_pt_cut;
    json _jcuts;
    Cutflows _cutflows;
    std::vector<std::string> _hist_names;

    /// @}

  };


  RIVET_DECLARE_PLUGIN(WZ_lllv);

}
