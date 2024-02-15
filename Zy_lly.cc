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
  class Zy_lly : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(Zy_lly);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      
      std::string out_dir = getOption("OUTDIR");
      
      _docut = 0;
      if (out_dir.find("DOCUT_YES") != string::npos) _docut = 1;
      if (out_dir.find("DOCUT_NO") != string::npos) _docut = 0;
      std::cout << "++++++received outidir" << out_dir << "meaning _docut is " << _docut << "\n";

      int stop_prod_dec_dir_ind = out_dir.find("/user.okurdysh.MadGraph");
      std::string almost_prod_dec = out_dir.substr(0, stop_prod_dec_dir_ind);
      int st_prod_dec_dir_ind = almost_prod_dec.find("/eft_files/");
      std::string prod_dec = almost_prod_dec.substr(st_prod_dec_dir_ind + strlen("/eft_files/"));
      std::string jsonfilestr =  "/exp/atlas/kurdysh/vbs_cross_terms_study/plotting/" + prod_dec + "_cuts.json"; 
      std::cout << "++++++assume .json for this prod_dec" << prod_dec << "is " << jsonfilestr << "\n";
      std::ifstream json_file(jsonfilestr);
      
      _jcuts = json::parse(json_file);
      std::cout << "++++++ to check 1 var got photon pt min" << _jcuts["photon_pt"] << "\n";
      _electron_eta_cut = (Cuts::absetaIn(_jcuts["e_eta"][0][0], _jcuts["e_eta"][0][1])) || 
                                (Cuts::absetaIn(_jcuts["e_eta"][1][0], _jcuts["e_eta"][1][1]));
      _photon_eta_cut = (Cuts::absetaIn(_jcuts["photon_eta"][0][0], _jcuts["photon_eta"][0][1])) || 
                                (Cuts::absetaIn(_jcuts["photon_eta"][1][0], _jcuts["photon_eta"][1][1]));
      _muon_eta_cut = Cuts::absetaIn(0.0, _jcuts["muon_eta"]);
      _lepton2_pt_cut = Cuts::pT > dbl(_jcuts["lepton2_pt"])*GeV;      

      // The basic final-state projection:
      // all final-state particles within
      const FinalState fs;

      // photons as separate particles for final state, not dressing
      DirectFinalState photons(Cuts::abspid == PID::PHOTON);
      declare(photons, "photons");

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
      VetoedFinalState hadrons(FinalState(Cuts::absetaIn(0.0, _jcuts["jet_eta"])));
      hadrons.addVetoOnThisFinalState(dressed_e);
      hadrons.addVetoOnThisFinalState(dressed_mu);
      declare(hadrons, "hadrons");
      FastJets jetsfs(hadrons, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetsfs, "jets");

      // FS excluding the leading big pt photons, muons and neutrinos to calculate cone energy of bit pt photons
      // so basically electrons + jets left
      DirectFinalState bare_muon_for_ph_iso(Cuts::abspid == PID::MUON);
      DressedLeptons dressed_muons_for_ph_iso(photons_for_dressing, bare_muon_for_ph_iso, 0.1, _muon_eta_cut &&_lepton2_pt_cut);
      declare(dressed_muons_for_ph_iso, "muons_for_ph_iso");
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(photons);
      vfs.addVetoOnThisFinalState(dressed_muons_for_ph_iso);
      vfs.addVetoOnThisFinalState(InvisibleFinalState());
      declare(vfs, "ElAndJetsForPhotonIsoCalc");

      // Book histograms
      int n_nbins = 10;
      int n_pt = 200;
      double max_pt = 3000.0;
      int n_rap = 60;
      double max_rap = 6.0;
      int n_pid = 15;
      //jet plots
      book(_h["n_jet"], "n_jet", n_nbins, 0.0, n_nbins);
      book(_h["n_gap_jet"], "n_gap_jet", n_nbins, 0.0, n_nbins);
      book(_h["pt_tagjet1"], "pt_tagjet1", n_pt, 0.0, max_pt); 
      book(_h["pt_tagjet2"], "pt_tagjet2", n_pt, 0.0, max_pt);
      book(_h["m_tagjets"], "m_tagjets", n_pt, 0.0, max_pt);
      book(_h["dy_tagjets"], "dy_tagjets", n_rap, 0, max_rap);
      book(_h["eta_tagjets"], "eta_tagjets", n_rap, -1*max_rap, max_rap);
      book(_h["dphi_tagjets"], "dphi_tagjets", n_rap, 0, max_rap);
      // //lepton plots
      book(_h2["leptons_pids"], "leptons_pids", 2*n_pid, -1*n_pid, n_pid, 2*n_pid, -1*n_pid, n_pid);
      book(_h["n_lepton_stable"], "n_lepton_stable", n_nbins, 0.0, n_nbins);
      book(_h["lepton_pt"], "lepton_pt", int(n_pt), 0.0, max_pt);
      book(_h["lepton_eta"], "lepton_eta", n_rap, -1*max_rap, max_rap);
      book(_h["m_ll"], "m_ll", int(n_pt/10), 0.0, max_pt/10);
      //photon plots
      book(_h["n_photons_iso"], "n_photons_iso", n_nbins, 0.0, n_nbins);
      book(_h["frac_photons_iso"], "frac_photons_iso", 40, -2.0, 2.0); // relative to all photons
      book(_h["cone_to_photon_frac"], "cone_to_photon_frac", 20, 0.0, 0.2);
      book(_h["photon_iso_eta"], "photon_iso_eta", n_rap, -1*max_rap, max_rap);
      book(_h["photon_iso_pt"], "photon_iso_pt", int(n_pt), 0.0, max_pt);
      // lly
      book(_h["m_lly"], "m_lly", int(n_pt), 0.0, max_pt);
      book(_h["centrality_lly"], "centrality_lly", int(n_nbins), 0.0, n_nbins);
      //other
      book(_c["pos_w_initial"],"pos_w_initial");
      book(_c["pos_w_final"],"pos_w_final");
      book(_c["neg_w_initial"],"neg_w_initial");
      book(_c["neg_w_final"],"neg_w_final");
      // Cut-flows
      _cutflows.addCutflow("Zy_lly_selections", {"n_lep_ok_pt2_eta", "lep_pid_charge", "lep_pt1", "m_ll", "have_iso_photons_ok_pt_eta",
                                  "m_ll_plus_m_lly", "n_jets", "jet_pt","m_tagjets","dy_tagjets","centrality_lly",
                                  "n_gap_jets"});
      // setup for  file used for drawing images
      if (_docut==1){
        std::ofstream pic_csv (out_dir + "/Rivet.csv", std::ofstream::out);
        pic_csv << "ev_num;tagjet1_eta;tagjet2_eta;lep1_eta;lep2_eta;gamma_eta;tagjet1_pt;tagjet2_pt;lep1_pt;lep2_pt;gamma_pt \n";
        pic_csv.close();
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      double ev_nominal_weight =  event.weights()[0];
      // MSG_INFO("filling weight" << ev_nominal_weight << "\n");
      if (ev_nominal_weight>=0){_c["pos_w_initial"]->fill();} // dont need anything in bracket as this will be weight on weight
      else {_c["neg_w_initial"]->fill();}

      _cutflows.fillinit();

      // Retrieve dressed leptons, sorted by pT
      Particles e_stable;
      Particles mu_stable;
      if (_docut==1){
        e_stable = apply<FinalState>(event, "e_stable").particlesByPt(_electron_eta_cut && _lepton2_pt_cut);
        mu_stable = apply<FinalState>(event, "mu_stable").particlesByPt(_muon_eta_cut && _lepton2_pt_cut);
      }
      else{
        e_stable = apply<FinalState>(event, "e_stable").particlesByPt();
        mu_stable = apply<FinalState>(event, "mu_stable").particlesByPt();
      }
      // will ahve either e or mu pair so don't need to sort combined e+mu array
      Particles leptons_stable = e_stable + mu_stable; 

      int nlep_stable = leptons_stable.size();
      if (nlep_stable!=_jcuts["n_lep"])  vetoEvent; // meaning both are e,mu and not tau
      _cutflows.fillnext();

      const Particle& lep1 = leptons_stable[0];
      const Particle& lep2 = leptons_stable[1];
      if (lep1.pid()+lep2.pid()!=0) vetoEvent; // want opposite charge leptons of same fravour
      _cutflows.fillnext();
      if (_docut==1 && (lep1.pT() < dbl(_jcuts["lepton1_pt"])*GeV)) vetoEvent;
      _cutflows.fillnext();
      const FourMomentum fourvec_ll = lep1.mom() + lep2.mom(); 
      const double m_ll = fourvec_ll.mass()/GeV;
      if (_docut==1 && m_ll < dbl(_jcuts["m_ll"])*GeV) vetoEvent;
      _cutflows.fillnext();

      //photons
      Particles photons = apply<FinalState>(event, "photons").particlesByPt(_photon_eta_cut);
      if (photons.empty())  vetoEvent;
      //photon cone calculation and photon OR with leptons 
      Particles isolated_photons;
      std::vector<double> cone_to_photon_fracs = {};
      Particles cone_sum_particles = apply<VetoedFinalState>(event, "ElAndJetsForPhotonIsoCalc").particles();
      for (const Particle &i_high_pt_photon : photons) {
        // check photon isolation
        double i_high_pt_photon_cone_E = 0.0;
        for (const Particle &i_p_cone : cone_sum_particles) {
          if (deltaR(i_high_pt_photon, i_p_cone) < 0.2) { // etcone20
            i_high_pt_photon_cone_E += i_p_cone.Et();
          }
        }
        double i_cone_to_photon_frac = i_high_pt_photon_cone_E / i_high_pt_photon.pT(); 
        if (i_cone_to_photon_frac > _jcuts["cone_to_photon_frac"])  continue;
        if (any(leptons_stable, deltaRLess(i_high_pt_photon, 0.4))) continue;
        isolated_photons += i_high_pt_photon;
        cone_to_photon_fracs.push_back(i_cone_to_photon_frac);
      }
      if (isolated_photons.empty())  vetoEvent;
      const Particle& lead_iso_photon = isolated_photons[0];
      if (_docut==1 && lead_iso_photon.pT() < dbl(_jcuts["photon_pt"])*GeV) vetoEvent;
      _cutflows.fillnext();
      
      const FourMomentum fourvec_lly = lep1.mom() + lep2.mom() + lead_iso_photon.mom();
      const double m_lly = fourvec_lly.mass();
      if (_docut==1 && (m_ll + m_lly)<=dbl(_jcuts["m_lly"])*GeV)  vetoEvent;
      _cutflows.fillnext();

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 20*GeV); // then will do cut on two leading pt>50
      // Remove all jets within certain dR of a dressed lepton
      idiscardIfAnyDeltaRLess(jets, leptons_stable, 0.3);
      idiscardIfAnyDeltaRLess(jets, isolated_photons, 0.4);

      int n_jets = jets.size();
      if (_docut==1 && n_jets < _jcuts["n_jets"]) vetoEvent;
      _cutflows.fillnext();

      const FourMomentum tag1_jet = jets[0].mom();
      const FourMomentum tag2_jet = jets[1].mom();
      if (_docut==1 && (tag1_jet.pT()<dbl(_jcuts["jet_pt"])*GeV || tag2_jet.pT()<dbl(_jcuts["jet_pt"])*GeV)) vetoEvent; 
      _cutflows.fillnext();

      const double m_tagjets = (tag1_jet + tag2_jet).mass()/GeV;
      if (_docut==1 && m_tagjets<dbl(_jcuts["m_tagjets"])*GeV) vetoEvent;
      _cutflows.fillnext();

      const double dy_tagjets =  fabs(deltaRap(tag1_jet, tag2_jet));
      if (_docut==1 && dy_tagjets<_jcuts["dy_tagjets"]) vetoEvent;
      _cutflows.fillnext();

      const double centrality_lly = fabs(0.5 * (fourvec_lly.rap() - (tag1_jet.rap()+tag2_jet.rap())/2) / (tag1_jet.rap()-tag2_jet.rap()));
      if (_docut==1 && centrality_lly > _jcuts["centrality_lly"])  vetoEvent;
      _cutflows.fillnext();

      int n_gap_jets = 0;
      for (int i = 0; i < n_jets; i++) {
        const double i_jet_rap = jets[i].rap();
        if ((i_jet_rap < tag1_jet.rap() && i_jet_rap > tag2_jet.rap()) || (i_jet_rap < tag2_jet.rap() && i_jet_rap > tag1_jet.rap()))  ++n_gap_jets;
      }
      if (_docut==1 && n_gap_jets > _jcuts["n_gap_jets"])  vetoEvent;
      _cutflows.fillnext();

      //jet plots
      _h["n_jet"]->fill(n_jets);
      _h["n_gap_jet"]->fill(n_gap_jets);
      _h["pt_tagjet1"]->fill(tag1_jet.pt());
      _h["pt_tagjet2"]->fill(tag2_jet.pt());
      _h["m_tagjets"]->fill(m_tagjets);
      _h["dy_tagjets"]->fill(dy_tagjets);
      _h["eta_tagjets"]->fill(tag1_jet.eta()); _h["eta_tagjets"]->fill(tag2_jet.eta());
      _h["dphi_tagjets"]->fill(deltaPhi(tag1_jet,tag2_jet));
      //lepton plots
      _h2["leptons_pids"]->fill(lep1.pid(),lep2.pid());
      _h["n_lepton_stable"]->fill(nlep_stable);
      _h["lepton_pt"]->fill(lep1.pT()); _h["lepton_pt"]->fill(lep2.pT()); 
      _h["lepton_eta"]->fill(lep1.eta()); _h["lepton_eta"]->fill(lep2.eta());
      _h["m_ll"]->fill(m_ll);
      //photon plots
      _h["n_photons_iso"]->fill(isolated_photons.size());
      _h["frac_photons_iso"]->fill(photons.size()/isolated_photons.size());
      for (long unsigned int i = 0; i < cone_to_photon_fracs.size(); i++) {_h["cone_to_photon_frac"]->fill(cone_to_photon_fracs[i]);}
      for (long unsigned int i = 0; i < isolated_photons.size(); i++) {
        _h["photon_iso_eta"]->fill(isolated_photons[i].eta());
        _h["photon_iso_pt"]->fill(isolated_photons[i].pT());
        }
      // lly
      _h["m_lly"]->fill(m_lly);
      _h["centrality_lly"]->fill(centrality_lly);
      // other
      if (ev_nominal_weight>=0){_c["pos_w_final"]->fill();}
      else {_c["neg_w_final"]->fill();}

       // file used for drawing images
      if (_docut==1){
        int ev_num =  _c["pos_w_initial"]->numEntries() + _c["neg_w_initial"]->numEntries();
        // later to get average only makes sense todo it on +,- jets separately
        // similarly for leps where one will have slighyl bigger eta
        int ind_bigger_eta_tagjet = (tag1_jet.eta() >  tag2_jet.eta()) ? 0 : 1;
        int ind_bigger_eta_lep = (lep1.eta() >  lep2.eta()) ? 0 : 1;
        // find other index through ugly negating via boolean
        int ind_smaller_eta_tagjet = static_cast<int>(!static_cast<bool>(ind_bigger_eta_tagjet));
        int ind_smaller_eta_lep = static_cast<int>(!static_cast<bool>(ind_bigger_eta_lep));
        // pulling file into common with init() _fout didn't work so re-open
        std::ofstream pic_csv (getOption("OUTDIR") + "/Rivet.csv", std::ofstream::app); 
        pic_csv << ev_num << ";";
        //etas
        pic_csv << jets[ind_bigger_eta_tagjet].eta() << ";";
        pic_csv << jets[ind_smaller_eta_tagjet].eta() << ";";
        pic_csv << leptons_stable[ind_bigger_eta_lep].eta() << ";";
        pic_csv << leptons_stable[ind_smaller_eta_lep].eta() << ";";
        pic_csv << lead_iso_photon.eta() << ";";
        //pts
        pic_csv << jets[ind_bigger_eta_tagjet].pt() << ";";
        pic_csv << jets[ind_smaller_eta_tagjet].pt() << ";";
        pic_csv << leptons_stable[ind_bigger_eta_lep].pt() << ";";
        pic_csv << leptons_stable[ind_smaller_eta_lep].pt() << ";";
        pic_csv << lead_iso_photon.pt() << "\n";
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      std::string cut_str = _cutflows.str();
      std::string cutflow_file = getOption("OUTDIR") + "/cutflow.txt";
      std::ofstream ofs (cutflow_file, std::ofstream::out); // same for all variations so can overwrite and last will be correct
      ofs << cut_str;
      ofs.close();

      double pos_w_sum_initial = dbl(*_c["pos_w_initial"]); // from which also number of entries can be obtained
      double neg_w_sum_initial = dbl(*_c["neg_w_initial"]);
      double pos_w_sum_final = dbl(*_c["pos_w_final"]);
      double neg_w_sum_final = dbl(*_c["neg_w_final"]);
      MSG_INFO("\n pos weights initial final ratio " << pos_w_sum_initial <<" " << pos_w_sum_final <<" "<< pos_w_sum_final/pos_w_sum_initial << "\n" );
      MSG_INFO("\n neg weights initial final ratio " << neg_w_sum_initial <<" " << neg_w_sum_final <<" "<< neg_w_sum_final/neg_w_sum_initial << "\n" );
    
    }

    /// @}

    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, Histo2DPtr> _h2;
    map<string, CounterPtr> _c;
    int _docut;
    Cut _electron_eta_cut;
    Cut _photon_eta_cut;
    Cut _muon_eta_cut;
    Cut _lepton2_pt_cut;
    json _jcuts;
    std::string _jsonfilestr;
    Cutflows _cutflows;

    /// @}


  };


  RIVET_DECLARE_PLUGIN(Zy_lly);

}
