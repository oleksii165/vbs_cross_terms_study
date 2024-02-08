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

      _docut = 0;
      if (getOption("DOCUT") =="YES") _docut = 1;
      if (getOption("DOCUT") =="NO") _docut = 0;
      std::cout << "received docut " << _docut <<"\n";

      // The basic final-state projection:
      // all final-state particles within
      const FinalState fs;

      // photons as separate particles for final state, not dressing
      DirectFinalState photons(Cuts::abspid == PID::PHOTON && Cuts::abseta < 2.37 && Cuts::pT > 25*GeV);
      declare(photons, "photons");

      // FinalState of direct photons and bare muons and electrons in the event - ignore taus but if want to include use TauFinder
      DirectFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);
      DirectFinalState photons_for_dressing(Cuts::abspid == PID::PHOTON);
      // Dress the bare direct leptons with direct photons within dR < 0.1,
      // and apply some fiducial cuts on the dressed leptons depending on param passed
      Cut lepton_cuts = Cuts::abseta < 2.5 && Cuts::pT > 20.0*GeV; 
      DressedLeptons dressed_leps(photons_for_dressing, bare_leps, 0.1, lepton_cuts);
      declare(dressed_leps, "leptons_stable");

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering, also veto electrons(+muons but this is redundant) there
      VetoedFinalState hadrons(FinalState(Cuts::abseta < 4.4));
      hadrons.addVetoOnThisFinalState(dressed_leps);
      declare(hadrons, "hadrons");
      FastJets jetsfs(hadrons, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetsfs, "jets");

      // FS excluding the leading big pt photons, muons and neutrinos to calculate cone energy of bit pt photons
      // so basically electrons + jets left
      DirectFinalState bare_muon_for_ph_iso(Cuts::abspid == PID::MUON);
      DressedLeptons dressed_muons_for_ph_iso(photons_for_dressing, bare_muon_for_ph_iso, 0.1, lepton_cuts);
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
      book(_c["found_VBS_pair"],"found_VBS_pair");

      // Cut-flows
      _cutflows.addCutflow("sel", {"n_lep", "lep_pid_charge", "lep_pt", "m_ll", "have_photons", "have_iso_photons",
                                  "m_ll_plus_m_lly", "n_jets", "jet_pt","m_tagjets","dy_tagjets","centrality_lly",
                                  "n_gap_jets"});

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      _cutflows.fillinit();

      // Retrieve dressed leptons, sorted by pT
      Particles leptons_stable = apply<FinalState>(event, "leptons_stable").particlesByPt();
      int nlep_stable = leptons_stable.size();
      if (nlep_stable!=2)  vetoEvent; // meaning both are e,mu and not tau
      _cutflows.fillnext();

      const Particle& lep1 = leptons_stable[0];
      const Particle& lep2 = leptons_stable[1]; 
      if (lep1.pid()+lep2.pid()!=0) vetoEvent; // want opposite charge leptons of same fravour
      _cutflows.fillnext();
      if (_docut==1 && (lep1.pT() < 30.0 || lep2.pT() < 20.0)) vetoEvent; 
      _cutflows.fillnext();
      const FourMomentum fourvec_ll = lep1.mom() + lep2.mom(); 
      const double m_ll = fourvec_ll.mass()/GeV;
      if (_docut==1 && m_ll < 40.0) vetoEvent;
      _cutflows.fillnext();

      //photons
      Particles photons = apply<FinalState>(event, "photons").particlesByPt();
      if (photons.empty())  vetoEvent;
      _cutflows.fillnext();
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
        if (i_cone_to_photon_frac > 0.07)  continue;
        if (any(leptons_stable, deltaRLess(i_high_pt_photon, 0.4))) continue;
        isolated_photons += i_high_pt_photon;
        cone_to_photon_fracs.push_back(i_cone_to_photon_frac);
      }
      if (isolated_photons.empty())  vetoEvent;
      _cutflows.fillnext();

      const Particle& lead_iso_photon = isolated_photons[0];
      const FourMomentum fourvec_lly = lep1.mom() + lep2.mom() + lead_iso_photon.mom();
      const double m_lly = fourvec_lly.mass();
      if (_docut==1 && (m_ll + m_lly)<=182*GeV)  vetoEvent;
      _cutflows.fillnext();

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 20*GeV); // then will do cut on two leading pt>50
      // Remove all jets within certain dR of a dressed lepton
      idiscardIfAnyDeltaRLess(jets, leptons_stable, 0.3);
      idiscardIfAnyDeltaRLess(jets, isolated_photons, 0.4);

      int n_jets = jets.size();
      if (_docut==1 && n_jets < 2)  vetoEvent;
      _cutflows.fillnext();

      const FourMomentum tag1_jet = jets[0].mom();
      const FourMomentum tag2_jet = jets[1].mom();
      if (_docut==1 && (tag1_jet.pT()<50.0 || tag2_jet.pT()<50.0)) vetoEvent; 
      _cutflows.fillnext();

      const double m_tagjets = (tag1_jet + tag2_jet).mass()/GeV;
      if (_docut==1 && m_tagjets<500.0) vetoEvent;
      _cutflows.fillnext();

      const double dy_tagjets =  fabs(deltaRap(tag1_jet, tag2_jet));
      if (_docut==1 && dy_tagjets<1.0) vetoEvent;
      _cutflows.fillnext();

      const double centrality_lly = fabs(0.5 * (fourvec_lly.rap() - (tag1_jet.rap()+tag2_jet.rap())/2) / (tag1_jet.rap()-tag2_jet.rap()));
      if (_docut==1 && centrality_lly > 5.0)  vetoEvent;
      _cutflows.fillnext();

      int n_gap_jets = 0;
      for (int i = 0; i < n_jets; i++) {
        const double i_jet_rap = jets[i].rap();
        if ((i_jet_rap < tag1_jet.rap() && i_jet_rap > tag2_jet.rap()) || (i_jet_rap < tag2_jet.rap() && i_jet_rap > tag1_jet.rap()))  ++n_gap_jets;
      }
      if (_docut==1 && n_gap_jets > 0)  vetoEvent;
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
      _c["found_VBS_pair"]->fill();
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      double veto_survive_sumW = dbl(*_c["found_VBS_pair"]);
      double veto_survive_frac = veto_survive_sumW / sumW();
      std::cout << "survived veto, will norm m_tagjets to this: " << veto_survive_frac << "\n";
      double norm_to = veto_survive_frac*crossSection()/picobarn; // norm to generated cross-section in pb (after cuts)
      normalize(_h["m_tagjets"], norm_to);

      MSG_INFO("Cut-flow:\n" << _cutflows);
    }

    /// @}


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, Histo2DPtr> _h2;
    // map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    int _docut;
    Cutflows _cutflows;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(Zy_lly);

}
