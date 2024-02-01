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

namespace Rivet {


  /// @brief Add a short analysis description here
  class Wmy_lvy : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(Wmy_lvy);


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
      DirectFinalState photons(Cuts::abspid == PID::PHOTON && Cuts::abseta < 2.37 && Cuts::pT > 22*GeV);
      declare(photons, "photons");

      // FinalState of direct photons and bare muons and electrons in the event - ignore taus but if want to include use TauFinder
      DirectFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);
      DirectFinalState photons_for_dressing(Cuts::abspid == PID::PHOTON);
      // Dress the bare direct leptons with direct photons within dR < 0.1,
      // and apply some fiducial cuts on the dressed leptons depending on param passed
      Cut lepton_cuts = Cuts::abseta < 2.47 && Cuts::pT > 30.0*GeV; 
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

      declare(MissingMomentum(), "METFinder");

      // Book histograms
      int n_nbins = 10;
      int n_pt = 200;
      double max_pt = 3000.0;
      int n_rap = 60;
      double max_rap = 6.0;
      //jet plots
      book(_h["n_jet"], "n_jet", n_nbins, 0.0, n_nbins); // @TODO think which hists want 
      book(_h["n_gap_jet"], "n_gap_jet", n_nbins, 0.0, n_nbins);
      book(_h["pt_tagjet1"], "pt_tagjet1", n_pt, 0.0, max_pt); 
      book(_h["pt_tagjet2"], "pt_tagjet2", n_pt, 0.0, max_pt);
      book(_h["m_tagjets"], "m_tagjets", n_pt, 0.0, max_pt);
      book(_h["dy_tagjets"], "dy_tagjets", n_rap, 0, max_rap);
      book(_h["eta_tagjets"], "eta_tagjets", n_rap, -1*max_rap, max_rap);
      book(_h["dphi_tagjets"], "dphi_tagjets", n_rap, 0, max_rap);
      // //lepton plots
      book(_h["n_lepton_stable"], "n_lepton_stable", n_nbins, 0.0, n_nbins);
      book(_h["lepton_pt"], "lepton_pt", int(n_pt), 0.0, max_pt);
      book(_h["lepton_eta"], "lepton_eta", n_rap, -1*max_rap, max_rap);
      //photon plots
      book(_h["n_photons_iso"], "n_photons_iso", n_nbins, 0.0, n_nbins);
      book(_h["photon_iso_eta"], "photon_iso_eta", n_rap, -1*max_rap, max_rap);
      book(_h["photon_iso_pt"], "photon_iso_pt", int(n_pt), 0.0, max_pt);
      // lly
      book(_h["m_ly"], "m_ly", int(n_pt), 0.0, max_pt);
      book(_h["centrality_ly"], "centrality_ly", int(n_nbins), 0.0, n_nbins);
      //other
      book(_c["found_VBS_pair"],"found_VBS_pair");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Retrieve dressed leptons, sorted by pT
      Particles leptons_stable = apply<FinalState>(event, "leptons_stable").particlesByPt();
      int nlep_stable = leptons_stable.size();
      if (nlep_stable!=1)  vetoEvent; // meaning both are e,mu and not tau

      const FourMomentum lep = leptons_stable[0].mom();

      const MissingMomentum& METfinder = apply<MissingMomentum>(event, "METFinder");
      const double scalar_MET = METfinder.missingPt()/GeV;
      if (_docut==1 && scalar_MET<30.0) vetoEvent;

      const FourMomentum fourvec_MET = METfinder.missingMomentum();
      // const double mW_T = (fourvec_MET + lep).mass()/GeV;
      const double d_phi_MET_lep = deltaPhi(fourvec_MET.phi(), lep.phi());
      const double mW_T = sqrt( 2*lep.pT()*scalar_MET*(1 - cos(d_phi_MET_lep)) );
      if (_docut==1 && mW_T<30.0) vetoEvent;

      //photons
      Particles photons = apply<FinalState>(event, "photons").particlesByPt();
      if (photons.empty())  vetoEvent;
      //photon cone calculation and photon OR with leptons 
      Particles isolated_photons;
      Particles cone_sum_particles = apply<VetoedFinalState>(event, "ElAndJetsForPhotonIsoCalc").particles();
      for (const Particle &i_high_pt_photon : photons) {
        // check photon isolation
        double i_high_pt_photon_cone_E = 0.0;
        for (const Particle &i_p_cone : cone_sum_particles) {
          if (deltaR(i_high_pt_photon, i_p_cone) < 0.4) { // etcone40
            i_high_pt_photon_cone_E += i_p_cone.Et();
          }
        }
        if (i_high_pt_photon_cone_E > 0.2*i_high_pt_photon.pT())  continue;
        if (any(leptons_stable, deltaRLess(i_high_pt_photon, 0.4))) continue;
        isolated_photons += i_high_pt_photon;
      }
      if (isolated_photons.empty())  vetoEvent;

      const Particle& lead_iso_photon = isolated_photons[0];

      const FourMomentum fourvec_ly = lep + lead_iso_photon.mom();
      const double m_ly = fourvec_ly.mass();
      double m_z = 91.18; 
      if (_docut==1 && fabs(m_ly-m_z)<10*GeV)  vetoEvent;

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 20*GeV); // then will do cut on two leading pt>50
      Jets btagging_jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::absrap < 2.5 && Cuts::pT > 20*GeV);
      // Remove all jets within certain dR of a dressed lepton
      idiscardIfAnyDeltaRLess(jets, leptons_stable, 0.4);
      idiscardIfAnyDeltaRLess(jets, isolated_photons, 0.4);
      //
      idiscardIfAnyDeltaRLess(btagging_jets, leptons_stable, 0.4);
      idiscardIfAnyDeltaRLess(btagging_jets, isolated_photons, 0.4);

      int n_jets = jets.size();
      if (_docut==1 && n_jets < 2)  vetoEvent;  

      int nbtags = count(btagging_jets, hasBTag());
      if (_docut==1 && nbtags>0) vetoEvent;

      const FourMomentum tag1_jet = jets[0].mom();
      const FourMomentum tag2_jet = jets[1].mom();
      if (_docut==1 && (tag1_jet.pT()<50.0 || tag2_jet.pT()<50.0)) vetoEvent; 
      if (_docut==1 && deltaR(tag1_jet, tag2_jet) < 0.4) vetoEvent; 
      double dphi_tag1_MET = fabs(deltaPhi(tag1_jet,fourvec_MET));
      double dphi_tag2_MET = fabs(deltaPhi(tag2_jet,fourvec_MET));
      if (_docut==1 && (dphi_tag1_MET<0.4 || dphi_tag2_MET<0.4) ) vetoEvent; 

      const double m_tagjets = (tag1_jet + tag2_jet).mass()/GeV;
      if (_docut==1 && m_tagjets<1000.0) vetoEvent;

      const double dy_tagjets =  fabs(deltaRap(tag1_jet, tag2_jet));
      if (_docut==1 && dy_tagjets<2.0) vetoEvent;

      const double centrality_ly = fabs(0.5 * (fourvec_ly.rap() - (tag1_jet.rap()+tag2_jet.rap())/2) / (tag1_jet.rap()-tag2_jet.rap()));
      if (_docut==1 && centrality_ly > 0.35)  vetoEvent;

      int n_gap_jets = 0;
      for (int i = 0; i < n_jets; i++) {
        const double i_jet_rap = jets[i].rap();
        if ((i_jet_rap < tag1_jet.rap() && i_jet_rap > tag2_jet.rap()) || (i_jet_rap < tag2_jet.rap() && i_jet_rap > tag1_jet.rap()))  ++n_gap_jets;
      }
      if (_docut==1 && n_gap_jets > 0)  vetoEvent;

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
      _h["n_lepton_stable"]->fill(nlep_stable);
      _h["lepton_pt"]->fill(lep.pT()); 
      _h["lepton_eta"]->fill(lep.eta());      
      //photon plots
      _h["n_photons_iso"]->fill(isolated_photons.size());
      for (long unsigned int i = 0; i < isolated_photons.size(); i++) {
        _h["photon_iso_eta"]->fill(isolated_photons[i].eta());
        _h["photon_iso_pt"]->fill(isolated_photons[i].pT());
        }
      // ly
      _h["m_ly"]->fill(m_ly);
      _h["centrality_ly"]->fill(centrality_ly);
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


  RIVET_DECLARE_PLUGIN(Wmy_lvy);

}
