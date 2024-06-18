// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Math/Constants.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"
#include "Rivet/AnalysisHandler.hh"
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace Rivet {


  /// @brief Wyjj production at 13 TeV
  class Wy_John : public Analysis {
  public:

      /// Constructor
      DEFAULT_RIVET_ANALYSIS_CTOR(Wy_John);


      /// @name Analysis methods
      ///@{

      /// Book histograms and initialise projections before the run
      void init() {

        _cut_mode = getOption("cut");
        // Initialise and register projections

        // The basic final-state projection:
        // all final-state particles within
        // the given eta acceptance
        const FinalState fs(Cuts::abseta < 5.0);

        //prompt photons
        const PromptFinalState photons_fs(Cuts::abspid == PID::PHOTON && Cuts::pT > 12 * GeV && Cuts::abseta < 2.37);
        declare(photons_fs, "photons");

        // Prompt leptons
        const PromptFinalState allphoton_fs(Cuts::abspid == PID::PHOTON, true,
                                            true); // photons used for lepton dressing
        const PromptFinalState bare_mu_fs(Cuts::abspid == PID::MUON, false);
        const PromptFinalState bare_e_fs(Cuts::abspid == PID::ELECTRON, false);
        const DressedLeptons dressedelectron_fs(allphoton_fs, bare_e_fs, 0.1);

        const DressedLeptons dressedmuon_fs(allphoton_fs, bare_mu_fs, 0.1);


        declare(dressedmuon_fs, "muons_dressed");
        declare(dressedelectron_fs, "electrons_dressed");

        //Get Leptons for clustering
        // Muons
        const PromptFinalState photon_cluster_fs(Cuts::abspid == PID::PHOTON, false,
                                                 false); // photons used for leptons to ignore in clustering
        PromptFinalState bare_mu(Cuts::abspid == PID::MUON, false, false); // true = use muons from prompt tau decays
        DressedLeptons all_dressed_mu(photon_cluster_fs, bare_mu, 0.1, Cuts::abseta < 5, false, false);
        //
        // // Electrons
        PromptFinalState bare_el(Cuts::abspid == PID::ELECTRON, false,
                                 false); // true = use electrons from prompt tau decays
        DressedLeptons all_dressed_el(photon_cluster_fs, bare_el, 0.1, Cuts::abseta < 5, false, false);

        PromptFinalState invis(Cuts::abspid == PID::NU_E ||
                               Cuts::abspid == PID::NU_MU ||
                               Cuts::abspid == PID::NU_TAU);


        VetoedFinalState hadrons(FinalState(Cuts::abseta < 5.0));
        hadrons.addVetoOnThisFinalState(all_dressed_el);
        hadrons.addVetoOnThisFinalState(all_dressed_mu);
        hadrons.addVetoOnThisFinalState(invis);
        declare(hadrons, "hadrons");

        // Project jets
        FastJets jets(hadrons, FastJets::ANTIKT, 0.4, FastJets::Muons::ALL, FastJets::Invisibles::DECAY);
        declare(jets, "jets");

        // Get MET from generic invisibles
        FinalState invis_met(Cuts::pT > 1e-9 * MeV && (Cuts::abspid == PID::NU_E ||
                                                       Cuts::abspid == PID::NU_MU ||
                                                       Cuts::abspid == PID::NU_TAU));
        declare(invis_met, "invis_met");

        // FS excluding the leading photon
        VetoedFinalState vfs;
        vfs.vetoNeutrinos();
        declare(vfs, "isolatedFS");

        // Book histograms
        const std::vector<double> lep_pt_bins = {30, 43, 60, 85, 130, 550};
        // last bin included only to exclude it later from integral (to have normalization same as they have)
        const std::vector<double> lep_pt_bins_with_underoverflow = {0, 30, 43, 60, 85, 130, 550, 20000};
        const std::vector<double> ly_m_bins = {0, 125, 180, 245, 350, 1200};
        const std::vector<double> lepgam_dphi_signed_bins = {-3.1416, -2.55, -2.0, -1.1, 0.0, 1.1, 2.0, 2.55, 3.1416};
        const std::vector<double> jj_m_bins = {1000, 1400, 1700, 2100, 5300};
        const std::vector<double> jj_pt_bins = {0, 80, 125, 180, 270, 840};
        const std::vector<double> jj_dphi_signed_bins = {-3.1416, -2.7489, -2.4435, 0.0, 2.4435, 2.7489, 3.1416};

        book(_hist_xs, "hist_xs", 1, 0., 1.);
        book(_hist_total_xs, "hist_total_xs", 1, 0., 1.);
        book(_hist_lep_pt, "hist_lep_pt", lep_pt_bins);
        book(_pt_lepton_clip_inf, "pt_lepton_clip_inf", lep_pt_bins_with_underoverflow);
        book(_pt_lepton_clip_3000, "pt_lepton_clip_3000", lep_pt_bins_with_underoverflow);
        book(_pt_lepton_clip_2000, "pt_lepton_clip_2000", lep_pt_bins_with_underoverflow);
        book(_pt_lepton_clip_1500, "pt_lepton_clip_1500", lep_pt_bins_with_underoverflow);
        book(_pt_lepton_clip_1000, "pt_lepton_clip_1000", lep_pt_bins_with_underoverflow);
        book(_pt_lepton_clip_700, "pt_lepton_clip_700", lep_pt_bins_with_underoverflow);
        book(_m_Wy, "m_Wy", 600, 0, 10000);
        book(_hist_ly_m, "hist_ly_m", ly_m_bins);
        book(_hist_lepgam_dphi_signed, "hist_lepgam_dphi_signed", lepgam_dphi_signed_bins);
        book(_hist_jj_m, "hist_jj_m", jj_m_bins);
        book(_hist_jj_pt, "hist_jj_pt", jj_pt_bins);
        book(_hist_jj_dphi_signed, "hist_jj_dphi_signed", jj_dphi_signed_bins);
        book(_hist_cutflow, "hist_cutflow", 25, -0.5, 24.5);

        //counter for efficiency
        book(_c["pos_w_initial"], "pos_w_initial");
        book(_c["neg_w_initial"], "neg_w_initial");
        for (std::string &i_clip: _clips) {
          std::string i_pos_c_name = "pos_w_final_clip_" + i_clip;
          std::string i_neg_c_name = "neg_w_final_clip_" + i_clip;
          book(_c[i_pos_c_name], i_pos_c_name);
          book(_c[i_neg_c_name], i_neg_c_name);
        }

      }


      /// Perform the per-event analysis
      void analyze(const Event &event) {
        double ev_nominal_weight = event.weights()[0];
        if (ev_nominal_weight >=
            0) { _c["pos_w_initial"]->fill(); } // dont need anything in bracket as this will be weight on weight
        else { _c["neg_w_initial"]->fill(); }

        const size_t numWeights = event.weights().size();
        for (size_t m = 0; m < numWeights; ++m) {
          _hist_cutflow.get()->_getPersistent(m)->fill(0, 1.0);
        }

        _hist_total_xs->fill(0.5);

        const Particles &photons = apply<PromptFinalState>(event, "photons").particlesByPt();
        const Jets &jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::absrap < 4.4 && Cuts::pT > 25 * GeV);
        const Particles &electrons = apply<ParticleFinder>(event, "electrons_dressed").particlesByPt();
        const Particles &muons = apply<ParticleFinder>(event, "muons_dressed").particlesByPt();

        Particles pre_selectedPh;
        Particles noniso_Ph;
        Particles pre_selectedElectrons;
        Particles pre_selectedMuons;
        Jets pre_selectedJets;
        Particles selectedPh;
        Particles final_selectedElectrons;
        Particles final_selectedMuons;
        Jets prefinal_selectedJets;
        Jets final_selectedJets;
        Particles final_selectedPh;

        // Lepton preselection
        for (const Particle &lep: electrons) {
          double eta = std::abs(lep.eta());
          if (eta < 2.47 && !(eta > 1.37 && eta < 1.52) && lep.pt() > 7.0) {
            pre_selectedElectrons.push_back(lep);
          }
        }
        for (const Particle &lep: muons) {
          double eta = std::abs(lep.eta());
          if (eta < 2.7 && lep.pt() > 7.0) {
            pre_selectedMuons.push_back(lep);
          }
        }

        Particles fs = apply<VetoedFinalState>(event, "isolatedFS").particles();

        // Photon preselection
        for (const Particle &ph: photons) {
          // check photon isolation
          double eta = std::abs(ph.eta());
          if (eta > 2.37 || (eta > 1.37 && eta < 1.52) || ph.pt() < 12) continue;
          float coneEnergy_scaler = 0;
          for (const Particle &p: fs) {
            // Do not count the photons energy
            if (isSame(ph, p) || !isStable(p)) {
              continue;
            }
            if (deltaR(ph, p) < 0.4) {
              coneEnergy_scaler += p.pt();
            }
          }
          if ((coneEnergy_scaler > 0.2 * ph.Et())) continue;
          pre_selectedPh.push_back(ph);

        }


        int init = 0;
        if (pre_selectedPh.size() > 0) {
          for (size_t m = 0; m < numWeights; ++m) {
            _hist_cutflow.get()->_getPersistent(m)->fill(1, 1.0);
          }
          if (pre_selectedMuons.size() + pre_selectedElectrons.size() > 0) {
            for (size_t m = 0; m < numWeights; ++m) {
              _hist_cutflow.get()->_getPersistent(m)->fill(2, 1.0);
            }
            if (jets.size() > 1) {
              for (size_t m = 0; m < numWeights; ++m) {
                _hist_cutflow.get()->_getPersistent(m)->fill(3, 1.0);
                init = 1;
              }
            }
          }
        }
        if (!init) vetoEvent;


        // Overlap removal
        // remove jets against photon and electrons
        for (const Jet &jet: jets) {
          bool overlap = false;
          for (const Particle &ph: pre_selectedPh) {
            if (deltaR(jet, ph, RAPIDITY) < 0.4) {
              overlap = true;
            }
          }
          for (const Particle &lep: pre_selectedElectrons) {
            if (deltaR(jet, lep, RAPIDITY) < 0.2) {
              overlap = true;
            }
          }
          if (!overlap) pre_selectedJets.push_back(jet);
        }
        // remove photons against electrons and muons
        for (const Particle &ph: pre_selectedPh) {
          bool overlap = false;
          for (const Particle &lep: pre_selectedElectrons) {
            if (deltaR(ph, lep) < 0.4) {
              overlap = true;
            }
          }
          for (const Particle &lep: pre_selectedMuons) {
            if (deltaR(ph, lep) < 0.4) {
              overlap = true;
            }
          }
          if (!overlap) selectedPh.push_back(ph);
        }
        // reject jets against electrons
        for (const Jet &jet: pre_selectedJets) {
          bool overlap = false;
          for (const Particle &lep: pre_selectedElectrons) {
            if (deltaR(jet, lep, RAPIDITY) < 0.2) {
              overlap = true;
            }
          }
          if (!overlap) prefinal_selectedJets.push_back(jet);
        }
        // reject electrons against jets
        for (const Particle &lep: pre_selectedElectrons) {
          bool overlap = false;
          for (const Jet &jet: prefinal_selectedJets) {
            if (deltaR(jet, lep, RAPIDITY) < 0.4) {
              overlap = true;
            }
          }
          if (!overlap) final_selectedElectrons.push_back(lep);
        }

        // reject jets against muons
        for (const Jet &jet: prefinal_selectedJets) {
          bool overlap = false;
          for (const Particle &lep: pre_selectedMuons) {
            if (deltaR(jet, lep, RAPIDITY) < 0.2) {
              overlap = true;
            }
          }
          if (!overlap) final_selectedJets.push_back(jet);
        }
        // reject muons against jets
        for (const Particle &lep: pre_selectedMuons) {
          bool overlap = false;
          for (const Jet &jet: final_selectedJets) {
            if (deltaR(jet, lep, RAPIDITY) < 0.4) {
              overlap = true;
            }
          }
          if (!overlap) final_selectedMuons.push_back(lep);
        }

        // remove photons against jets
        for (const Particle &ph: selectedPh) {
          bool overlap = false;
          for (const Jet &jet: final_selectedJets) {
            if (deltaR(jet, ph, RAPIDITY) < 0.4) {
              overlap = true;
            }
          }
          if (!overlap) final_selectedPh.push_back(ph);
        }

        int pass_or = 0;
        if (final_selectedPh.size() > 0) {
          for (size_t m = 0; m < numWeights; ++m) {
            _hist_cutflow.get()->_getPersistent(m)->fill(4, 1.0);
          }
          if (final_selectedMuons.size() + final_selectedElectrons.size() > 0) {
            for (size_t m = 0; m < numWeights; ++m) {
              _hist_cutflow.get()->_getPersistent(m)->fill(5, 1.0);
            }
            if (final_selectedJets.size() > 1) {
              for (size_t m = 0; m < numWeights; ++m) {
                _hist_cutflow.get()->_getPersistent(m)->fill(6, 1.0);
                pass_or = 1;
              }
            }
          }
        }
        if (!pass_or) vetoEvent;


        if (final_selectedPh.size() < 1) vetoEvent;
        if (final_selectedPh[0].pt() < 22) vetoEvent;
        for (size_t m = 0; m < numWeights; ++m) {
          _hist_cutflow.get()->_getPersistent(m)->fill(7, 1.0);
        }
        Particle final_selectedLep;
        if (final_selectedMuons.size() == 1) {
          final_selectedLep = final_selectedMuons[0];
        } else if (final_selectedElectrons.size() == 1) {
          final_selectedLep = final_selectedElectrons[0];
        }
        if (std::abs(final_selectedLep.eta()) > 2.5) vetoEvent;
        if (final_selectedLep.pt() < 30) vetoEvent;
        for (size_t m = 0; m < numWeights; ++m) {
          _hist_cutflow.get()->_getPersistent(m)->fill(8, 1.0);
        }


        // Apply mass cut on lepton-photon system
        double zMass = 91.1876; // Z boson mass
        double lep_ph_mass = (final_selectedLep.momentum() + final_selectedPh[0].momentum()).mass();
        if (std::abs(lep_ph_mass - zMass) < 10) vetoEvent;
        for (size_t m = 0; m < numWeights; ++m) {
          _hist_cutflow.get()->_getPersistent(m)->fill(9, 1.0);
        }

        const Particles &invis_met = apply<FinalState>(event, "invis_met").particlesByPt();
        double nu_px = 0;
        double nu_py = 0;
        for (const Particle &p: invis_met) {
          nu_px += p.px();
          nu_py += p.py();
        }
        double met = sqrt(nu_px * nu_px + nu_py * nu_py);
        Vector3 met_tvec(nu_px, nu_py, 0);
        Vector3 lep_tvec(final_selectedLep.px(), final_selectedLep.py(), 0);

        const double mt = sqrt(
                2.0 * (final_selectedLep.pt() / GeV) * (met / GeV) * (1. - cos(deltaPhi(lep_tvec, met_tvec))));
        if (met <= 30 * GeV) vetoEvent;
        for (size_t m = 0; m < numWeights; ++m) {
          _hist_cutflow.get()->_getPersistent(m)->fill(10, 1.0);
        }
        if (mt <= 30 * GeV) vetoEvent;
        for (size_t m = 0; m < numWeights; ++m) {
          _hist_cutflow.get()->_getPersistent(m)->fill(11, 1.0);
        }

        size_t nbtags = count(final_selectedJets, hasBTag(Cuts::abseta < 2.5));

        if (nbtags > 0) vetoEvent;
        for (size_t m = 0; m < numWeights; ++m) {
          _hist_cutflow.get()->_getPersistent(m)->fill(12, 1.0);
        }
        if (final_selectedJets.size() < 2) vetoEvent;
        for (size_t m = 0; m < numWeights; ++m) {
          _hist_cutflow.get()->_getPersistent(m)->fill(13, 1.0);
        }
        double dijetMass = (final_selectedJets[0].momentum() + final_selectedJets[1].momentum()).mass();
        double dijetpT = (final_selectedJets[0].momentum() + final_selectedJets[1].momentum()).pT();
        if (final_selectedJets[0].pt() < 50.0 || final_selectedJets[1].pt() < 50.0) vetoEvent;
        for (size_t m = 0; m < numWeights; ++m) {
          _hist_cutflow.get()->_getPersistent(m)->fill(14, 1.0);
        }
        if (dijetMass < 1000) vetoEvent;
        for (size_t m = 0; m < numWeights; ++m) {
          _hist_cutflow.get()->_getPersistent(m)->fill(15, 1.0);
        }
        // Additional dR cuts
        if (deltaR(final_selectedJets[0], final_selectedJets[1], RAPIDITY) < 0.4) vetoEvent;
        if (deltaR(final_selectedJets[0], final_selectedLep, RAPIDITY) < 0.4) vetoEvent;
        if (deltaR(final_selectedJets[1], final_selectedLep, RAPIDITY) < 0.4) vetoEvent;
        if (deltaR(final_selectedPh[0], final_selectedJets[0], RAPIDITY) < 0.4) vetoEvent;
        if (deltaR(final_selectedPh[0], final_selectedJets[1], RAPIDITY) < 0.4) vetoEvent;
        if (deltaR(final_selectedPh[0], final_selectedLep, RAPIDITY) < 0.4) vetoEvent;
        for (size_t m = 0; m < numWeights; ++m) {
          _hist_cutflow.get()->_getPersistent(m)->fill(16, 1.0);
        }
        if (deltaPhi(met_tvec, final_selectedJets[0]) < 0.4) vetoEvent;
        if (deltaPhi(met_tvec, final_selectedJets[1]) < 0.4) vetoEvent;
        for (size_t m = 0; m < numWeights; ++m) {
          _hist_cutflow.get()->_getPersistent(m)->fill(17, 1.0);
        }

        int n_gapjets = 0;
        double jj_drap = abs(final_selectedJets[0].rapidity() - final_selectedJets[1].rapidity());
        for (unsigned int i = 2; i < final_selectedJets.size(); i++) {
          double _dRapj0 = abs(final_selectedJets[0].rapidity() - final_selectedJets[i].rapidity());
          double _dRapj1 = abs(final_selectedJets[1].rapidity() - final_selectedJets[i].rapidity());
          if ((_dRapj0 < jj_drap) && (_dRapj1 < jj_drap)) ++n_gapjets;
        }
        if (jj_drap < 2) vetoEvent;
        for (size_t m = 0; m < numWeights; ++m) {
          _hist_cutflow.get()->_getPersistent(m)->fill(18, 1.0);
        }

        FourMomentum lepgam = final_selectedLep.momentum() + final_selectedPh[0].momentum();
        double lepgam_xi =
                abs(lepgam.rapidity() - 0.5 * (final_selectedJets[0].rapidity() + final_selectedJets[1].rapidity())) /
                jj_drap;
        double ly_m = lepgam.mass();

        if (n_gapjets > 0) vetoEvent;
        for (size_t m = 0; m < numWeights; ++m) {
          _hist_cutflow.get()->_getPersistent(m)->fill(19, 1.0);
        }
        if (lepgam_xi > 0.35) vetoEvent;
        for (size_t m = 0; m < numWeights; ++m) {
          _hist_cutflow.get()->_getPersistent(m)->fill(20, 1.0);
        }

        //calculate dphi signed
        double phi_fwd_jj = ((final_selectedJets[0].momentum()).rapidity() >
                             (final_selectedJets[1].momentum()).rapidity()) ? (final_selectedJets[0].momentum()).phi()
                                                                            : (final_selectedJets[1].momentum()).phi();
        double phi_bwd_jj = ((final_selectedJets[0].momentum()).rapidity() <=
                             (final_selectedJets[1].momentum()).rapidity()) ? (final_selectedJets[0].momentum()).phi()
                                                                            : (final_selectedJets[1].momentum()).phi();
        double jj_dphi_signed = mpi_pi(phi_fwd_jj - phi_bwd_jj);

        double phi_fwd_lepgam = ((final_selectedLep.momentum()).rapidity() >
                                 (final_selectedPh[0].momentum()).rapidity()) ? (final_selectedLep.momentum()).phi()
                                                                              : (final_selectedPh[0].momentum()).phi();
        double phi_bwd_lepgam = ((final_selectedLep.momentum()).rapidity() <=
                                 (final_selectedPh[0].momentum()).rapidity()) ? (final_selectedLep.momentum()).phi()
                                                                              : (final_selectedPh[0].momentum()).phi();
        double lepgam_dphi_signed = mpi_pi(phi_fwd_lepgam - phi_bwd_lepgam);
        double lep_pt = final_selectedLep.Et() / GeV;

        _hist_xs->fill(0.5);

        _hist_lep_pt->fill(lep_pt);
        _hist_ly_m->fill(ly_m);
        _hist_lepgam_dphi_signed->fill(lepgam_dphi_signed);
        _hist_jj_m->fill(dijetMass);
        _hist_jj_pt->fill(dijetpT);
        _hist_jj_dphi_signed->fill(jj_dphi_signed);

        // do clipping - in case with 2+ bosons to avoid much work take two w with highest pt
        std::vector <FourMomentum> hs_bosons = {};
        const Particles all_particles = event.allParticles();
        for (const Particle &p_rivet: all_particles) {
          ConstGenParticlePtr p_hepmc = p_rivet.genParticle();
          int status = p_hepmc->status();
          if (abs(status) == 23 or abs(status) == 22) {
            int i_pid = p_hepmc->pid();
            FourMomentum i_mom = p_hepmc->momentum();
            if (abs(i_pid) == 24 or abs(i_pid) == 22) { hs_bosons.push_back(i_mom); }
          }
        }
        std::sort(hs_bosons.begin(), hs_bosons.end(),
                  [](FourMomentum const &a, FourMomentum const &b) { return a.pT() > b.pT(); });
        bool have_two_hs_bosons = false;
        double hs_diboson_mass = 0.0;
        int n_hs = hs_bosons.size();
        if (n_hs > 1) {
          hs_diboson_mass = (hs_bosons[0] + hs_bosons[1]).mass() / GeV;
          have_two_hs_bosons = true;
        }
        if (!have_two_hs_bosons) vetoEvent; // just in case reject events where dont have wz somehow

        _m_Wy->fill(hs_diboson_mass);
        if (ev_nominal_weight >= 0) { _c["pos_w_final_clip_inf"]->fill(); }
        else { _c["neg_w_final_clip_inf"]->fill(); }
        for (std::string &i_clip: _clips) {
          if (i_clip == "inf") continue; // as done above without cut
          int i_clip_num = std::stoi(i_clip);
          std::string i_pos_c_name = "pos_w_final_clip_" + i_clip;
          std::string i_neg_c_name = "neg_w_final_clip_" + i_clip;
          if (hs_diboson_mass < i_clip_num) {
            if (ev_nominal_weight >= 0) { _c[i_pos_c_name]->fill(); }
            else { _c[i_neg_c_name]->fill(); }
          }
        }

        _pt_lepton_clip_inf->fill(lep_pt);
        if (hs_diboson_mass < 3000) { _pt_lepton_clip_3000->fill(lep_pt); }
        if (hs_diboson_mass < 2000) { _pt_lepton_clip_2000->fill(lep_pt); }
        if (hs_diboson_mass < 1500) { _pt_lepton_clip_1500->fill(lep_pt); }
        if (hs_diboson_mass < 1000) { _pt_lepton_clip_1000->fill(lep_pt); }
        if (hs_diboson_mass < 700) { _pt_lepton_clip_700->fill(lep_pt); }
      }



    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_hist_xs,1.0); // last argument is inclusion or not of overflow to normalization
      normalize(_hist_total_xs,1.0);
      normalize(_hist_lep_pt,1.0);
      normalize(_hist_ly_m,1.0);
      normalize(_hist_lepgam_dphi_signed,1.0);
      normalize(_hist_jj_m,1.0);
      normalize(_hist_jj_pt,1.0);
      normalize(_hist_jj_dphi_signed,1.0);

      normalize(_pt_lepton_clip_inf,1.0);
      normalize(_pt_lepton_clip_3000,1.0);
      normalize(_pt_lepton_clip_2000,1.0);
      normalize(_pt_lepton_clip_1500,1.0);
      normalize(_pt_lepton_clip_1000,1.0);
      normalize(_pt_lepton_clip_700,1.0);
      normalize(_m_Wy,1.0);

      }


    double mpi_pi(double angle){

      while (angle >= pi) angle -= 2*pi;
      while (angle < -1*pi) angle += 2*pi;

      return angle;

    }


    private:
    
      /// Histograms
      Histo1DPtr _hist_xs,_hist_total_xs,_hist_lep_pt,_pt_lepton_clip_inf,_pt_lepton_clip_3000,_pt_lepton_clip_2000,_pt_lepton_clip_1500,_pt_lepton_clip_1000,_pt_lepton_clip_700,_m_Wy,_hist_ly_m, _hist_lepgam_dphi_signed,_hist_jj_m,_hist_jj_pt,_hist_mtw,_hist_jj_dphi_signed,_hist_cutflow;
      std::string _cut_mode;
      map<string, CounterPtr> _c;
      std::vector<std::string> _clips {"inf", "3000", "2000", "1500", "1000", "700"};

  };


  DECLARE_RIVET_PLUGIN(Wy_John);

}

