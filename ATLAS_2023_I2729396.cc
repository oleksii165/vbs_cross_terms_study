// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"


namespace Rivet {

  /// @brief Electroweak same-sign WW production at 13 TeV
  class ATLAS_2023_I2729396 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2023_I2729396);

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
        _cut_mode = getOption("cut");


      // Projection for dressed electrons and muons
      FinalState photons(Cuts::abspid == PID::PHOTON);

      PromptFinalState electrons(Cuts::abspid == PID::ELECTRON, false);
      DressedLeptons dressed_electrons(photons, electrons, 0.1, Cuts::open(), false, false);
      declare(dressed_electrons, "DressedElectrons");

      PromptFinalState muons(Cuts::abspid == PID::MUON, false); // true = use muons from prompt tau decays
      DressedLeptons dressed_muons(photons, muons, 0.1, Cuts::open(), false, false);
      declare(dressed_muons, "DressedMuons");

      // Jets, anti-kt 0.4
      VetoedFinalState fsJets; //final state for jet finding: veto leptons and neutrinos
      // match AntiKt4TruthJets:
      // remove neutrinos and undressed muons
      declare(FastJets(fsJets, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE), "Jets");

      // Missing momentum from invisibles
      declare(InvisibleFinalState(), "MET");
      // Histograms --- enforce to have the sames in SR and CR

      book(_h["mll"],      1, 1, 1);
      book(_h["mt"],       2, 1, 1);
      book(_h["mjj"],      3, 1, 1);
      book(_h["ngapjets"], 4, 1, 1);
      book(_h["zeppj3"],   5, 1, 1);

      book(_h["m_diboson"],"m_diboson", 600, 0.0, 10000.0);
      // last edge kind of irrelevant as overflow will be merged there
        const std::vector<double> m_ll_SM_bins = {20, 80, 130, 170, 220, 320, 500};
        book(_h["m_ll_SM_clip_inf"],"m_ll_SM_clip_inf", m_ll_SM_bins);
        book(_h["m_ll_SM_clip_3000"],"m_ll_SM_clip_3000", m_ll_SM_bins);
        book(_h["m_ll_SM_clip_2000"],"m_ll_SM_clip_2000", m_ll_SM_bins);
        book(_h["m_ll_SM_clip_1500"],"m_ll_SM_clip_1500", m_ll_SM_bins);
        book(_h["m_ll_SM_clip_1000"],"m_ll_SM_clip_1000", m_ll_SM_bins);
        book(_h["m_ll_SM_clip_700"],"m_ll_SM_clip_700", m_ll_SM_bins);
        const std::vector<double> m_ll_bins = {0, 250, 500, 750, 1000, 1500};
        book(_h["m_ll_clip_inf"],"m_ll_clip_inf", m_ll_bins);
        book(_h["m_ll_clip_3000"],"m_ll_clip_3000", m_ll_bins);
        book(_h["m_ll_clip_2000"],"m_ll_clip_2000", m_ll_bins);
        book(_h["m_ll_clip_1500"],"m_ll_clip_1500", m_ll_bins);
        book(_h["m_ll_clip_1000"],"m_ll_clip_1000", m_ll_bins);
        book(_h["m_ll_clip_700"],"m_ll_clip_700", m_ll_bins);


        book(_c["pos_w_initial"],"pos_w_initial");
        book(_c["neg_w_initial"],"neg_w_initial");
        for (std::string& i_clip : _clips){
            std::string i_pos_c_name = "pos_w_final_clip_" + i_clip;
            std::string i_neg_c_name = "neg_w_final_clip_" + i_clip;
            book(_c[i_pos_c_name], i_pos_c_name);
            book(_c[i_neg_c_name], i_neg_c_name);
        }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
        // save weights before cuts
        double ev_nominal_weight =  event.weights()[0];
        if (ev_nominal_weight>=0){_c["pos_w_initial"]->fill();}
        else {_c["neg_w_initial"]->fill();}

      const float mZ = 91.1876*GeV;

      // MET retrieving
      Vector3 met_vec;
      for (const Particle& p : apply<InvisibleFinalState>(event,"MET").particles()) {
        met_vec += p.mom().perpVec();
      }

      // MET selection
      if (met_vec.mod() < 30*GeV) vetoEvent;

      // Get jets
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 4.5);
      Jets btagging_jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 20*GeV && Cuts::abseta < 2.5);

      size_t nbtags = count(btagging_jets, hasBTag(Cuts::pT>5*GeV));
      if (nbtags > 0) vetoEvent;

      // Get electrons (sorted by pT)
      Particles electrons = apply<DressedLeptons>(event, "DressedElectrons").particlesByPt();
      // Get muons (sorted by pT)
      Particles muons = apply<DressedLeptons>(event, "DressedMuons").particlesByPt();

      // 3rd lepton veto:
      size_t n_leptonsVeto = select(electrons, Cuts::pT > 3*GeV && Cuts::abseta < 2.5).size();
      n_leptonsVeto += select(muons, Cuts::pT > 3*GeV && Cuts::abseta < 2.5).size();
      if (n_leptonsVeto > 2) vetoEvent;

      // Remove jet if lepton too close
      idiscardIfAny(jets, electrons, [](const Jet& jet, const Particle& lep) {
        return deltaR(lep, jet) < 0.4 && (lep.pT() / jet.pT() > 0.5);
      });

      // Get signal electrons
      iselect(electrons, Cuts::pT > 27*GeV && Cuts::abseta < 2.5);
      idiscardIfAny(electrons, jets, [](const Particle& lep, const Jet& jet) {
        return deltaR(lep, jet) < 0.4 && (lep.pT() / jet.pT() < 0.5);
      });

      // Get signal muons
      iselect(muons, Cuts::pT > 27*GeV && Cuts::abseta < 2.5);

      // Veto event if we don't have enough leptons
      if (electrons.size() + muons.size() < 2) vetoEvent;

      Particles leptons = electrons + muons;
      const Particle& lep1 = leptons[0];
      const Particle& lep2 = leptons[1];

      //electron rapidity in ee is < 1.37
      if (electrons.size() == 2 && (electrons[0].abseta() > 1.37 || electrons[1].abseta() > 1.37)) vetoEvent;

      // Require same-sign leading leptons
      if (lep1.charge()*lep2.charge() < 0.)  vetoEvent;

      // Delta R,ll cut
      if (deltaR(leptons[0], leptons[1]) < 0.2)  vetoEvent;

      // Mll cut
      if ((lep1.momentum() + lep2.momentum()).mass() < 20*GeV) vetoEvent;

      //Z veto in ee channel:
      if(electrons.size() == 2 && abs((lep1.mom() + lep2.mom()).mass() - mZ) < 15.*GeV) vetoEvent;

      auto dilep = lep1.mom() + lep2.mom();

      // Now jet cuts
      if (jets.size() < 2) vetoEvent;
      if (jets[0].momentum().pt() < 65.*GeV) vetoEvent;
      if (jets[1].momentum().pt() < 35.*GeV) vetoEvent;

      // VBS phase space
      const Jet& jet1 = jets[0];
      const Jet& jet2 = jets[1];
      double mjj = (jet1.mom() + jet2.mom()).mass();
      if (mjj < 500*GeV) vetoEvent;

      double deltayjj = deltaRap(jet1.mom(),jet2.mom());
      if (deltayjj < 2) vetoEvent;


        // do clipping - in case with 2+ w to avoid much work take two w with highest pt
        std::vector<FourMomentum> hs_bosons = {}; // also can be WZ for WZ CR
        const Particles all_particles = event.allParticles();
        for(const Particle& p_rivet : all_particles){
            ConstGenParticlePtr p_hepmc = p_rivet.genParticle();
            int status = p_hepmc->status();
            if (abs(status)==23 or abs(status)==22){
                int i_pid = p_hepmc->pid();
                FourMomentum i_mom = p_hepmc->momentum();
                if (abs(i_pid) == 24 or abs(i_pid) == 23){
                    hs_bosons.push_back(i_mom);
                }
            }
        }
        std::sort(hs_bosons.begin(), hs_bosons.end(), [](FourMomentum const &a, FourMomentum const &b) {return a.pT() > b.pT(); }); // biggest pT will be first in array
        bool have_two_hs_bosons = false;
        double hs_diboson_mass = 0.0;
        int n_z_hs = hs_bosons.size();
        if (n_z_hs > 1){
              hs_diboson_mass = (hs_bosons[0] + hs_bosons[1]).mass() / GeV;
              have_two_hs_bosons = true;
        }
        if (!have_two_hs_bosons) vetoEvent; // just in case reject events where dont have ww somehow

        _h["m_diboson"]->fill(hs_diboson_mass);
        _h["m_ll_clip_inf"]->fill(dilep.mass()/GeV);
        _h["m_ll_SM_clip_inf"]->fill(dilep.mass()/GeV);
        // _h["m_tagjets_clip_inf"]->fill(m_tagjets);
        if (ev_nominal_weight>=0){_c["pos_w_final_clip_inf"]->fill();}
        else {_c["neg_w_final_clip_inf"]->fill();}
        for (std::string& i_clip : _clips){
            if (i_clip=="inf") continue; // as done above without cut
            int i_clip_num = std::stoi(i_clip);
            std::string i_pos_c_name = "pos_w_final_clip_" + i_clip;
            std::string i_neg_c_name = "neg_w_final_clip_" + i_clip;
            std::string i_hist_name_m_ll_SM = "m_ll_SM_clip_" + i_clip;
            std::string i_hist_name_m_ll = "m_ll_clip_" + i_clip;
            // std::string i_hist_name_m_tagjets = "m_tagjets_clip_" + i_clip;
            if (hs_diboson_mass < i_clip_num) {
                _h[i_hist_name_m_ll_SM]->fill(dilep.mass()/GeV);
                _h[i_hist_name_m_ll]->fill(dilep.mass()/GeV);
                // _h[i_hist_name_m_tagjets]->fill(m_tagjets);
                if (ev_nominal_weight>=0){_c[i_pos_c_name]->fill();}
                else {_c[i_neg_c_name]->fill();}
            }
        }

      // original Histograms
      // mll
      _h["mll"]->fill(dilep.mass()/GeV);
      // mt
      double etll = sqrt(dilep.pt2() + dilep.mass2());
      double mt = sqrt(sqr(etll + met_vec.mod()) - (dilep.perpVec() + met_vec.perpVec()).mod2());
      _h["mt"]->fill(mt/GeV);
      // m_jj
      _h["mjj"]->fill(mjj/GeV);
      // zeppj3
      double avgj12 = 0.5*(jet1.eta()+jet2.eta());
      double norm12 = abs(jet1.eta()-jet2.eta());
      if (jets.size() > 2) {
        double z = abs(jets[2].eta() - avgj12)/norm12;
        _h["zeppj3"]->fill(z);
      }
      // ngapjets
      double ngapjets = 0;
      for (size_t i=2; i<jets.size(); i++) {
        if( (jets[i].eta() > jet1.eta() && jets[i].eta() < jet2.eta()) ||
            (jets[i].eta() > jet2.eta() && jets[i].eta() < jet1.eta()) ) {
          ++ngapjets;
        }
      }
      _h["ngapjets"]->fill(ngapjets);
    }

    /// Normalise histograms etc after the run
    void finalize() {
      for (auto& item : _h) {
        scale(item.second, crossSection()/sumOfWeights());
        // merge overflow bin into the last bin
        item.second->fillBin(item.second->numBins()-1, item.second->overflow().sumW());
        item.second->overflow().reset();
      }

      // in the paper, the differential cross section is per 10 GeV for mll and mt and per 100 GeV for mjj
      scale(_h["mll"], 10);
      scale(_h["mt"], 10);
      scale(_h["mjj"], 100);
    }

    //@}

    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, CounterPtr> _c;
    std::vector<std::string> _clips {"inf", "3000", "2000", "1500", "1000", "700"};
    std::string _cut_mode;
    //@}
  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2023_I2729396);
}
