#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"

namespace Rivet {

  /// @brief VBS Zy production at 13 TeV
  class ATLAS_2023_I2663725 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2023_I2663725);

    /// @name Analysis methods
    //@{
    /// Book histograms and initialise projections before the run
    void init() {
      _cut_mode = getOption("cut");

      FinalState fs;

      // Prompt photons
      const PromptFinalState photon_fs(Cuts::abspid == PID::PHOTON && Cuts::pT > 25 * GeV && Cuts::abseta < 2.37);
      declare(photon_fs, "Photons");

      // Prompt leptons
      const PromptFinalState bareelectron_fs = Cuts::abspid == PID::ELECTRON;
      const PromptFinalState baremuon_fs = Cuts::abspid == PID::MUON;

      // Dressed leptons
      const FinalState allphoton_fs(Cuts::abspid == PID::PHOTON);
      const Cut leptoncut = Cuts::pT > 20 * GeV && Cuts::abseta < 2.47;
      const DressedLeptons dressedelectron_fs(allphoton_fs, bareelectron_fs, 0.1, leptoncut, true); // use *all* photons for lepton dressing
      const DressedLeptons dressedmuon_fs(allphoton_fs, baremuon_fs, 0.1, leptoncut, true);         // use *all* photons for lepton dressing

      declare(dressedelectron_fs, "Electrons");
      declare(dressedmuon_fs, "Muons");

      // FS excluding the leading photon
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(photon_fs);
      vfs.addVetoOnThisFinalState(dressedmuon_fs);
      vfs.addVetoOnThisFinalState(InvisibleFinalState());
      declare(vfs, "isolatedFS");

      VetoedFinalState hadrons(FinalState(Cuts::abseta < 4.4));
      hadrons.addVetoOnThisFinalState(dressedelectron_fs);
      hadrons.addVetoOnThisFinalState(dressedmuon_fs);
      declare(hadrons, "hadrons");

      // Jets
      FastJets jets(hadrons, FastJets::ANTIKT, 0.4, JetAlg::Muons::ALL, JetAlg::Invisibles::DECAY);
      declare(jets, "jets");

      // bins for hists
      const std::vector<double> Dphi_lly_jj_EWK_bins = {0, 2.96, 3.2};
      const std::vector<double> pT_y1_EWK_bins = {25,40,100,500};
      const std::vector<double> pT_lly_EWK_bins = {0,100,700};

      // Histograms
      //note:
      // - all the hist without _EWK denotes dist. in Extented SR (Centrality < 0.4)
      // - all the hists with _EWK denotes dist. in SR (Centrality < 0.4 && mjj>500)
      book(_h["NJets"], "NJets", {0,1,2,3,4,5,10});

      book(_h["pT_j1"], "pT_j1", {50,100,175,300,500,1000});
      book(_h["pT_j1_EWK"], "pT_j1_EWK", {50,150,250,1000});//
      book(_h["pT_j2"], "pT_j2", {50,100,175,300,500,1000});
      book(_h["pT_y1"], "pT_y1", {25,35,50,75,150,500});
      book(_h["pT_y1_EWK"], "pT_y1_EWK", pT_y1_EWK_bins);//
      book(_h["pT_l1"], "pT_l1", {30,65,100,150,250,500});
      book(_h["pT_l1_EWK"], "pT_l1_EWK", {30,150,1000});//
      book(_h["pT_l2"], "pT_l2", {30,65,100,150,250,500});
      book(_h["pT_ll"], "pT_ll", {0,50,100,150,250,800});
      book(_h["pT_ll_EWK"], "pT_ll_EWK", {0,50,100,150,250,800});//
      book(_h["pT_lly"], "pT_lly", {0,50,100,150,200,500});
      book(_h["pT_lly_EWK"], "pT_lly_EWK", pT_lly_EWK_bins);//

      book(_h["eta_j1"], "eta_j1", {0,0.5,1,1.5,2,3,4,5});
      book(_h["eta_j2"], "eta_j2", {0,0.5,1,1.5,2,3,4,5});
      book(_h["eta_y1"], "eta_y1", {0,0.5,1,1.5,2,3,4,5});
      book(_h["eta_l1"], "eta_l1", {0,0.5,1,1.5,2,3,4,5});
      book(_h["eta_l2"], "eta_l2", {0,0.5,1,1.5,2,3,4,5});

      book(_h["phi_j1"], "phi_j1", {0,0.5,1,1.5,2,2.5,3,3.5});
      book(_h["phi_j2"], "phi_j2", {0,0.5,1,1.5,2,2.5,3,3.5});
      book(_h["phi_y1"], "phi_y1", {0,0.5,1,1.5,2,2.5,3,3.5});
      book(_h["phi_l1"], "phi_l1", {0,0.5,1,1.5,2,2.5,3,3.5});
      book(_h["phi_l2"], "phi_l2", {0,0.5,1,1.5,2,2.5,3,3.5});

      book(_h["njGap"], "njGap", {0,1,2,3,5});
      book(_h["Centrality"], "Centrality", {0,0.075,0.15,0.225,0.3,0.4});
      book(_h["Centrality_EWK"], "Centrality_EWK", {0,0.075,0.15,0.225,0.3,0.4});//
      book(_h["Dyjj"], "Dyjj", {1,1.75,2.5,3.5,5,9});
      book(_h["Dyjj_EWK"], "Dyjj_EWK", {1,3.5,9});//
      book(_h["mjj"], "mjj", {150,250,500,1000,1600,2400,5000});
      book(_h["mjj_EWK"], "mjj_EWK", {450,500,1000,1600,2400,5000});// special since 450-500 in the underflow bin considered in analysis
      book(_h["mll"], "mll", {40,50,60,70,80,90,100,110,150,200});
      book(_h["mlly"], "mlly", {100,150,200,250,350,700});
      book(_h["mlly_EWK"], "mlly_EWK", {100,150,200,250,350,700});//
      book(_h["j1rap"], "j1rap", {0,0.5,1,1.5,2,3,4,5});
      book(_h["j2rap"], "j2rap", {0,0.5,1,1.5,2,3,4,5});
      book(_h["Dphi_lly_jj"], "Dphi_lly_jj", {0,2,2.95,3.2});
      book(_h["Dphi_lly_jj_EWK"], "Dphi_lly_jj_EWK", Dphi_lly_jj_EWK_bins);

      book(_h["m_Zy"], "m_Zy", 600, 0, 10000.0);

      // hists for EFT per clip
      for (std::string &i_clip: _clips) {
        std::string i_name = "Dphi_lly_jj_EWK_clip_" + i_clip;
        book(_h[i_name], i_name, Dphi_lly_jj_EWK_bins);
        //
        i_name = "pT_y1_EWK_clip_" + i_clip;
        book(_h[i_name], i_name, pT_y1_EWK_bins);
        //
        i_name = "pT_lly_EWK_clip_" + i_clip;
        book(_h[i_name], i_name, pT_lly_EWK_bins);
        }

      //counter for efficiency - in this ana not doing hists myself so can ignore clippings
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
      // save weights before cuts
      double ev_nominal_weight = event.weights()[0];
      if (ev_nominal_weight >=
          0) { _c["pos_w_initial"]->fill(); } // dont need anything in bracket as this will be weight on weight
      else { _c["neg_w_initial"]->fill(); }

      // Get objects
      Particles electrons = apply<DressedLeptons>(event, "Electrons").particlesByPt();
      Particles muons = apply<DressedLeptons>(event, "Muons").particlesByPt();
      Particles photons = apply<PromptFinalState>(event, "Photons").particlesByPt();

      if (photons.empty())  vetoEvent;
      if (electrons.size() < 2 && muons.size() < 2)  vetoEvent;

      Particles lep;
      if (electrons.size() >= 2) {
        lep.push_back(electrons[0]);
        lep.push_back(electrons[1]);
      }
      else {
        lep.push_back(muons[0]);
        lep.push_back(muons[1]);
      }

      if (lep[0].pT() < 30*GeV) vetoEvent;

      const double mll = (lep[0].momentum() + lep[1].momentum()).mass();
      if (mll < 40 * GeV)  vetoEvent;

      Particles selectedPh;
      Particles fs = apply<VetoedFinalState>(event, "isolatedFS").particles();
      for (const Particle &ph : photons) {
        // check photon isolation
        double coneEnergy = 0.0;
        for (const Particle &p : fs) {
          if (deltaR(ph, p) < 0.2) { // etcone20
            coneEnergy += p.Et();
          }
        }
        if (coneEnergy / ph.pT() > 0.07)  continue;
        if (any(electrons, deltaRLess(ph, 0.4))) continue;
        if (any(muons, deltaRLess(ph, 0.4))) continue;
        selectedPh += ph;
      }
      if (selectedPh.empty())  vetoEvent;

      FourMomentum lly = lep[0].mom() + lep[1].mom() + selectedPh[0].mom();
      const double mlly = lly.mass();

      if (mll + mlly <= 182*GeV)  vetoEvent;

      // Get jets
      const Cut jetscut = (Cuts::pT > 25*GeV && Cuts::absrap < 4.4);
      Jets jetsMix = apply<FastJets>(event, "jets").jetsByPt(jetscut);
      idiscardIfAnyDeltaRLess(jetsMix, photons, 0.4);
      idiscardIfAnyDeltaRLess(jetsMix, electrons, 0.3);
      idiscardIfAnyDeltaRLess(jetsMix, muons, 0.3);

      // Jet multiplicity
      size_t njets = jetsMix.size();
      if (njets < 2)  vetoEvent;

      if (jetsMix[1].pT() < 50 * GeV)  vetoEvent;

      // other complex cuts
      FourMomentum dijet = jetsMix[0].mom() + jetsMix[1].mom();
      const double mjj = dijet.mass();
      if (mjj < 150*GeV)  vetoEvent;

      double Dyjj = fabs(deltaRap(jetsMix[0], jetsMix[1]));
      if (Dyjj < 1)  vetoEvent;

      const double j1rap = jetsMix[0].rapidity();
      const double j2rap = jetsMix[1].rapidity();

      const double Centrality = fabs(0.5*(lly.rap() - (j1rap + j2rap)) / (j1rap - j2rap));
      if (Centrality > 5)  vetoEvent;

      size_t njGap = 0;
      for (size_t i=2; i < njets; ++i) {
        const double jrap = jetsMix[i].rap();
        if ((jrap < j1rap && jrap > j2rap) || (jrap < j2rap && jrap > j1rap))  ++njGap;
      }
      if (njGap > 0)  vetoEvent;

      // now start to log hist
      // first we select the extended SR
      if (Centrality > 0.4) vetoEvent;
      //// basic variables for check
      _h["NJets"]->fill(njets);
      _h["pT_j1"]->fill(jetsMix[0].pT() / GeV);
      _h["pT_j2"]->fill(jetsMix[1].pT() / GeV);
      _h["pT_y1"]->fill(selectedPh[0].pT() / GeV);
      _h["pT_l1"]->fill(lep[0].pT() / GeV);
      _h["pT_l2"]->fill(lep[0].pT() / GeV);
      _h["pT_ll"]->fill((lep[0].mom() + lep[1].mom()).pT() / GeV);
      _h["pT_lly"]->fill(lly.pT() / GeV);

      _h["eta_j1"]->fill(jetsMix[0].eta());
      _h["eta_j2"]->fill(jetsMix[1].eta());
      _h["eta_y1"]->fill(selectedPh[0].eta());
      _h["eta_l1"]->fill(lep[0].eta());
      _h["eta_l2"]->fill(lep[0].eta());

      _h["phi_j1"]->fill(jetsMix[0].phi());
      _h["phi_j2"]->fill(jetsMix[1].phi());
      _h["phi_y1"]->fill(selectedPh[0].phi());
      _h["phi_l1"]->fill(lep[0].phi());
      _h["phi_l2"]->fill(lep[0].phi());

      //// some cut variables

      _h["njGap"]->fill(njGap);
      _h["Centrality"]->fill(Centrality);
      _h["Dyjj"]->fill(Dyjj);
      _h["mjj"]->fill(mjj);
      _h["mll"]->fill(mll);
      _h["mlly"]->fill(mlly);
      _h["j1rap"]->fill(j1rap);
      _h["j2rap"]->fill(j2rap);

      //// other unfolded variables //
      _h["Dphi_lly_jj"]->fill(fabs(deltaPhi(lly, dijet)));

			// deal with EWK SR500
			if (mjj < 450*GeV)  vetoEvent;
			_h["mjj_EWK"]->fill(mjj);

			if (mjj < 500*GeV)  vetoEvent;

      double my_dphi = fabs(deltaPhi(lly, dijet));
      double my_pty = fabs(selectedPh[0].pT() / GeV);
      double my_ptZy = fabs(lly.pT() / GeV);

			_h["pT_j1_EWK"]->fill(jetsMix[0].pT() / GeV);
			_h["pT_y1_EWK"]->fill(my_pty);
			_h["pT_y1_EWK_clip_inf"]->fill(my_pty);
			_h["pT_l1_EWK"]->fill(lep[0].pT() / GeV);
			_h["pT_ll_EWK"]->fill((lep[0].mom() + lep[1].mom()).pT() / GeV);
			_h["pT_lly_EWK"]->fill(my_ptZy);
			_h["pT_lly_EWK_clip_inf"]->fill(my_ptZy);
			_h["Centrality_EWK"]->fill(Centrality);
			_h["Dyjj_EWK"]->fill(Dyjj);
			_h["mlly_EWK"]->fill(mlly);
			_h["Dphi_lly_jj_EWK"]->fill(my_dphi);
			_h["Dphi_lly_jj_EWK_clip_inf"]->fill(my_dphi);


      // do clipping - sometimes there are two Z and one gamma - in this case to avoid much work take Z with highest pt and gamma
      std::vector<FourMomentum> hs_bosons_z = {};
      std::vector<FourMomentum> hs_bosons_y = {};
      for(const Particle& p_rivet : event.allParticles()){
        ConstGenParticlePtr p_hepmc = p_rivet.genParticle();
        int status = p_hepmc->status();
        if (abs(status)==23 or abs(status)==22){
          int i_pid = p_hepmc->pid();
          FourMomentum i_mom = p_hepmc->momentum();
          if (abs(i_pid) == 23){hs_bosons_z.push_back(i_mom);}
          else if (abs(i_pid) == 22){hs_bosons_y.push_back(i_mom);}
        }
      }
      std::sort(hs_bosons_z.begin(), hs_bosons_z.end(), [](FourMomentum const &a, FourMomentum const &b) {return a.pT() > b.pT(); }); // biggest pT will be first in array
      std::sort(hs_bosons_y.begin(), hs_bosons_y.end(), [](FourMomentum const &a, FourMomentum const &b) {return a.pT() > b.pT(); });
      bool have_two_hs_bosons = false;
      double hs_diboson_mass = 0.0;
      if (hs_bosons_z.size()>0 && hs_bosons_y.size()>0){
        hs_diboson_mass = (hs_bosons_z[0]+hs_bosons_y[0]).mass()/GeV;
        have_two_hs_bosons = true;
      }
      if (!have_two_hs_bosons) vetoEvent; // just in case reject events where dont have z+y

      _h["m_Zy"]->fill(hs_diboson_mass);
      if (ev_nominal_weight >= 0) { _c["pos_w_final_clip_inf"]->fill(); }
      else { _c["neg_w_final_clip_inf"]->fill(); }
      // fill clipped
      for (std::string &i_clip: _clips) {
        if (i_clip == "inf") continue; // as done above without cut
        int i_clip_num = std::stoi(i_clip);
        std::string i_pos_c_name = "pos_w_final_clip_" + i_clip;
        std::string i_neg_c_name = "neg_w_final_clip_" + i_clip;
        std::string i_hist_name_dphi = "Dphi_lly_jj_EWK_clip_" + i_clip;
        std::string i_hist_name_pty = "pT_y1_EWK_clip_" + i_clip;
        std::string i_hist_name_ptZy = "pT_lly_EWK_clip_" + i_clip;
        if (hs_diboson_mass < i_clip_num) {
          _h[i_hist_name_dphi]->fill(my_dphi);
          _h[i_hist_name_pty]->fill(my_pty);
          _h[i_hist_name_ptZy]->fill(my_ptZy);
          if (ev_nominal_weight >= 0) { _c[i_pos_c_name]->fill(); }
          else { _c[i_neg_c_name]->fill(); }
        }
      }

    } // end of analyze

    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSection() / femtobarn / sumOfWeights();
      scale(_h, sf);

      // normalize all to 1 since in case of mostly negative weights not clear what it will do
      const std::vector<std::string> h_names = {"Dphi_lly_jj_EWK_clip_", "pT_y1_EWK_clip_", "pT_lly_EWK_clip_"};
      for (auto & i_name : h_names){
        for (std::string &i_clip: _clips) {
        std::string  i_clip_name = i_name + i_clip;
        std::cout << "normalizeing hist " << i_clip_name <<" to 1; " ;
        normalize(_h[i_clip_name], 1.0);
        }
      }

    } // end of finalize

  private:

    /// Histograms
    map<string, Histo1DPtr> _h;
    map <string, CounterPtr> _c;
    std::string _cut_mode;
    std::vector<std::string> _clips{"inf", "3000", "2000", "1500", "1000", "700"};

    // Data members like post-cuts event weight counters go here
    size_t _mode;
  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2023_I2663725);

}
