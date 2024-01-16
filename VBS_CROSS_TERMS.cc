// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class VBS_CROSS_TERMS : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(VBS_CROSS_TERMS);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 4.5);

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs, "jets");
      
      // FinalState of direct photons and bare muons and electrons in the event
      DirectFinalState photons(Cuts::abspid == PID::PHOTON);
      DirectFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);
      // Dress the bare direct leptons with direct photons within dR < 0.1,
      // and apply some fiducial cuts on the dressed leptons
      Cut lepton_cuts = Cuts::abseta < 2.5 && Cuts::pT > 20*GeV;
      DressedLeptons dressed_leps(photons, bare_leps, 0.1, lepton_cuts);
      declare(dressed_leps, "leptons");

      // Book histograms
      // specify custom binning
      book(_h["njet"], "njet", 10, 0.0, 10.0);
      book(_h["pt_jet1"], "pt_jet1", 200, 0.0, 3000.0);
      book(_h["mjj"], "mjj", 200, 0.0, 3000.0);
      book(_h["dyjj"], "dyjj", 50, 0.0, 6.0);
      book(_h2["mjj_dyjj"], "mjj_dyjj", 200, 0.0, 3000.0, 50, 0.0, 6.0);
      book(_c["found_VBS_pair"],"found_VBS_pair");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Retrieve dressed leptons, sorted by pT
      Particles leptons = apply<FinalState>(event, "leptons").particles();

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 30*GeV);
      // Remove all jets within dR < 0.2 of a dressed lepton
      idiscardIfAnyDeltaRLess(jets, leptons, 0.2);

      int njets = jets.size();
      if (njets < 2)  vetoEvent;  

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

      if (!foundVBSJetPair)  vetoEvent;

      FourMomentum tag1_jet_4vec = jets[tag1_jet_index].mom();
      FourMomentum tag2_jet_4vec = jets[tag2_jet_index].mom();
      const double mjj = (tag1_jet_4vec + tag2_jet_4vec).mass()/GeV;
      const double dyjj = fabs(tag1_jet_4vec.rap() - tag2_jet_4vec.rap());

      _h["njet"]->fill(jets.size());
      _h["pt_jet1"]->fill(jets[0].pt());
      _h["mjj"]->fill(mjj);
      _h["dyjj"]->fill(dyjj);
      _h2["mjj_dyjj"]->fill(mjj,dyjj);
      _c["found_VBS_pair"]->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double veto_survive_sumW = dbl(*_c["found_VBS_pair"]);
      double veto_survive_frac = veto_survive_sumW / sumW();
      std::cout << "survived veto, will norm to this: " << veto_survive_frac << "\n";
      double norm_to = veto_survive_frac*crossSection()/picobarn; // norm to generated cross-section in pb (after cuts)
      normalize(_h["njet"], norm_to);
      normalize(_h["pt_jet1"], norm_to);
      normalize(_h["mjj"], norm_to); 
      normalize(_h["dyjj"], norm_to);
      normalize(_h2["mjj_dyjj"], norm_to);

    }

    /// @}


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, Histo2DPtr> _h2;
    // map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(VBS_CROSS_TERMS);

}
