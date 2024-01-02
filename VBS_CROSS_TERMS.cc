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

      if (jets.size() < 2)  vetoEvent;  

      FourMomentum jet_lead = jets[0].mom(); // which assume to be 1st tagging jet
      FourMomentum jet_tag2;
      bool foundVBSJetPair = false; // look in opposite hemispheres
      for (const Jet& jet : jets) {
        if(jet.eta()*jets[0].eta() < 0.) {
          jet_tag2 = jet.mom();
          foundVBSJetPair = true;
          break;
        }
      }
      if (!foundVBSJetPair)  vetoEvent;
     
      const double mjj = (jet_lead + jet_tag2).mass()/GeV;
      const double dyjj = fabs(jet_lead.rap() - jet_tag2.rap());

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
      double norm_to = veto_survive_frac*crossSection()/picobarn; // norm to generated cross-section in pb (after cuts)
      normalize(_h["njet"], norm_to);
      normalize(_h["pt_jet1"], norm_to);
      normalize(_h["mjj"], norm_to); 
      normalize(_h["dyjj"], norm_to);
      normalize(_h2["mjj_dyjj"], norm_to);
      std::cout << "xsec incoming in pb"<<  crossSection()/picobarn <<"frac of w surviving vetos " << veto_survive_frac  <<" integral mjj after scaling (one for syst variation) " << _h["mjj"]->integral() <<"\n";

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
