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
#include "Rivet/Math/MathUtils.hh"

namespace Rivet {


    /// @brief Add a short analysis description here
    class ZZ_llll : public Analysis {
    public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ZZ_llll);


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
        const FinalState fs(Cuts::abseta < 4.5);

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
        if (_docut==1){lepton_cuts= Cuts::abseta < 2.5 && Cuts::abseta > 0.1 && Cuts::pT > 7.0*GeV;} 
        else{lepton_cuts= Cuts::abseta < 10.0 && Cuts::pT > 0.001*GeV;}
        
        DressedLeptons dressed_leps(photons, bare_leps, 0.1, lepton_cuts);
        declare(dressed_leps, "leptons_stable");

        declare(MissingMomentum(), "METFinder");

        // Book histograms
        int n_nbins = 10;
        int n_pt = 200;
        double max_pt = 3000.0;
        int n_rap = 60;
        double max_rap = 6.0;
        //jet plots
        book(_h["n_jet"], "n_jet", n_nbins, 0.0, n_nbins);
        book(_h["pt_tagjet1"], "pt_tagjet1", n_pt, 0.0, max_pt); 
        book(_h["pt_tagjet2"], "pt_tagjet2", n_pt, 0.0, max_pt);
        book(_h["m_tagjets"], "m_tagjets", n_pt, 0.0, max_pt);
        book(_h["dy_tagjets"], "dy_tagjets", n_rap, 0, max_rap);
        book(_h["eta_tagjets"], "eta_tagjets", n_rap, -1*max_rap, max_rap);
        book(_h["deta_tagjets"], "deta_tagjets", n_rap, 0, max_rap);
        book(_h["dphi_tagjets"], "dphi_tagjets", n_rap, 0, max_rap);
        // //lepton plots
        book(_h["leptons_sum_abs_pids"], "leptons_sum_abs_pids", n_nbins*6, 0, n_nbins*6);
        book(_h["n_lepton_stable"], "n_lepton_stable", n_nbins, 0.0, n_nbins);
        book(_h["lepton_pt"], "lepton_pt", int(n_pt), 0.0, max_pt);
        book(_h["lepton_eta"], "lepton_eta", n_rap, -1*max_rap, max_rap);
        book(_h["all_lep_pairs_m_ll"], "all_lep_pairs_m_ll", int(n_pt), 0.0, max_pt);
        book(_h["m_ll_of_pairs_best_quadruplet"], "m_ll_of_pairs_best_quadruplet", int(n_pt), 0.0, max_pt);
        book(_h["m_4l"], "m_4l", int(n_pt), 0.0, max_pt);
        //other
        book(_h["lep_pairs_dR"], "lep_pairs_dR", 10, 0.0, 1); 
        book(_c["found_VBS_pair"],"found_VBS_pair");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
        double m_z = 91.18; 
        // Retrieve dressed leptons, sorted by pT
        Particles leptons_stable = apply<FinalState>(event, "leptons_stable").particlesByPt();
        int nlep_stable = leptons_stable.size();
        if (nlep_stable!=4)  vetoEvent; // meaning all e,mu and not tau

        const Particle& lep1 = leptons_stable[0];
        const Particle& lep2 = leptons_stable[1]; 
        const Particle& lep3 = leptons_stable[2]; 
        const Particle& lep4 = leptons_stable[3];
        // std::cout << "pt of leptons expect to be high-low" <<  lep1.pT() <<" " << lep2.pT() << " " << lep3.pT() << " " <<lep4.pT() <<"\n"; 
        if (_docut==1 && (lep1.pT()<20.0 || lep2.pT()<20.0 || lep3.pT()<10.0)) vetoEvent; 
        double m_4l = (lep1.mom() + lep2.mom() + lep3.mom() + lep4.mom()).mass()/GeV;

        int sum_charges = lep1.charge() + lep2.charge() + lep3.charge() + lep4.charge();  
        if (sum_charges!=0) vetoEvent; // want two opposiet-sign pairs

        int sum_pids = lep1.pid() + lep2.pid() + lep3.pid() + lep4.pid();
        if (sum_pids!=0) vetoEvent; // want two same flavour pairs (anitimuon negative pid)
        int leptons_sum_abs_pids = fabs(lep1.pid()) + fabs(lep2.pid()) + fabs(lep3.pid()) + fabs(lep4.pid());
        

        // deltaR between leptons and m_ll against quarkonia
        for (int i = 0; i < nlep_stable; i++) {
            const Particle& i_lep = leptons_stable[i];
            for (int j = 0; j < nlep_stable; j++) {
                if (i!=j){
                    const Particle& j_lep = leptons_stable[j];
                    double i_dR = deltaR(i_lep.rap(),i_lep.phi(), j_lep.rap(),j_lep.phi());
                    if (_docut==1 && i_dR<0.2) vetoEvent;
                    if ((i_lep.pid() + j_lep.pid())==0){ 
                        double i_m_ll = (i_lep.mom() + j_lep.mom()).mass()/GeV;
                        if (_docut==1 && i_m_ll<10.0) vetoEvent;
                    }
                }
            }
        }
        // save to plot all m_ll pairs and deltaR
        for (int i = 0; i < nlep_stable; i++) {
            const Particle& i_lep = leptons_stable[i];
            for (int j = 0; j < nlep_stable; j++) {
                if (i!=j){
                    const Particle& j_lep = leptons_stable[j];
                    double i_dR = deltaR(i_lep.rap(),i_lep.phi(), j_lep.rap(),j_lep.phi());
                    _h["lep_pairs_dR"]->fill(i_dR);
                    if ((i_lep.pid() + j_lep.pid())==0){ 
                        double i_m_ll = (i_lep.mom() + j_lep.mom()).mass()/GeV;
                        if (_docut==1 && i_m_ll<10.0) vetoEvent;
                        _h["all_lep_pairs_m_ll"]->fill(i_m_ll);
                    }

                }
            }
        }
        
        ///////////
        //make quadruplet with closest to two Z peaks
        ///////////
        std::vector<int> quadruplet1_pair1_indexes = {-1, -1};
        for (int i = 0; i < nlep_stable; i++) {
            int i_pid = leptons_stable[i].pid();
            for (int j = 0; j < nlep_stable; j++) {
                if(i==j) continue;
                int j_pid = leptons_stable[j].pid();
                if(i_pid+j_pid!=0) continue; // skip if e-e-
                quadruplet1_pair1_indexes[0] = i;
                quadruplet1_pair1_indexes[1] = j;
                break;
                }
            }
        // std::cout << "found 1/2 first pair l+l- pair indexes" << quadruplet1_pair1_indexes[0] <<" " << quadruplet1_pair1_indexes[1] << "\n"; 
        if (quadruplet1_pair1_indexes[0]==-1 || quadruplet1_pair1_indexes[1]==-1) vetoEvent;
        // std::cout << "these particles have pdg" << leptons_stable[quadruplet1_pair1_indexes[0]].pid() <<" " << leptons_stable[quadruplet1_pair1_indexes[1]].pid() << "\n"; 
        std::vector<int> quadruplet1_pair2_indexes={};
        for (int i = 0; i < nlep_stable; i++) {
            if(i==quadruplet1_pair1_indexes[0] || i==quadruplet1_pair1_indexes[1]) continue;
            quadruplet1_pair2_indexes.push_back(i);
            } 
        // std::cout << "found 1/2 second pair l+l- pair indexes" << quadruplet1_pair2_indexes[0] <<" " << quadruplet1_pair2_indexes[1] << "\n";
        // std::cout << "these particles have pdg" << leptons_stable[quadruplet1_pair2_indexes[0]].pid() <<" " << leptons_stable[quadruplet1_pair2_indexes[1]].pid() << "\n"; 

        double q_1_m_pair1 = (leptons_stable[quadruplet1_pair1_indexes[0]].mom() + leptons_stable[quadruplet1_pair1_indexes[1]].mom()).mass()/GeV;
        double q_1_m_pair2 = (leptons_stable[quadruplet1_pair2_indexes[0]].mom() + leptons_stable[quadruplet1_pair2_indexes[1]].mom()).mass()/GeV;
        double q_1_Z_dist = fabs(q_1_m_pair1-m_z) + fabs(q_1_m_pair2-m_z); 
        // std::cout << "distance to Z in q1 "<< q_1_Z_dist << "\n";
        // find second multiplet, if eeee or mmmm
        std::vector<int> pair1_indexes={-1,-1};
        std::vector<int> pair2_indexes={-1,-1};
        if (leptons_sum_abs_pids!=48){ //for e+e mu+mu- there is only one quadruplet possible, here seach for two
            int index_neg2_lep=-1; //start from same particle  index_pos1_lep but now find second pair where it can go
            int pid_1 = leptons_stable[quadruplet1_pair1_indexes[0]].pid();
            for (int i = 0; i < nlep_stable; i++) {
                if(i==quadruplet1_pair1_indexes[0] || i==quadruplet1_pair1_indexes[1]) continue;
                int pid_i = leptons_stable[i].pid();
                if(pid_1+pid_i!=0) continue; // skip if l-l-, want opposite charges
                index_neg2_lep = i;
                }
            if (index_neg2_lep==-1) vetoEvent;
            std::vector<int> quadruplet2_pair1_indexes = {quadruplet1_pair1_indexes[0], index_neg2_lep};
            // std::cout << "found 2/2 first pair l+l- pair indexes" << quadruplet2_pair1_indexes[0] <<" " << quadruplet2_pair1_indexes[1] << "\n";
            // std::cout << "these particles have pdg" << leptons_stable[quadruplet2_pair1_indexes[0]].pid() <<" " << leptons_stable[quadruplet2_pair1_indexes[1]].pid() << "\n";
            std::vector<int> quadruplet2_pair2_indexes={};
            for (int i = 0; i < nlep_stable; i++) {
                if (i==quadruplet2_pair1_indexes[0] || i==quadruplet2_pair1_indexes[1]) continue;
                quadruplet2_pair2_indexes.push_back(i);
                }
            // std::cout << "found 2/2 second pair l+l- pair indexes" << quadruplet2_pair2_indexes[0] <<" " << quadruplet2_pair2_indexes[1] << "\n";
            // std::cout << "these particles have pdg" << leptons_stable[quadruplet2_pair2_indexes[0]].pid() <<" " << leptons_stable[quadruplet2_pair2_indexes[1]].pid() << "\n";
            double q_2_m_pair1 = (leptons_stable[quadruplet2_pair1_indexes[0]].mom() + leptons_stable[quadruplet2_pair1_indexes[1]].mom()).mass()/GeV;
            double q_2_m_pair2 = (leptons_stable[quadruplet2_pair2_indexes[0]].mom() + leptons_stable[quadruplet2_pair2_indexes[1]].mom()).mass()/GeV;
            double q_2_Z_dist = fabs(q_2_m_pair1-m_z) + fabs(q_2_m_pair2-m_z); 
            // std::cout << "distance to Z in q1 "<< q_1_Z_dist<< "and q2   "<< q_2_Z_dist << "\n"; 
            if (q_2_Z_dist < q_1_Z_dist){
                pair1_indexes[0] = quadruplet2_pair1_indexes[0];
                pair1_indexes[1] = quadruplet2_pair1_indexes[1];
                pair2_indexes[0] = quadruplet2_pair2_indexes[0];
                pair2_indexes[1] = quadruplet2_pair2_indexes[1];
            } 
            else {
                pair1_indexes[0] = quadruplet1_pair1_indexes[0];
                pair1_indexes[1] = quadruplet1_pair1_indexes[1];
                pair2_indexes[0] = quadruplet1_pair2_indexes[0];
                pair2_indexes[1] = quadruplet1_pair2_indexes[1];
            }
        }
        else{
            pair1_indexes[0] = quadruplet1_pair1_indexes[0];
            pair1_indexes[1] = quadruplet1_pair1_indexes[1];
            pair2_indexes[0] = quadruplet1_pair2_indexes[0];
            pair2_indexes[1] = quadruplet1_pair2_indexes[1];
        }
        // std::cout << "in the end indexes of first pair"<< pair1_indexes[0]<< " "<< pair1_indexes[1]<< "second pair"<< pair2_indexes[0]<< " "<< pair2_indexes[1] << "\n"; 
        ////////////
        ////// end of best quadruplet searching
        //////////////

        double m_pair1 = (leptons_stable[pair1_indexes[0]].mom() + leptons_stable[pair1_indexes[1]].mom()).mass()/GeV;
        double m_pair2 = (leptons_stable[pair2_indexes[0]].mom() + leptons_stable[pair2_indexes[1]].mom()).mass()/GeV;
        // std::cout << "m of best Z dist pair1 "<< m_pair1 << "of pair 2 " << m_pair2 << "\n";  
        if (_docut==1 && (m_pair1<60.0 || m_pair1>120.0)) vetoEvent; 
        if (_docut==1 && (m_pair2<60.0 || m_pair2>120.0)) vetoEvent; 

        // Retrieve clustered jets, sorted by pT, with a minimum pT cut
        Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 25*GeV);
        // Remove all jets within dR < 0.2 of a dressed lepton
        idiscardIfAnyDeltaRLess(jets, leptons_stable, 0.2);

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
        if (tag2_jet_index < tag1_jet_index) swap(tag1_jet_index, tag2_jet_index); // organize tag jets by pt  
        if (!foundVBSJetPair)  vetoEvent;

        const FourMomentum tag1_jet = jets[tag1_jet_index].mom();
        const FourMomentum tag2_jet = jets[tag2_jet_index].mom();
        if (_docut==1 && (tag1_jet.pT()<40.0 || tag2_jet.pT()<40.0)) vetoEvent; 

        const double m_tagjets = (tag1_jet + tag2_jet).mass()/GeV;
        if (_docut==1 && m_tagjets<300.0) vetoEvent;

        const double dy_tagjets = fabs(tag1_jet.rap() - tag2_jet.rap());
        if (_docut==1 && dy_tagjets<2.0) vetoEvent;
        
        //jet plots
        _h["n_jet"]->fill(njets);
        _h["pt_tagjet1"]->fill(tag1_jet.pt());
        _h["pt_tagjet2"]->fill(tag2_jet.pt());
        // above this worked before
        _h["m_tagjets"]->fill(m_tagjets);
        _h["dy_tagjets"]->fill(dy_tagjets);
        _h["eta_tagjets"]->fill(tag1_jet.eta()); _h["eta_tagjets"]->fill(tag2_jet.eta()); // fill both to the same hists
        _h["deta_tagjets"]->fill(deltaEta(tag1_jet,tag2_jet));
        _h["dphi_tagjets"]->fill(deltaPhi(tag1_jet,tag2_jet));
        //lepton plots
        _h["leptons_sum_abs_pids"]->fill(leptons_sum_abs_pids);
        _h["n_lepton_stable"]->fill(nlep_stable);
        _h["lepton_pt"]->fill(lep1.pT()); _h["lepton_pt"]->fill(lep2.pT());  _h["lepton_pt"]->fill(lep3.pT()); _h["lepton_pt"]->fill(lep4.pT());
        _h["lepton_eta"]->fill(lep1.eta()); _h["lepton_eta"]->fill(lep2.eta()); _h["lepton_eta"]->fill(lep3.eta()); _h["lepton_eta"]->fill(lep4.eta());
        // all_lep_pairs_m_ll filled in loop
        _h["m_ll_of_pairs_best_quadruplet"]->fill(m_pair1); _h["m_ll_of_pairs_best_quadruplet"]->fill(m_pair2);
        _h["m_4l"]->fill(m_4l);
        // other
        // lep_pairs_dR filled in loop
        _c["found_VBS_pair"]->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {

        double veto_survive_sumW = dbl(*_c["found_VBS_pair"]);
        double veto_survive_frac = veto_survive_sumW / sumW();
        std::cout << "survived veto, will norm to this: " << veto_survive_frac << "\n";
        double norm_to = veto_survive_frac*crossSection()/picobarn; // norm to generated cross-section in pb (after cuts)
        
        std::vector<std::string> hist_names_1d = {"n_jet","pt_tagjet1","pt_tagjet2","m_tagjets",
        // "dy_tagjets","eta_tagjet1","eta_tagjet2","eta_tagjets", "deta_tagjets", 
        // "phi_tagjet1","phi_tagjet2","dphi_tagjets",
        // "n_lepton_stable","lepton_pt","lepton_eta",
        // "m_ll","MET","m_T","centrality","jet3_centrality","tagjet1_index","tagjet2_index","jet3_index"
        };       
        for (auto&& i_name : hist_names_1d){ normalize(_h[i_name], norm_to);}
        // also norm few 2d
        // normalize(_h2["m_dy_tagjets"], norm_to);
        // normalize(_h2["leptons_pids"],norm_to);
        // normalize(_h2["tagjets_index"],norm_to);

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


    RIVET_DECLARE_PLUGIN(ZZ_llll);

}
