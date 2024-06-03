import lib_utils as lu
import os
import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--gen_prod_dec", default = "Zy_vvy")
parser.add_option("--routine", default = "Zy_vvy")
opts, _ = parser.parse_args()
cut = "SR"

_, top_files_dir = lu.find_prod_dec_and_dir(f"user.okurdysh.MadGraph_{opts.gen_prod_dec}_FM0_SM")
base_plot_dir, routine_dir = lu.get_plotdir(opts.gen_prod_dec, opts.routine, cut)

# names are like  M1_quad_pt_photon_clip_1000 or T2_quad_pt_lepton_clip_3000
clips = ["inf","3000","2000","1500","1000","700"] if opts.gen_prod_dec=="Zy_vvy" else ["inf","3000","2000","1500","1000"]

fit_plot_str, fit_plot_bins = lu.get_fitted_plot(opts.routine, cut)
plots_to_save = [lu.get_clip_hist_name(fit_plot_str, i_clip) for i_clip in clips]
print("plots to save", plots_to_save)

ws_hists_name = "Zvvy" if opts.gen_prod_dec=="Zy_vvy" else "Wy"
ws_hist_file = ROOT.TFile(f"../../fits/ws_extracted_hists/{ws_hists_name}.root", "read")
ws_hist_list = [ih.GetName() for ih in list(ws_hist_file.GetListOfKeys())]
for op_dir in [i_obj for i_obj in os.listdir(top_files_dir) if os.path.isdir(top_files_dir + "/" + i_obj)]:
    full_op_dir = os.path.join(top_files_dir,op_dir,routine_dir,"")
    op, order, _ = lu.get_op_from_dir(op_dir, opts.gen_prod_dec)
    if order=="CROSS":
        continue
    hists_file = full_op_dir + "hists.root"
    op_hists_default_bin = lu.read_hists(hists_file, plots_to_save)
    for i_clip in clips:
        i_plot_dir = f"{base_plot_dir}/ws_rivet_hists/clip_{i_clip}/"
        if not os.path.exists(i_plot_dir): os.makedirs(i_plot_dir)
        # print("working on clip", i_clip)
        i_clip_hist_name = lu.get_clip_hist_name(fit_plot_str, i_clip)
        if i_clip_hist_name not in op_hists_default_bin.keys():
            continue 
        i_clip_hist = op_hists_default_bin[i_clip_hist_name]
        i_fid_xsec_file = full_op_dir+f"xsec_times_frac_fb_clip_{i_clip}.txt" 
        with open(i_fid_xsec_file, 'r') as f: i_fid_xsec_fb = float(f.read())
        i_full_hist_name = f"{op[0]}_{order}_{i_clip_hist_name}"
        i_clip_hist_dressed = lu.dress_hist(i_clip_hist, "rivet_"+i_full_hist_name, 2, i_fid_xsec_fb*139 / (2.4**2 * 1.55),
                                            # re_bins=fit_plot_bins, 
                                            # re_overflow=0
                                            )
        # print(i_clip_hist_dressed)
        #
        rivet_h_name = i_full_hist_name.replace("QUAD", "quad").replace("INT", "lin")[1:] # 1: for FT1->T1
        if rivet_h_name not in ws_hist_list:
            continue # crosses are not there and some ops like T1 are not
        print("getting name from ws", rivet_h_name)
        rivet_hist = ws_hist_file.Get(rivet_h_name)
        print(rivet_hist)
        #
        stack = ROOT.THStack()
        stack.Add(i_clip_hist_dressed)
        stack.Add(rivet_hist)
        c = ROOT.TCanvas()
        stack.Draw("nostack")
        c.BuildLegend()
        # ROOT.gPad.SetLogy()
        c.Modified()
        c.Update()
        c.Show()
        c.SaveAs(i_plot_dir + i_full_hist_name + ".pdf")
ws_hist_file.Close()