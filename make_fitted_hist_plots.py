import lib_utils as lu
import os
import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--tProd", default = "Zy")
parser.add_option("--tDec", default = "vvy")
parser.add_option("--runWithCuts", default = "yes")
opts, _ = parser.parse_args()

prod_dec = f"{opts.tProd}_{opts.tDec}"
_, top_files_dir = lu.find_prod_dec_and_dir(f"user.okurdysh.MadGraph_{prod_dec}_FM0_SM")
docut_dir = "DOCUT_YES" if opts.runWithCuts=="yes" else "DOCUT_NO"
plot_dir = lu.get_plotdir(prod_dec, docut_dir)

clips = ["","3000","2000", "1500"]
def get_clip_hist_name(orig_name, clip):
    clip_name = orig_name+f"_clip{clip}" if len(clip)>0  else orig_name 
    return clip_name
fit_plot_str, fit_plot_bins = lu.get_fitted_plot(prod_dec)
plots_to_save = [get_clip_hist_name(fit_plot_str, i_clip) for i_clip in clips]
print("plots to save", plots_to_save)

for op_dir in [i_obj for i_obj in os.listdir(top_files_dir) if os.path.isdir(top_files_dir + "/" + i_obj)]:
    if not any(["FT8_QUAD" in op_dir,"FT9_QUAD" in op_dir, "FT5_QUAD" in op_dir, "FT0_QUAD" in op_dir]): continue
    full_op_dir = os.path.join(top_files_dir,op_dir,docut_dir,"")
    op, order = lu.get_op_from_dir(op_dir, prod_dec)
    hists_file = full_op_dir + "hists.root"
    op_hists_default_bin = lu.read_hists(hists_file, plots_to_save)
    dressed_hists_file = ROOT.TFile(plot_dir+"/hists_for_ws_comparison.root", "UPDATE") # cannot pull it outside of outer loop somehow
    for i_clip in clips:
        # print("working on clip", i_clip)
        i_clip_hist_name = get_clip_hist_name(fit_plot_str, i_clip)
        i_clip_hist = op_hists_default_bin[i_clip_hist_name]
        i_fid_xsec_file = full_op_dir+f"xsec_times_frac_fb_clip{i_clip}.txt" if len(i_clip)>0 else full_op_dir+f"xsec_times_frac_fb.txt" 
        with open(i_fid_xsec_file, 'r') as f: i_fid_xsec_fb = float(f.read())
        i_clip_hist_dressed = lu.dress_hist(i_clip_hist, f"{op[0]}_{order}_{i_clip_hist_name}", 1, i_fid_xsec_fb*139, re_bins=fit_plot_bins)
        i_clip_hist_dressed.Write("", ROOT.TObject.kOverwrite)
    dressed_hists_file.Close()