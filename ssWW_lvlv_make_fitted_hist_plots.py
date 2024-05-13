import lib_utils as lu
import os
import ROOT
import pandas as pd
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--runWithCuts", default = "yes")
opts, _ = parser.parse_args()

prod_dec = "ssWW_lvlv"
_, top_files_dir = lu.find_prod_dec_and_dir(f"user.okurdysh.MadGraph_{prod_dec}_FM0_SM")
docut_dir = "DOCUT_YES" if opts.runWithCuts=="yes" else "DOCUT_NO"
base_plot_dir = lu.get_plotdir(prod_dec, docut_dir)

clips = ["inf","1500"]

fit_plot_str, fit_plot_bins = lu.get_fitted_plot(prod_dec)
def get_clip_hist_name(fitvar,clip):
    return f"{fitvar}_clip_{clip}"
plots_to_save = [get_clip_hist_name(fit_plot_str, i_clip) for i_clip in clips]
print("plots to save", plots_to_save)

excel = {}
for clip in clips:
    excel[clip] = pd.read_excel(base_plot_dir+f"ws_ww_yields_{clip}.ods", engine="odf", index_col=0)
ws_hists = {}
for i_clip,i_df in excel.items():
    for i_hist_name,i_row_bin in i_df.iterrows():
        arr_counts = list(i_row_bin)
        i_hist = lu.arr_to_hist1d(arr_counts, i_hist_name, fit_plot_bins)
        ws_hists[(i_clip,i_hist_name)] = i_hist
ws_hist_list = [i_key[1] for i_key in ws_hists.keys()]
ws_hist_ops_orders = [i_name[:i_name.find("_m_ll")] for i_name in ws_hist_list]
# final wilson = 8 * their wilson. table gives those numbers. theis is brough up for some reason
wilson_corr_coef = pd.read_excel(base_plot_dir+"ws_ww_wilson_correction.ods", engine="odf", index_col=0)
def get_wilson_coef(ops_arr):
    if order == "QUAD":
        mult_factor_on_my_xsec = wilson_corr_coef.at[ops_arr[0], "corr"] ** 2
    elif order == "INT":
        mult_factor_on_my_xsec = wilson_corr_coef.at[ops_arr[0], "corr"]
    else:  # cross
        mult_factor_on_my_xsec = wilson_corr_coef.at[ops_arr[0], "corr"] * wilson_corr_coef.at[ops_arr[1], "corr"]
    return mult_factor_on_my_xsec

for op_dir in [i_obj for i_obj in os.listdir(top_files_dir) if os.path.isdir(top_files_dir + "/" + i_obj)]:
    full_op_dir = os.path.join(top_files_dir,op_dir,docut_dir,"")
    ops_arr, order, op_str = lu.get_op_from_dir(op_dir, prod_dec)
    if f"{op_str}_{order}" not in ws_hist_ops_orders:
        continue
    mult_factor_on_my_xsec = get_wilson_coef(ops_arr)
    rivet_hists_file = full_op_dir + "hists.root"
    rivet_op_hists_default_bin = lu.read_hists(rivet_hists_file, plots_to_save)
    for i_clip in clips:
        i_plot_dir = f"{base_plot_dir}/ws_rivet_hists/clip_{i_clip}/"
        if not os.path.exists(i_plot_dir): os.makedirs(i_plot_dir)
        i_rivet_clip_hist_name = get_clip_hist_name(fit_plot_str, i_clip)
        i_rivet_clip_hist = rivet_op_hists_default_bin[i_rivet_clip_hist_name]
        i_rivet_fid_xsec_file = full_op_dir + f"xsec_times_frac_fb_clip_{i_clip}.txt"
        with open(i_rivet_fid_xsec_file, 'r') as f: i_rivet_fid_xsec_fb = float(f.read())
        h_name = f"{op_str}_{order}_{i_rivet_clip_hist_name}"
        i_rivet_norm = i_rivet_fid_xsec_fb * 139 * mult_factor_on_my_xsec
        i_rivet_clip_hist_dressed = lu.dress_hist(i_rivet_clip_hist, "rivet_" + h_name, 2, i_rivet_norm,
                                                  re_bins=fit_plot_bins, re_overflow=1)
        print("getting name from ws", h_name)
        ws_hist = ws_hists[(i_clip,h_name)]
        print(ws_hist)
        #
        stack = ROOT.THStack()
        stack.Add(i_rivet_clip_hist_dressed)
        stack.Add(ws_hist)
        c = ROOT.TCanvas()
        stack.Draw("nostack")
        c.BuildLegend()
        # ROOT.gPad.SetLogy()
        c.Modified()
        c.Update()
        c.Show()
        c.SaveAs(i_plot_dir + h_name + ".pdf")