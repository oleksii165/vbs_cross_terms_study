import lib_utils as lu
import os
import ROOT
import json
import pandas as pd
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--cut", default = "SR")
parser.add_option("--gen_prod_dec", default = "ssWW_lvlv") #WZ_lllv
opts, _ = parser.parse_args()
assert opts.cut in ["SR", "WZCR", "LowmjjCR"]
substr_cut_rivet_to_ws_dict = {"SR": "CutSR", "LowmjjCR": "CutLowMjj", "WZCR": "CutCRWZ3Lep"}
overflow_dict = {"SR": 1, "LowmjjCR": 1, "WZCR": 1}
include_overflow = overflow_dict[opts.cut]

clips = ["inf","1500"]
routine = "ssWW_lvlv"
fit_plot_str, fit_plot_bins = lu.get_fitted_plot(routine, opts.cut)
def get_clip_hist_name(fitvar,clip):
    return f"{fitvar}_clip_{clip}"
plots_to_save = [get_clip_hist_name(fit_plot_str, i_clip) for i_clip in clips]
print("plots to save", plots_to_save)
def get_ws_hist(ws_json_clip,op,order_capital, i_hist_out_name):
    eft_substr = "EFTWZ_" if opts.gen_prod_dec=="WZ_lllv" else "EFT_"
    op_substr = op[1:] if "F" in op else op
    arr_counts = -1
    # to select certain region it's "cutSR"
    index_region_dict= {"SR": 3, "WZCR": 0, "LowmjjCR": 2}
    sr_samples = ws_json_clip['distributions'][index_region_dict[opts.cut]]['samples']
    for i_hist_json in sr_samples:
        i_name = i_hist_json['name'] # like 'EFTWZ_M1_28_quad_CutSR_all_overallSyst_x_StatUncert'
        if eft_substr not in i_name: continue
        if op_substr not in i_name: continue
        if order_capital.lower() not in i_name: continue
        if substr_cut_rivet_to_ws_dict[opts.cut] not in i_name: continue
        arr_counts = i_hist_json['data']['contents']
    i_hist = lu.arr_to_hist1d(arr_counts, i_hist_out_name, fit_plot_bins) if arr_counts!=-1 else -1
    return i_hist


_, top_files_dir = lu.find_prod_dec_and_dir(f"user.okurdysh.MadGraph_{opts.gen_prod_dec}_FM0_SM")
base_plot_dir, routine_dir = lu.get_plotdir(opts.gen_prod_dec, routine, opts.cut)
# final wilson = 8 * their wilson. table gives those numbers. theis is brough up for some reason
wilson_corr_coef = {}
wilson_corr_coef["ssWW_lvlv"] = {"FS02":8.0,"FS1":20.0,"FM0":6.0,"FM1":10.0,"FM7":13.0,"FT0":0.6,"FT1":0.3,"FT2":1.0}
wilson_corr_coef["WZ_lllv"] = {"FS02":82,"FS1":128,"FM0":27,"FM1":28,"FM7":30,"FT0":2.4,"FT1":1.6,"FT2":5.5}
def get_wilson_coef(ops_arr):
    mult_factor_on_my_xsec = -1
    op_search = ops_arr[0]
    if op_search in wilson_corr_coef[opts.gen_prod_dec].keys():
        coef = wilson_corr_coef[opts.gen_prod_dec][op_search]
        if order == "QUAD":
            mult_factor_on_my_xsec = coef ** 2
        elif order == "INT":
            mult_factor_on_my_xsec = coef
    return mult_factor_on_my_xsec

workspaces = {}
for clip in clips:
    with open(f'{os.path.dirname(base_plot_dir[:-1])}/ws_{clip}.json') as f:
        ws = json.load(f)
        workspaces[clip] = ws

for op_dir in [i_obj for i_obj in os.listdir(top_files_dir) if os.path.isdir(top_files_dir + "/" + i_obj)]:
    full_op_dir = os.path.join(top_files_dir,op_dir,routine_dir,"")
    ops_arr, order, op_str = lu.get_op_from_dir(op_dir, opts.gen_prod_dec)
    mult_factor_on_my_xsec = get_wilson_coef(ops_arr)
    rivet_hists_file = full_op_dir + "hists.root"
    rivet_op_hists_default_bin = lu.read_hists(rivet_hists_file, plots_to_save)
    if len(rivet_op_hists_default_bin.keys())==0: continue
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
                                                  re_bins=fit_plot_bins, re_overflow=include_overflow)
        print("getting name from ws", h_name, op_str, order, opts.gen_prod_dec, "clip", i_clip)
        ws_hist = get_ws_hist(workspaces[i_clip], op_str,order, h_name)
        if ws_hist==-1:
            print("not able to find in workspace hist for", opts.gen_prod_dec)
            continue
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
        print("hi")