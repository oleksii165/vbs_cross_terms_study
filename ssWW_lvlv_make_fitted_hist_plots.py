import lib_utils as lu
import os
import ROOT
import json
import pandas as pd
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)

cut = "SR"
clips = ["inf","1500"]
gen_prod_dec_array = ["ssWW_lvlv", "WZ_lllv"]
routine = "ssWW_lvlv"
fit_plot_str, fit_plot_bins = lu.get_fitted_plot(routine)
def get_clip_hist_name(fitvar,clip):
    return f"{fitvar}_clip_{clip}"
plots_to_save = [get_clip_hist_name(fit_plot_str, i_clip) for i_clip in clips]
print("plots to save", plots_to_save)
def get_ws_hist(ws_json_clip,op,order_capital,gen_prod_dec, i_hist_out_name):
    eft_substr = "EFTWZ_" if gen_prod_dec=="WZ_lllv" else "EFT_"
    op_substr = op[1:] if "F" in op else op
    arr_counts = -1
    # to select certain region it's "cutSR"
    sr_samples = ws_json_clip['distributions'][3]['samples']
    for i_hist_json in sr_samples:
        i_name = i_hist_json['name'] # like 'EFTWZ_M1_28_quad_CutSR_all_overallSyst_x_StatUncert'
        if eft_substr not in i_name: continue
        if op_substr not in i_name: continue
        if order_capital.lower() not in i_name: continue
        arr_counts = i_hist_json['data']['contents']
    i_hist = lu.arr_to_hist1d(arr_counts, i_hist_out_name, fit_plot_bins) if arr_counts!=-1 else -1
    return i_hist

for gen_prod_dec in gen_prod_dec_array:
    include_overflow = 1
    # include_overflow = 1 if gen_prod_dec=="ssWW_lvlv" else 0 # get giat peak in last bin for WZ otherwise?
    _, top_files_dir = lu.find_prod_dec_and_dir(f"user.okurdysh.MadGraph_{gen_prod_dec}_FM0_SM")
    base_plot_dir, routine_dir = lu.get_plotdir(gen_prod_dec, routine, cut)
    # final wilson = 8 * their wilson. table gives those numbers. theis is brough up for some reason
    wilson_corr_coef = pd.read_excel(base_plot_dir+"ws_ww_wilson_correction.ods", engine="odf", index_col=0)
    def get_wilson_coef(ops_arr):
        mult_factor_on_my_xsec = -1
        op_search = ops_arr[0]
        if op_search in list(wilson_corr_coef.index):
            if order == "QUAD":
                mult_factor_on_my_xsec = wilson_corr_coef.at[ops_arr[0], "corr"] ** 2
            elif order == "INT":
                mult_factor_on_my_xsec = wilson_corr_coef.at[ops_arr[0], "corr"]
        return mult_factor_on_my_xsec
    workspaces = {}
    for clip in clips:
        with open(base_plot_dir+f'ws_{clip}.json') as f:
            ws = json.load(f)
            workspaces[clip] = ws

    for op_dir in [i_obj for i_obj in os.listdir(top_files_dir) if os.path.isdir(top_files_dir + "/" + i_obj)]:
        full_op_dir = os.path.join(top_files_dir,op_dir,routine_dir,"")
        ops_arr, order, op_str = lu.get_op_from_dir(op_dir, gen_prod_dec)
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
                                                      re_bins=fit_plot_bins, re_overflow=include_overflow)
            print("getting name from ws", h_name, op_str, order, gen_prod_dec, "clip", i_clip)
            ws_hist = get_ws_hist(workspaces[i_clip], op_str,order, gen_prod_dec, h_name)
            if ws_hist==-1:
                print("not able to find in workspace hist for", gen_prod_dec, )
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