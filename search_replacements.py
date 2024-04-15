import lib_utils as lu
import os
import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gROOT.ProcessLine( "gErrorIgnoreLevel = 1001;") # dont print when every plots is saved
ROOT.gROOT.LoadMacro("~/atlasstyle/AtlasStyle.C")
ROOT.SetAtlasStyle()
ROOT.gROOT.LoadMacro("~/atlasstyle/AtlasLabels.C")
ROOT.gROOT.LoadMacro("~/atlasstyle/AtlasUtils.C")
ROOT.gStyle.SetLegendBorderSize(0)
ROOT.gStyle.SetLegendFillColor(0)
# ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetLegendTextSize(0.04)
import matplotlib.pyplot as plt
plt.rcParams.update({'text.usetex': True,
                     "backend": 'PDF'})
from matplotlib.backends.backend_pdf import PdfPages
import math
import numpy as np
import pandas as pd
from pandas.plotting import table
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--tProd", default = "Zy")
parser.add_option("--tDec", default = "vvy")
parser.add_option("--runWithCuts", default = "yes")
opts, _ = parser.parse_args()

prod_dec = f"{opts.tProd}_{opts.tDec}"
_, top_files_dir = lu.find_prod_dec_and_dir(f"user.okurdysh.MadGraph_{prod_dec}_FM0_SM")
docut_dir = "DOCUT_YES" if opts.runWithCuts=="yes" else "DOCUT_NO"
base_plot_dir = lu.get_plotdir(prod_dec, docut_dir)

clips = ["inf","3000","2000"]  # search based on this, avoid 1000,700 for rep search as they are super low stat
all_clips = clips + ["1500", "1000", "700"] # to check later what search gives
search_order = "QUAD"
all_orders = [search_order, "INT"]
fit_plot_str, fit_plot_bins = lu.get_fitted_plot(prod_dec)
xlabel = lu.get_var_latex(fit_plot_str)
def plotname(i_clip):
    return f"{fit_plot_str}_clip_{i_clip}"

# organize hists and xsec into dicts (save all search+ check)
hists = {}
fid_xsecs = {}
for i_clip in all_clips: # separate seach per clipping
    for op_dir in [i_obj for i_obj in os.listdir(top_files_dir) if os.path.isdir(top_files_dir + "/" + i_obj)]:
        full_op_dir = os.path.join(top_files_dir,op_dir,docut_dir,"")
        _, order, op = lu.get_op_from_dir(op_dir, prod_dec)
        if order not in all_orders: 
            continue
        dictkey = (i_clip, order, op)
        # get check that hists and fid xsec file are ok
        hists_file = full_op_dir + "hists.root"
        op_hists_default_bin = lu.read_hists(hists_file, [plotname(i_clip)])
        if len(op_hists_default_bin)==0: 
            print("skip as didn't find hist", dictkey)
            continue
        xsec_file = full_op_dir + f"xsec_times_frac_fb_clip_{i_clip}.txt"
        if not os.path.exists(xsec_file): 
            print("skio as didnt find xsec", dictkey)
            continue
        # get hist
        i_hist_in, i_hist_name = op_hists_default_bin[plotname(i_clip)], plotname(i_clip)
        i_hist_dressed = lu.dress_hist(i_hist_in, i_hist_name, 1, 1, re_bins=fit_plot_bins)
        hists[dictkey] = i_hist_dressed
        # get xsec
        with open(xsec_file, 'r') as f: 
            fid_xsec_fb = float(f.read())
        fid_xsecs[dictkey] = fid_xsec_fb
print("collected hists and fid xsecs for", hists.keys())

chi2_dfs = {} # per search clipping for all pairs
have_ops = sorted(list(set([ihist[2] for ihist in hists.keys()])))
for i_clip in all_clips:
    rchi2_arr = []
    pairs_arr = []    
    for order in all_orders:
        print("saveing paired hists and maybe doing search for clip", i_clip)
        pairs_dir = base_plot_dir + f"/pair_ratios_{order}_clip_{i_clip}/"
        if not os.path.exists(pairs_dir): os.makedirs(pairs_dir)
        # get compatibility between all pairs per clipping and save plots
        for i_key1 in have_ops:
            for i_key2 in have_ops:
                if i_key1==i_key2: continue # "AvsB, BvsA both will be saved intentionally for easier latex batch plotting
                # with same normalization derive test
                i_hist1_raw = hists[(i_clip,order,i_key1)]
                i_hist2_raw = hists[(i_clip,order,i_key2)]
                i_hist1 = lu.dress_hist(i_hist1_raw, f"{i_key1}_{order}_{i_clip}", 1) 
                i_hist2 = lu.dress_hist(i_hist2_raw, f"{i_key2}_{order}_{i_clip}", 2) 
                ratio_plot, rchi2, _ = lu.get_ratio_plot_tests(i_hist1, i_hist2)
                i_stack = ROOT.THStack()
                i_stack.Add(i_hist1)
                i_stack.Add(i_hist2)
                lu.draw_stack_with_ratio(i_stack, ratio_plot, xlabel, pairs_dir+f"/{i_key1}vs{i_key2}.pdf")
                if i_clip in clips and order==search_order: # only for search stuff save shape test res
                    pairs_arr.append([i_key1, i_key2])
                    rchi2_arr.append(rchi2)
    # for each clip
    # only for testted clip - save table with tests for all pairs
    df = pd.DataFrame(index=have_ops, columns=have_ops)
    for i_pair, i_rchi2 in zip(pairs_arr, rchi2_arr):
        i_op1, i_op2 = i_pair[0], i_pair[1]
        df.at[i_op1, i_op2] = round(i_rchi2,2)
    lu.save_df(df, f"{base_plot_dir}/{search_order}_tests_table_clip_{i_clip}.pdf")
    chi2_dfs[i_clip] = df
sum_chi2_df = chi2_dfs[clips[0]].add(chi2_dfs[clips[1]])
for i_clip_add in clips[2:]:
    sum_chi2_df = sum_chi2_df.add(chi2_dfs[i_clip_add])
sum_chi2_df.astype('float64').round(2)
lu.save_df(sum_chi2_df, f"{base_plot_dir}/{search_order}_tests_table_clip_sum.pdf")

# save table for reshuffling and replacement and also of renormalizations
missing_ops_ana = lu.get_missing_ops(prod_dec)
missing_ops = [i_op for i_op in missing_ops_ana if i_op in list(sum_chi2_df.keys())] # like can miss because irrelevant for process
existing_ops = [i_op for i_op in list(sum_chi2_df.keys()) if i_op not in missing_ops]
cols_rep_df = ["toreplace", "replacement","sumchi2"]
def make_df(to_be_replaced_ops, search_within_ops, out_name):
    rep_df = pd.DataFrame(columns=cols_rep_df)
    for i_op_to_rep in to_be_replaced_ops:
        rep_op, resh_chi2 = lu.get_replacement(sum_chi2_df, i_op_to_rep, existing_ops)
        print("for op", i_op_to_rep, "replacement op is", rep_op, "with sum of chi2 across clips", resh_chi2)
        rep_df.loc[len(rep_df.index)] = [i_op_to_rep, rep_op, resh_chi2] 
    lu.save_df(rep_df, out_name, save_csv=True, aspect=(6,6))
    # find renormalizations 
    df_norms = [pd.DataFrame(index=list(rep_df["toreplace"]),columns=["rep"]+all_clips) for order in all_orders] # one per order
    for _, rep_row in rep_df.iterrows():
        to_replace = rep_row["toreplace"]
        rep = rep_row["replacement"]
        for i_clip in all_clips:
            for order,df_norm in zip(all_orders, df_norms): 
                xsec_to_replace = fid_xsecs[(i_clip,order,to_replace)]
                xsec_rep = fid_xsecs[(i_clip,order,rep)]
                norm = xsec_to_replace/xsec_rep
                df_norm.at[to_replace, i_clip] = norm
                df_norm.at[to_replace, "rep"] = rep
    for df_norm, i_order in zip(df_norms,all_orders):
        lu.save_df(df_norm, out_name+f"_norms_{i_order}.pdf", save_csv=True, aspect=(6,6))
    return rep_df
print("##### reshuffling") # find good reps for non-closure
df_resh = make_df(existing_ops, existing_ops, f"{base_plot_dir}/reshuffling_ws_table.pdf")
print("##### replace missing") # find good reps for missing
df_rep = make_df(missing_ops, existing_ops, f"{base_plot_dir}/replacement_ws_table.pdf")