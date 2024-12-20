import lib_utils as lu
import os
import json
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
parser.add_option("--tGenProd", default = "Zy")
parser.add_option("--tGenDec", default = "vvy")
parser.add_option("--routine", default = "Zy_vvy")
parser.add_option("--cut", default = "SR")
parser.add_option("--searchClips", default = "inf,3000,2000")
parser.add_option("--additionalClips", default = "1500,1000,700") # add 700 even if know it's wrong
parser.add_option("--doReshuffling", default = 1, type="int")
parser.add_option("--doReplacement", default = 0, type="int")
parser.add_option("--doINT", default = 0, type="int")
parser.add_option("--doCROSS", default = 0, type="int")
parser.add_option("--savePdf", default=1, type="int")
opts, _ = parser.parse_args()

prod_dec = f"{opts.tGenProd}_{opts.tGenDec}"
_, top_files_dir = lu.find_prod_dec_and_dir(f"user.okurdysh.MadGraph_{prod_dec}_FM0_SM")
base_plot_dir, routine_dir = lu.get_plotdir(prod_dec, opts.routine, opts.cut)

search_clips = opts.searchClips.split(",")  # search based on this, avoid too low for rep search as they are super low stat
all_clips = search_clips + opts.additionalClips.split(",") # to check later what search gives
fit_plot_str, fit_plot_bins, _, _ = lu.get_fitted_plot(opts.routine,opts.cut)
xlabel = lu.get_var_latex(fit_plot_str)
def plotname(i_clip):
    return f"{fit_plot_str}_clip_{i_clip}"

# organize hists and xsec into dicts (save all search+ check)
hists = {}
for i_clip in all_clips: # separate seach per clipping
    for op_dir in [i_obj for i_obj in os.listdir(top_files_dir) if os.path.isdir(top_files_dir + "/" + i_obj)]:
        full_op_dir = os.path.join(top_files_dir,op_dir,routine_dir,"")
        _, order, op = lu.get_op_from_dir(op_dir, prod_dec)
        dictkey = (i_clip, order, op)
        # get check that hists and fid xsec file are ok
        hists_file = full_op_dir + "hists_norm_run2.root"
        op_hists_default_bin = lu.read_hists(hists_file, [plotname(i_clip)])
        if len(op_hists_default_bin)==0: 
            print("skip as didn't find hist", dictkey)
            continue
        # get hist
        i_hist_in, i_hist_name = op_hists_default_bin[plotname(i_clip)], plotname(i_clip)
        hists[dictkey] = i_hist_in #i_hist_dressed

print("collected hists and fid xsecs  and eff with uncert for", hists.keys())

def save_pairs_find_rep(arr_ops_1, arr_ops_2, order1, order2, make_rep_df=False):
    chi2_dfs = {} # per search clipping for all pairs
    for i_clip in all_clips:
        rchi2_arr = []
        pairs_arr = []    
        print("saveing paired hists and maybe doing search for clip order", i_clip, order2)
        pairs_dir = base_plot_dir + f"/pair_ratios_{order2}_clip_{i_clip}/"
        if not os.path.exists(pairs_dir): os.makedirs(pairs_dir)
        # get compatibility between all pairs per clipping and save plots
        for i_key1 in arr_ops_1:
            for i_key2 in arr_ops_2:
                if i_key1==i_key2: continue # "AvsB, BvsA both will be saved intentionally for easier latex batch plotting
                # with same normalization derive test
                i_hist1_raw = hists[(i_clip,order1,i_key1)]
                i_hist2_raw = hists[(i_clip,order2,i_key2)]
                i_hist1 = lu.dress_hist(i_hist1_raw, f"{i_key1}_{order1}_{i_clip}", 1) 
                i_hist2 = lu.dress_hist(i_hist2_raw, f"{i_key2}_{order2}_{i_clip}", 2) 
                ratio_plot, rchi2, _ = lu.get_ratio_plot_tests(i_hist1, i_hist2)
                i_stack = ROOT.THStack()
                i_stack.Add(i_hist1)
                i_stack.Add(i_hist2)
                lu.draw_stack_with_ratio(i_stack, ratio_plot, xlabel, pairs_dir+f"/{i_key1}vs{i_key2}.pdf")
                if make_rep_df:
                    pairs_arr.append([i_key1, i_key2])
                    rchi2_arr.append(rchi2)
        # for each clip
        # only for testted clip - save table with tests for all pairs
        if make_rep_df:
            df = pd.DataFrame(index=arr_ops_1, columns=arr_ops_2)
            for i_pair, i_rchi2 in zip(pairs_arr, rchi2_arr):
                i_op1, i_op2 = i_pair[0], i_pair[1]
                df.at[i_op1, i_op2] = round(i_rchi2,2)
            lu.save_df(df, f"{base_plot_dir}/{order2}_tests_table_clip_{i_clip}.pdf",
                       save_csv=True, save_pdf=opts.savePdf, save_latex=1)
            chi2_dfs[i_clip] = df
    if make_rep_df:
        sum_chi2_df = chi2_dfs[search_clips[0]].add(chi2_dfs[search_clips[1]])
        for i_clip_add in search_clips[2:]:
            sum_chi2_df = sum_chi2_df.add(chi2_dfs[i_clip_add])
        # sum_chi2_df.astype('float64').round(2)
        lu.save_df(sum_chi2_df, f"{base_plot_dir}/{order2}_tests_table_clip_sum.pdf",
                   save_csv=True, save_pdf=opts.savePdf, save_latex=1)
        return sum_chi2_df
    else:
        return 0

# find replacements
have_ops_int_quad = sorted(list(set([ihist[2] for ihist in hists.keys() if ihist[1]!="CROSS"])))
sum_chi2_df_quad = save_pairs_find_rep(have_ops_int_quad, have_ops_int_quad, "QUAD", "QUAD", make_rep_df=True)
if opts.doINT:
    save_pairs_find_rep(have_ops_int_quad, have_ops_int_quad, "INT", "INT") # dont search based on INT, apply what was found in QUAD selected terms - just save plots
if opts.doCROSS:
    have_ops_cross = sorted(list(set([ihist[2] for ihist in hists.keys() if ihist[1]=="CROSS"]))) # ["FM0vsFM1", "FT0vsFT5", "FT8vsFT9"]
    sum_chi2_df_cross = save_pairs_find_rep(have_ops_int_quad, have_ops_cross, "QUAD", "CROSS", make_rep_df=True)
    missing_ops_cross = list(sum_chi2_df_cross.columns)

# save table for reshuffling and replacement and also of renormalizations
missing_ops_ana = lu.get_missing_ops(opts.routine)
missing_ops_quad = [i_op for i_op in missing_ops_ana if i_op in list(sum_chi2_df_quad.keys())] # like can miss because irrelevant for process
existing_ops_quad = [i_op for i_op in list(sum_chi2_df_quad.keys()) if i_op not in missing_ops_quad]
print("have existing ops", existing_ops_quad)
print("have missing ops", missing_ops_quad)
cols_rep_df = ["toreplace", "replacement","sumchi2"]
def make_df(df_to_search, to_be_replaced_ops, search_within_ops, order_to_replace, order_rep, out_name):
    rep_df = pd.DataFrame(columns=cols_rep_df)
    for i_op_to_rep in to_be_replaced_ops:
        rep_op, resh_chi2 = lu.get_replacement(df_to_search, i_op_to_rep, search_within_ops)
        print("for op", i_op_to_rep, "replacement op is", rep_op, "with sum of chi2 across clips", resh_chi2)
        rep_df.loc[len(rep_df.index)] = [i_op_to_rep, rep_op, resh_chi2] 
    lu.save_df(rep_df, out_name, save_csv=True, aspect=(6,6), save_pdf=opts.savePdf, save_latex=1)
    # find renormalizations 
    df_norm = pd.DataFrame(index=list(rep_df["toreplace"]),columns=["rep"]+all_clips)
    df_norm_unc = pd.DataFrame(index=list(rep_df["toreplace"]),columns=["rep"]+all_clips)
    for _, rep_row in rep_df.iterrows():
        to_replace = rep_row["toreplace"]
        rep = rep_row["replacement"]
        for i_clip in all_clips:
            xsec_to_replace = hists[(i_clip,order_to_replace,to_replace)].Integral() #fid_xsecs[(i_clip, order_to_replace, to_replace)]
            xsec_rep = hists[(i_clip,order_rep,rep)].Integral() #fid_xsecs[(i_clip, order_rep, rep)]
            norm = xsec_to_replace/xsec_rep
            df_norm.at[to_replace, i_clip] = norm
            df_norm.at[to_replace, "rep"] = rep
            
    lu.save_df(df_norm, out_name+f"_norms_{order_to_replace}.pdf", save_csv=True, aspect=(6,6),
               save_pdf=opts.savePdf, save_latex=1)
    return rep_df, df_norm

if opts.doReshuffling:
    print("##### reshuffling based on QUAD - find QUAD coeficienes and INT?", opts.doINT) # find good reps for non-closure
    df_resh_q, df_norms_resh_q = make_df(sum_chi2_df_quad, existing_ops_quad, existing_ops_quad, "QUAD", "QUAD", f"{base_plot_dir}/reshuffling_ws_table.pdf")
    if opts.doINT:
        df_resh_i, df_norms_resh_i = make_df(sum_chi2_df_quad, existing_ops_quad, existing_ops_quad, "INT", "INT", f"{base_plot_dir}/reshuffling_ws_table.pdf")
if opts.doReplacement:
    print("##### replace missing QUAD - find  QUAD coeficienes and INT?", opts.doINT) # find good reps for missing
    df_rep_q, df_norms_rep_q = make_df(sum_chi2_df_quad, missing_ops_quad, existing_ops_quad, "QUAD", "QUAD", f"{base_plot_dir}/replacement_ws_table.pdf")
    if opts.doINT:
        df_rep_i, df_norms_rep_i = make_df(sum_chi2_df_quad, missing_ops_quad, existing_ops_quad, "INT", "INT", f"{base_plot_dir}/replacement_ws_table.pdf")
if opts.doCROSS:
    print("##### replace missing CROSS") # find way to insert missing crosses
    df_rep_c, df_norms_rep_c = make_df(sum_chi2_df_cross, missing_ops_cross, existing_ops_quad, "CROSS", "QUAD", f"{base_plot_dir}/cross_ws_table.pdf")

# convert to counts for mathieu
if opts.doCROSS and opts.routine in ["ssWW_lvlv","ATLAS_2023_I2729396"]:
    clips_for_counts = ["inf","1500"] # currently only have this json from Mathieu
    gen_prod_dec = f"{opts.tGenProd}_{opts.tGenDec}"
    base_plot_dir, _ = lu.get_plotdir(gen_prod_dec, opts.routine, opts.cut)
    workspaces = {}
    for clip in clips_for_counts:
        with open(f'{os.path.dirname(base_plot_dir[:-1])}/ws_{clip}.json') as f:
            ws = json.load(f)
            workspaces[clip] = ws
    def get_original_ws_array(clip,op, order_capital):
        assert order_capital!="CROSS", "man you're doing something wrong machinery not adapted to get CROSS"
        clip_h_name = lu.get_clip_hist_name(fit_plot_str, clip)
        _, arr_counts = lu.ssWW_get_ws_hist(workspaces[clip], op, order_capital,
                                            f"{op}_{order}_{clip_h_name}", gen_prod_dec, opts.cut)
        return arr_counts

    def save_counts(clip, df_with_norms, order_rep,order_to_rep,out_file):
        f = open(out_file, "w")
        for to_replace, rep_row in df_with_norms.iterrows():
            rep = rep_row["rep"]
            rep_arr_no_coef = get_original_ws_array(clip,rep, order_rep)
            # print("arr from ws for QUAD", clip,rep, "is", rep_arr_no_coef)
            renorm = rep_row[clip]
            rep_arr_scaled = np.array(rep_arr_no_coef) * renorm
            # print("which is scaled by", renorm, "is", rep_arr_scaled, "and this is now", to_replace)
            f.write(f"{to_replace}_{order_to_rep}:{','.join([str(count) for count in rep_arr_scaled])}\n")
        f.close()

    for clip in clips_for_counts:
        if opts.doCROSS:
            save_counts(clip, df_norms_rep_c,
                        "QUAD", "CROSS",f"{base_plot_dir}/counts_CROSS_clip_{clip}.txt")