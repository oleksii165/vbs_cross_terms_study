import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gROOT.LoadMacro("~/atlasstyle/AtlasStyle.C")
ROOT.SetAtlasStyle()
ROOT.gROOT.LoadMacro("~/atlasstyle/AtlasLabels.C")
ROOT.gROOT.LoadMacro("~/atlasstyle/AtlasUtils.C")
ROOT.gStyle.SetLegendBorderSize(0)
ROOT.gStyle.SetLegendFillColor(0)
# ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetLegendTextSize(0.04)
import os
import lib_utils as lu
import matplotlib.pyplot as plt
from statistics import mean
plt.rcParams.update({'text.usetex': True,
                     "backend": 'PDF'})
from matplotlib.backends.backend_pdf import PdfPages
import math
import numpy as np
import pandas as pd
from pandas.plotting import table
import subprocess

from optparse import OptionParser
parser = OptionParser()
parser.add_option("--tProd", default = "WmWm")
parser.add_option("--tDec", default = "lvlv")
parser.add_option("--runSMAndFULL", default = "no")
parser.add_option("--runQUADAndCROSS", default = "no")
parser.add_option("--numJobsParallel", default = 15)
parser.add_option("--runWithCuts", default = "yes")
parser.add_option("--runAgain", default = "no")
parser.add_option("--SMOnSumPlots", default = "no")
parser.add_option("--drawCutflow", default = "no")
opts, _ = parser.parse_args()

prod_dec = f"{opts.tProd}_{opts.tDec}"
_, top_files_dir = lu.find_prod_dec_and_dir(f"user.okurdysh.MadGraph_{prod_dec}_FM0_SM")
docut_dir = "DOCUT_YES" if opts.runWithCuts=="yes" else "DOCUT_NO"

plot_dir = lu.get_plotdir(prod_dec, docut_dir)
big_pairs_cutoff = 1.05
ks_cut = 0.6
chi2_cut = 2.0
plots_no_logy=["n_","ph","dp","et","dy", "dR"]

plots_to_save_info_dict = lu.get_hists_bounds_cuts(prod_dec) 
plots_to_save = plots_to_save_info_dict.keys()
###########
# build summary plots out of files created above
################

plots_dict = {'SM':{}, 'FULL':{}, 'QUAD':{}, 'CROSS':{}, 'INT':{}}
for i_eft in plots_dict.keys():
    for i_plot in plots_to_save:
        plots_dict[i_eft][i_plot] = {}
xsec_dict = {'SM':{}, 'FULL':{}, 'QUAD':{}, 'CROSS':{}, 'INT':{}}
frac_dict = {'SM':{}, 'FULL':{}, 'QUAD':{}, 'CROSS':{}, 'INT':{}}
frac_er_dict = {'SM':{}, 'FULL':{}, 'QUAD':{}, 'CROSS':{}, 'INT':{}}
cutflow_file_dict = {'SM':{}, 'FULL':{}, 'QUAD':{}, 'CROSS':{}, 'INT':{}}
for op_dir in [i_obj for i_obj in os.listdir(top_files_dir) if os.path.isdir(top_files_dir + "/" + i_obj)]:
    full_op_dir = os.path.join(top_files_dir,op_dir,docut_dir,"")
    print("#### plotting new conf", op_dir, "full path", full_op_dir)
    _, log_file = lu.get_evnt_log_files(top_files_dir,op_dir[:op_dir.find("_EXT0")])
    ops_arr, regime = lu.get_op_from_dir(op_dir, prod_dec)
    print("this is for ops", ops_arr, "in regime", regime)
    product_file = full_op_dir + "xsec_times_frac_fb.txt"
    frac_file = full_op_dir + "frac_after_cuts.txt"
    frac_er_file = full_op_dir + "frac_after_cuts_error_bar.txt"
    hists_file = full_op_dir + "hists.root"
    cutflow_file = full_op_dir + "cutflow.txt"

    if not os.path.exists(product_file) or not os.path.exists(frac_file) or not os.path.exists(hists_file) or not os.path.exists(log_file):
        print("didnt find the txt xsec fid and/or root file or frac file or log file for", op_dir)
        continue 
    if not os.path.exists(cutflow_file) or not os.path.exists(frac_er_file):
        print("didnt find the txt cutflow  and/or frac error file", op_dir)
        continue

    # xsec organization
    with open(product_file, 'r') as f: fid_xsec_fb = float(f.read())
    print("reading back fid_xsec in fb", fid_xsec_fb)
    with open(frac_file, 'r') as f: frac = float(f.read())    
    print("reading back efficiency", frac)
    with open(frac_er_file, 'r') as f: frac_er = float(f.read())    
    print("reading back efficiency error", frac)
    if regime=="CROSS":
        my_op1, my_op2  = ops_arr[0], ops_arr[1]
        if my_op1 not in xsec_dict[regime].keys(): xsec_dict[regime][my_op1] = {}
        xsec_dict[regime][my_op1][my_op2] = fid_xsec_fb
        #
        if my_op1 not in frac_dict[regime].keys(): frac_dict[regime][my_op1] = {}
        frac_dict[regime][my_op1][my_op2] = frac
        #
        if my_op1 not in frac_er_dict[regime].keys(): frac_er_dict[regime][my_op1] = {}
        frac_er_dict[regime][my_op1][my_op2] = frac_er
        #
        if my_op1 not in cutflow_file_dict[regime].keys(): cutflow_file_dict[regime][my_op1] = {}
        cutflow_file_dict[regime][my_op1][my_op2] = cutflow_file
        
    else:
        my_op = ops_arr[0]
        xsec_dict[regime][my_op] = fid_xsec_fb
        frac_dict[regime][my_op] = frac
        frac_er_dict[regime][my_op] = frac_er
        cutflow_file_dict[regime][my_op] = cutflow_file
    # hists organization
    op_hists = lu.read_hists(hists_file, plots_to_save)
    for i_hist_name in plots_to_save:
        if i_hist_name not in op_hists.keys(): continue

        if regime=="CROSS":
            my_op1,my_op2 = ops_arr[0], ops_arr[1]
            if my_op1 not in plots_dict[regime][i_hist_name].keys(): 
                plots_dict[regime][i_hist_name][my_op1] = {}
            plots_dict[regime][i_hist_name][my_op1][my_op2] = op_hists[i_hist_name] 
        else:
            my_op = ops_arr[0]
            plots_dict[regime][i_hist_name][my_op] = op_hists[i_hist_name]
    
print("xsec SM", xsec_dict["SM"])
print("xsec FULL", xsec_dict["FULL"])
print("xsec QUAD", xsec_dict["QUAD"])
for i_key1 in xsec_dict["CROSS"].keys():
    print("xsec CROSS",i_key1,xsec_dict["CROSS"][i_key1])
#
print("saved hists SM", [[i_plot, plots_dict["SM"][i_plot].keys()] for i_plot in plots_to_save])
print("saved hists FULL", [[i_plot, plots_dict["FULL"][i_plot].keys()] for i_plot in plots_to_save])
print("saved hists QUAD", [[i_plot, plots_dict["QUAD"][i_plot].keys()] for i_plot in plots_to_save])
print("saved hists CROSS")
for i_plot in plots_to_save:
    print("for plot", i_plot)
    for i_key1 in plots_dict["CROSS"][i_plot].keys():
        print("have hists", plots_dict["CROSS"][i_plot][i_key1].keys())

all_ops =  sorted(list(xsec_dict["QUAD"].keys()))
print("all ops", all_ops)
nbins = len(all_ops)
if opts.runSMAndFULL=="yes" and len(xsec_dict["SM"])>=1 and len(xsec_dict["FULL"])>=1:
    SM_ref =  xsec_dict["SM"][list(xsec_dict["SM"].keys())[0]] 
    print("SM xsec in fb", SM_ref)
    lu.write_to_f(top_files_dir+"/SM_xsec_times_frac_fb.txt", SM_ref)
    ############# FULL and in comparsion with SM
    FULL_h = ROOT.TH2F("FULL_h","FULL_h", nbins,0,nbins,1,0,1)
    FULL_ratio_SM_h = ROOT.TH2F("FULL_ratio_SM_h","FULL_ratio_SM_h", nbins,0,nbins,1,0,1)
    for num_bin,i_op in enumerate(all_ops, start=1):
        FULL_h.GetXaxis().SetBinLabel(num_bin,i_op)
        FULL_ratio_SM_h.GetXaxis().SetBinLabel(num_bin,i_op)
        if i_op in xsec_dict["FULL"].keys(): 
            FULL_h.SetBinContent(num_bin, 1, xsec_dict["FULL"][i_op])
            FULL_ratio_SM_h.SetBinContent(num_bin, 1, xsec_dict["FULL"][i_op] / SM_ref)
    lu.save_plot(FULL_h,plot_dir + "FULL.pdf")
    lu.save_plot(FULL_ratio_SM_h,plot_dir + "FULL_ratio_SM_h.pdf")

if opts.runQUADAndCROSS=="yes":
    ##########
    # make tables and efficiency plot
    ###########
    QUAD_h = ROOT.TH2F("QUAD_h","QUAD_h", nbins,0,nbins,1,0,1)
    for num_bin,i_op in enumerate(all_ops, start=1):
        QUAD_h.GetXaxis().SetBinLabel(num_bin,i_op)
        if i_op in xsec_dict["QUAD"].keys(): QUAD_h.SetBinContent(num_bin, 1, xsec_dict["QUAD"][i_op])
    lu.save_plot(QUAD_h,plot_dir + "QUAD.pdf")

    CROSS_h = ROOT.TH2F("CROSS_h","CROSS_h", nbins,0,nbins,nbins,0,nbins)
    CROSS_geom_QUAD_h = ROOT.TH2F("CROSS_geom_QUAD_h","CROSS_geom_QUAD_h", nbins,0,nbins,nbins,0,nbins)
    CROSS_el_area_ratio_h = ROOT.TH2F("CROSS_el_area_ratio_h","CROSS_el_area_ratio_h", nbins,0,nbins,nbins,0,nbins)
    pairs_str_avalaible = []
    area_ratios = []
    fracs_quad1 = []
    fracs_quad2 = []
    fracs_cross = []
    fracs_ave = []
    fracs_envelope = []
    fracs_errors_quad1 = []
    fracs_errors_quad2 = []
    fracs_errors_cross = []
    fracs_errors_ave = []
    cutflow_plot_dir = plot_dir + "/cutflows_comp/"
    if not os.path.exists(cutflow_plot_dir): os.makedirs(cutflow_plot_dir)
    for i_op1 in all_ops:
        if i_op1 not in xsec_dict["CROSS"].keys(): continue

        bin_x = all_ops.index(i_op1) + 1
        CROSS_h.GetXaxis().SetBinLabel(bin_x,i_op1)
        CROSS_geom_QUAD_h.GetXaxis().SetBinLabel(bin_x,i_op1)
        CROSS_el_area_ratio_h.GetXaxis().SetBinLabel(bin_x,i_op1)
        for i_op2 in xsec_dict["CROSS"][i_op1]:
            ########### cutflow comparison
            if opts.drawCutflow=="yes":
                i_cutflow_quad1 = cutflow_file_dict["QUAD"][i_op1]
                i_cutflow_quad2 = cutflow_file_dict["QUAD"][i_op2]
                i_cutflow_cross = cutflow_file_dict["CROSS"][i_op1][i_op2]
                #
                cut_names, _, cut_incr_quad1 = lu.get_cutflow_arrays(i_cutflow_quad1)
                _, _, cut_incr_quad2 = lu.get_cutflow_arrays(i_cutflow_quad2)
                _, _, cut_incr_cross = lu.get_cutflow_arrays(i_cutflow_cross)

                lu.draw_cutflows(cut_names, [cut_incr_quad1,cut_incr_quad2,cut_incr_cross],
                                [f"incr {i_op1}",f"incr {i_op2}", f"incr {i_op1}vs{i_op2}"],
                                cutflow_plot_dir+f"/{i_op1}vs{i_op2}.pdf", prod_dec)

    # tables
    lu.save_plot(CROSS_h, plot_dir + "CROSS.pdf")
    lu.save_plot(CROSS_geom_QUAD_h, plot_dir + "CROSS_geom_QUAD.pdf")
    lu.save_plot(CROSS_el_area_ratio_h, plot_dir + "CROSS_el_area_ratio_h.pdf", draw_option = "text colz")

    root_file = ROOT.TFile(plot_dir + "/el_area_ratio.root","UPDATE")  # to divide between ana
    CROSS_el_area_ratio_h.Write("", ROOT.TObject.kOverwrite)
    root_file.Close()
    # fractions
    for i_pair, i_frac_quad1, i_frac_quad2, i_frac_cross in zip(pairs_str_avalaible,fracs_quad1,fracs_quad2,fracs_cross):
        print(f"for pair {i_pair} q1 q2 c fractions are {i_frac_quad1}, {i_frac_quad2}, {i_frac_cross}")

    def draw_efficiency(area_ratio_min):
        new_pairs_str = []
        new_fracs_quad1, new_fracs_errors_quad1 = [], []
        new_fracs_quad2, new_fracs_errors_quad2 = [], []
        new_fracs_cross, new_fracs_errors_cross = [], []
        new_fracs_ave, new_fracs_errors_ave = [], []
        new_fracs_envelope = []
        for i,i_ratio in enumerate(area_ratios):
            if i_ratio < area_ratio_min: continue
            new_pairs_str.append(pairs_str_avalaible[i])
            new_fracs_quad1.append(fracs_quad1[i])
            new_fracs_errors_quad1.append(fracs_errors_quad1[i])
            new_fracs_quad2.append(fracs_quad2[i])
            new_fracs_errors_quad2.append(fracs_errors_quad2[i])
            new_fracs_cross.append(fracs_cross[i])
            new_fracs_errors_cross.append(fracs_errors_cross[i])
            new_fracs_ave.append(fracs_ave[i])
            new_fracs_errors_ave.append(fracs_errors_ave[i])
            new_fracs_envelope.append(fracs_envelope[i])
        new_ops_num = range(len(new_pairs_str))
        #
        plt.clf()
        fig, (ax1, ax2) = plt.subplots(2, figsize=(16,9), gridspec_kw={'height_ratios': [2, 1]})
        ax1.errorbar(new_ops_num, new_fracs_quad1, yerr=new_fracs_errors_quad1, fmt='v', mfc="blue", mec="blue", ecolor="blue", label="QUAD1", alpha=0.7)
        ax1.errorbar(new_ops_num, new_fracs_quad2, yerr=new_fracs_errors_quad2, fmt='^', mfc="green", mec="green", ecolor="green", label="QUAD2", alpha=0.7)
        ax1.errorbar(new_ops_num, new_fracs_cross, yerr=new_fracs_errors_cross, fmt='*', mfc="black", mec="black", ecolor="black", label="CROSS", alpha=0.7)
        ax1.set_ylabel('fraction of weights')
        ax1.set_title(r"$A_{cross}/A_{nocross}>$"+str(area_ratio_min))
        ax1.legend()
        ax1.set_xticks(new_ops_num, minor=False)
        ax1.set_xticklabels(new_pairs_str, fontdict=None, minor=False,fontsize = 9)
        ax1.tick_params(axis='x', labelrotation=90)
        #
        color_ave = 'teal'
        ax2.set_ylabel(r"$C /<Q1,Q2>$", color=color_ave)
        ax2.errorbar(new_ops_num, new_fracs_ave, yerr=new_fracs_errors_ave, mfc=color_ave, mec=color_ave, ecolor=color_ave, fmt='o')
        ax2.tick_params(axis='y', labelcolor=color_ave)
        ax2.set_xticks(new_ops_num, minor=False)
        plt.setp(ax2.get_xticklabels(), visible=False)
        ax2.set_xlabel('operator pair')
        #
        ax3 = ax2.twinx()
        color_env = 'purple'
        ax3.set_ylabel('envelope', color=color_env)
        ax3.scatter(new_ops_num, new_fracs_envelope, color=color_env)
        ax3.tick_params(axis='y', labelcolor=color_env)
        #
        for i_ax in [ax1,ax2]:
            for xc in new_ops_num: i_ax.axvline(x=xc, color='0.3', linestyle='--', linewidth=0.3)
        fig.tight_layout()
        fig.savefig(plot_dir + f"/frac_w_after_cuts_above_{area_ratio_min}.pdf")

    draw_efficiency(0.99)
    draw_efficiency(big_pairs_cutoff)