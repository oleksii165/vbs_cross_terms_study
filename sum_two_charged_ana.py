# for ana like ssWW with w-w- and w+w+ component
# add histograms maitaining fraction between two component (but still normalize to 1)
# save xsec after cuts for all clippings = xsec+ + xsec-

import lib_utils as lu
import os
from glob import glob
import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--sumGenAna", default = "WpyWmy_lvy")
parser.add_option("--routine", default = "Wy_John")
parser.add_option("--cut", default = "SR")
opts, _ = parser.parse_args()
assert opts.sumGenAna in ["WpyWmy_lvy","WpZWmZ_lllv"], "not sure can work with this ana"
_, routine_cut_dir = lu.get_plotdir(opts.sumGenAna, opts.routine, opts.cut)
if opts.sumGenAna == "ssWW_lvlv":
    genAna_plus,genAna_minus = "WpWp_lvlv", "WmWm_lvlv"
elif opts.sumGenAna == "WpZWmZ_lllv":
    genAna_plus,genAna_minus = "WpZ_lllv", "WmZ_lllv"
elif opts.sumGenAna == "WpyWmy_lvy":
    genAna_plus,genAna_minus = "Wpy_lvy", "Wmy_lvy"
debug = False

_, top_files_dir_sum = lu.find_prod_dec_and_dir(f"user.okurdysh.MadGraph_{opts.sumGenAna}_FM0_SM")
_, top_files_dir_plus = lu.find_prod_dec_and_dir(f"user.okurdysh.MadGraph_{genAna_plus}_FM0_SM")
_, top_files_dir_minus = lu.find_prod_dec_and_dir(f"user.okurdysh.MadGraph_{genAna_minus}_FM0_SM")
#arr of elements like ["QUAD","FT7"]
confs_plus = [lu.get_op_from_dir(conf_dir,genAna_plus)[1:] for conf_dir in os.listdir(top_files_dir_plus)]
confs_minus = [lu.get_op_from_dir(conf_dir,genAna_minus)[1:] for conf_dir in os.listdir(top_files_dir_minus)]
confs_common = list(set(confs_plus).intersection(confs_minus))
print("configs plus with len",len(confs_plus),confs_plus)
print("configs minus of len",len(confs_minus),confs_minus )
print("configs in both minus and plus of len",len(confs_common),confs_common)

for i_order_op in confs_common:
    print("looking for input dirs of form", f"{top_files_dir_plus}/*{i_order_op[1]}_{i_order_op[0]}*EXT0/{routine_cut_dir}/")
    i_dir_plus_cand = glob(f"{top_files_dir_plus}/*{i_order_op[1]}_{i_order_op[0]}*EXT0/{routine_cut_dir}/")
    i_dir_minus_cand = glob(f"{top_files_dir_minus}/*{i_order_op[1]}_{i_order_op[0]}*EXT0/{routine_cut_dir}/")
    print("found input dirs", i_dir_plus_cand, i_dir_minus_cand)
    if len(i_dir_plus_cand)!=1 or len(i_dir_minus_cand)!=1: continue
    i_dir_plus, i_dir_minus = i_dir_plus_cand[0], i_dir_minus_cand[0]
    i_rootname_plus, i_rootname_minus = f"{i_dir_plus}/hists_norm_run2.root", f"{i_dir_minus}/hists_norm_run2.root"
    if not os.path.exists(i_rootname_plus) or not os.path.exists(i_rootname_minus): continue
    if debug: print("will work with +- root files", i_rootname_plus, i_rootname_minus)
    # save plots
    i_root_plus, i_root_minus = ROOT.TFile(i_rootname_plus, "read"),ROOT.TFile(i_rootname_minus, "read")
    i_plots_sum = {}
    for i_hist_name in [ih.GetName() for ih in list(i_root_plus.GetListOfKeys())]:
        if debug: print("---working with var", i_hist_name)
        i_hist_minus = i_root_minus.Get(i_hist_name).Clone()
        i_hist_minus.SetName(i_hist_name+"_m") # otherwise as name is the same and root navigates by name +,- will be the same
        i_hist_plus = i_root_plus.Get(i_hist_name).Clone()
        i_hist_plus.SetName(i_hist_name+"_p")
        # as they are already norm to run2 countes just do the sum
        if debug: print("xsec integral of plus and minus", i_hist_plus.Integral(), i_hist_minus.Integral())
        i_hist_sum = i_hist_plus.Clone()
        i_hist_sum.SetName(i_hist_name) # bring back original name
        i_hist_sum.Add(i_hist_minus)
        if debug: print("integral of sum", i_hist_sum.Integral())
        i_hist_sum.SetDirectory(0)
        i_plots_sum[i_hist_name] = i_hist_sum
    i_root_plus.Close()
    i_root_minus.Close()
    # push xsec and plots and sum dir
    i_dir_sum = i_dir_plus.replace(genAna_plus, opts.sumGenAna)
    if not os.path.exists(i_dir_sum): os.makedirs(i_dir_sum)
    i_root_sum = ROOT.TFile(i_dir_sum+"hists_norm_run2.root", "recreate")
    for i_hist in i_plots_sum.values():
        i_hist.Write("", ROOT.TObject.kOverwrite)
    i_root_sum.Close()
