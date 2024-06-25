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
parser.add_option("--sumGenAna", default = "ssWW_lvlv")
parser.add_option("--routine", default = "ssWW_lvlv")
parser.add_option("--cut", default = "SR")
opts, _ = parser.parse_args()
_, routine_cut_dir = lu.get_plotdir(opts.sumGenAna, opts.routine, opts.cut)
if opts.sumGenAna == "ssWW_lvlv":
    genAna_plus,genAna_minus = "WpWp_lvlv", "WmWm_lvlv"
elif opts.sumGenAna == "WZ_lllv":
    genAna_plus,genAna_minus = "WpZ_lllv", "WmZ_lllv"
elif opts.sumGenAna == "Wy_lvy":
    genAna_plus,genAna_minus = "Wpy_lvy", "Wmy_lvy"
clips_not_inf = ["3000","2000", "1500", "1000", "700"]
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

def read_xsec(xsec_file):
    with open(xsec_file, 'r') as f:
        fid_xsec_fb = float(f.read())
    return fid_xsec_fb

for i_order_op in confs_common:
    print("looking for input dirs of form", f"{top_files_dir_plus}/*{i_order_op[1]}_{i_order_op[0]}*EXT0/{routine_cut_dir}/")
    i_dir_plus_cand = glob(f"{top_files_dir_plus}/*{i_order_op[1]}_{i_order_op[0]}*EXT0/{routine_cut_dir}/")
    i_dir_minus_cand = glob(f"{top_files_dir_minus}/*{i_order_op[1]}_{i_order_op[0]}*EXT0/{routine_cut_dir}/")
    print("found input dirs", i_dir_plus_cand, i_dir_minus_cand)
    if len(i_dir_plus_cand)!=1 or len(i_dir_minus_cand)!=1: continue
    i_dir_plus, i_dir_minus = i_dir_plus_cand[0], i_dir_minus_cand[0]
    i_rootname_plus, i_rootname_minus = f"{i_dir_plus}/hists.root", f"{i_dir_minus}/hists.root"
    if not os.path.exists(i_rootname_plus) or not os.path.exists(i_rootname_minus): continue
    if debug: print("will work with +- root files", i_rootname_plus, i_rootname_minus)
    i_xsecname_plus_inf, i_xsecname_minus_inf = f"{i_dir_plus}/xsec_times_frac_fb_clip_inf.txt", f"{i_dir_minus}/xsec_times_frac_fb_clip_inf.txt"
    i_seleff_plus_inf, i_seleff_minus_inf = f"{i_dir_plus}/frac_after_cuts_clip_inf.txt", f"{i_dir_minus}/frac_after_cuts_clip_inf.txt"
    i_seleff_error_plus_inf, i_seleff_error_minus_inf = f"{i_dir_plus}/frac_after_cuts_error_bar_clip_inf.txt", f"{i_dir_minus}/frac_after_cuts_error_bar_clip_inf.txt"
    if not os.path.exists(i_xsecname_plus_inf) or not os.path.exists(i_xsecname_minus_inf): continue
    if not os.path.exists(i_seleff_plus_inf) or not os.path.exists(i_seleff_minus_inf): continue
    if not os.path.exists(i_seleff_error_plus_inf) or not os.path.exists(i_seleff_error_minus_inf): continue
    # save xsecs
    i_xsecs_plus, i_xsecs_minus = {},{}
    i_seleff_plus, i_seleff_minus = {},{}
    i_seleff_error_plus, i_seleff_error_minus = {},{}
    i_xsecs_plus["inf"],i_xsecs_minus["inf"] = read_xsec(i_xsecname_plus_inf), read_xsec(i_xsecname_minus_inf)
    i_seleff_plus["inf"],i_seleff_minus["inf"] = read_xsec(i_seleff_plus_inf), read_xsec(i_seleff_minus_inf)
    i_seleff_error_plus["inf"],i_seleff_error_minus["inf"] = read_xsec(i_seleff_error_plus_inf), read_xsec(i_seleff_error_minus_inf)
    for i_clip in clips_not_inf:
        i_xsecs_plus[i_clip] = read_xsec(i_xsecname_plus_inf.replace("clip_inf",f"clip_{i_clip}"))
        i_xsecs_minus[i_clip] = read_xsec(i_xsecname_minus_inf.replace("clip_inf",f"clip_{i_clip}"))
        i_seleff_plus[i_clip] = read_xsec(i_seleff_plus_inf.replace("clip_inf",f"clip_{i_clip}"))
        i_seleff_minus[i_clip] = read_xsec(i_seleff_minus_inf.replace("clip_inf",f"clip_{i_clip}"))
        i_seleff_error_plus[i_clip] = read_xsec(i_seleff_error_plus_inf.replace("clip_inf",f"clip_{i_clip}"))
        i_seleff_error_minus[i_clip] = read_xsec(i_seleff_error_minus_inf.replace("clip_inf",f"clip_{i_clip}"))

    if debug: print("got pos xsecs", i_xsecs_plus)
    if debug: print("got min xsecs",i_xsecs_minus)
    # save plots
    i_root_plus, i_root_minus = ROOT.TFile(i_rootname_plus, "read"),ROOT.TFile(i_rootname_minus, "read")
    i_plots_sum = {}
    for i_hist_name in [ih.GetName() for ih in list(i_root_plus.GetListOfKeys())]:
        if debug: print("---working with var", i_hist_name)
        i_hist_plus, i_hist_minus = i_root_plus.Get(i_hist_name).Clone(), i_root_minus.Get(i_hist_name).Clone()
        i_hist_plus_int, i_hist_minus_int = i_hist_plus.Integral(), i_hist_minus.Integral()
        if i_hist_plus_int==0 or i_hist_minus_int==0:
            print("################### skip as for var clip op", i_hist_name, i_clip, i_order_op, "have intergral zero in one of charges")
            continue
        else:
        # in the end will norm one but for now norm to xsec to get relative fraction right
        # assume fraction at differenct clippings is the same at clip=inf
            # as intial norm is typically  1 but can be not 1 like for Wy because of overflow bin in norm
            i_hist_plus.Scale(i_xsecs_plus["inf"]/i_hist_plus_int) 
            i_hist_minus.Scale(i_xsecs_minus["inf"]/i_hist_minus_int)
            if debug: print("after scaling by xsec integral of plus and minus", i_hist_plus.Integral(), i_hist_minus.Integral())
            i_hist_sum = i_hist_plus.Clone()
            i_hist_sum.Add(i_hist_minus)
            if debug: print("integral of sum", i_hist_sum.Integral())
            i_hist_sum.Scale(1/i_hist_sum.Integral())        
            if debug: print("integral of sum after norm to 1", i_hist_sum.Integral())
            i_hist_sum.SetDirectory(0)
            i_plots_sum[i_hist_name] = i_hist_sum
    i_root_plus.Close()
    i_root_minus.Close()
    # push xsec and plots and sum dir
    i_dir_sum = i_dir_plus.replace(genAna_plus, opts.sumGenAna)
    if not os.path.exists(i_dir_sum): os.makedirs(i_dir_sum)
    i_root_sum = ROOT.TFile(i_dir_sum+"hists.root", "recreate")
    for i_hist in i_plots_sum.values():
        i_hist.Write("", ROOT.TObject.kOverwrite)
    i_root_sum.Close()
    #
    for i_clip in i_xsecs_plus.keys():
        i_clip_xsec_sum = float(i_xsecs_plus[i_clip]) + float(i_xsecs_minus[i_clip])
        i_clip_xsecname_sum = f"{i_dir_sum}/xsec_times_frac_fb_clip_{i_clip}.txt"
        lu.write_to_f(i_clip_xsecname_sum, i_clip_xsec_sum)
    for i_clip in i_seleff_plus.keys():
        i_clip_seleff_ave = (float(i_seleff_plus[i_clip]) + float(i_seleff_minus[i_clip])) / 2
        lu.write_to_f(f"{i_dir_sum}/frac_after_cuts_clip_{i_clip}.txt", i_clip_seleff_ave)
        i_clip_seleff_error_ave = (float(i_seleff_error_plus[i_clip]) + float(i_seleff_error_minus[i_clip])) / 2
        lu.write_to_f(f"{i_dir_sum}/frac_after_cuts_error_bar_clip_{i_clip}.txt", i_clip_seleff_error_ave)






