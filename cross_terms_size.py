import os
import lib_utils as lu
import itertools
import ROOT
import myROOTDrawing.utils as mrd
import math
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--tProd", default = "Zy")
parser.add_option("--tDec", default = "vvy")
opts, _ = parser.parse_args()

prod_dec = f"{opts.tProd}_{opts.tDec}"
assert prod_dec in ["Zy_vvy", "Wmy_lvy", "WmZ_lllv"]
top_files_dir = f"/sps/atlas/kurdysh/vbs_cross_terms_study/eft_files/{prod_dec}_all_crosses/"
docut_dir = "DOCUT_YES"
plot_dir = f"/exp/atlas/kurdysh/vbs_cross_terms_study/plots/{prod_dec}/{docut_dir}/"

xsec_dict = {'SM':{}, 'FULL':{}, 'QUAD':{}, 'CROSS':{}, 'INT':{}}
for op_dir in [i_obj for i_obj in os.listdir(top_files_dir) if os.path.isdir(top_files_dir + "/" + i_obj)]:
    full_op_dir = os.path.join(top_files_dir,op_dir,docut_dir,"")
    print("#### plotting new conf", op_dir, "full path", full_op_dir)
    ops_arr, regime, _ = lu.get_op_from_dir(op_dir, prod_dec)
    print("this is for ops", ops_arr, "in regime", regime)
    product_file = full_op_dir + "xsec_times_frac_fb.txt"
    if not os.path.exists(product_file):
        print("didnt find the txt xsec fid and/or root file or frac file or log file for", op_dir)
        continue
    # xsec organization
    with open(product_file, 'r') as f: fid_xsec_fb = float(f.read())
    print("reading back fid_xsec in fb", fid_xsec_fb)
    if regime=="CROSS":
        my_op1, my_op2  = ops_arr[0], ops_arr[1]
        if my_op1 not in xsec_dict[regime].keys(): xsec_dict[regime][my_op1] = {}
        xsec_dict[regime][my_op1][my_op2] = fid_xsec_fb
    else:
        my_op = ops_arr[0]
        xsec_dict[regime][my_op] = fid_xsec_fb

all_ops =  sorted(list(xsec_dict["QUAD"].keys()))
print("all ops", all_ops)
nbins = len(all_ops)

CROSS_geom_QUAD_h = ROOT.TH2F("CROSS_geom_QUAD_h","CROSS_geom_QUAD_h", nbins,0,nbins,nbins,0,nbins)
CROSS_el_area_ratio_h = ROOT.TH2F("CROSS_el_area_ratio_h","CROSS_el_area_ratio_h", nbins,0,nbins,nbins,0,nbins)
pairs_filled = []
for i_op1 in all_ops:
    if i_op1 not in xsec_dict["CROSS"].keys(): continue

    bin_x = all_ops.index(i_op1) + 1
    CROSS_el_area_ratio_h.GetXaxis().SetBinLabel(bin_x,i_op1)
    for i_op2 in xsec_dict["CROSS"][i_op1]:
        if i_op2 not in all_ops: continue
        bin_y = all_ops.index(i_op2) + 1
        CROSS_el_area_ratio_h.GetYaxis().SetBinLabel(bin_y,i_op2)
        cross = xsec_dict["CROSS"][i_op1][i_op2]
        # fill geometric average

        if i_op1 not in xsec_dict["QUAD"].keys() or i_op2 not in xsec_dict["QUAD"].keys(): continue

        quad1 = xsec_dict["QUAD"][i_op1]
        quad2 = xsec_dict["QUAD"][i_op2]
        ##### ellipses
        lumi = 140 # lumi of run2 in 1/fb, expect xsec to be in fb
        el_q1_c = lumi * quad1 / 3 # dont square xsec since already squared from madgraph
        el_q2_c = lumi * quad2 / 3 # divide by 3 since formula for are expect factor to be 1
        el_cross_c = lumi * cross / 3
        den_with_cross2 = 4*el_q1_c*el_q2_c - el_cross_c**2
        area_ratio = 1.0 # placeholder
        if den_with_cross2 > 0:
            area_no_cross = 2 * math.pi / math.sqrt(4*el_q1_c*el_q2_c)
            area_with_cross = 2 * math.pi / math.sqrt(den_with_cross2)
            area_ratio = area_with_cross/area_no_cross
            print(f"for i_op1 i_op2 {i_op1} {i_op2} got areas no cross {area_no_cross:.2f} with cross {area_with_cross:.2f} ratio {area_ratio:.2f} used ellipse q1 q2 cross {el_q1_c:.2f} {el_q2_c:.2f} {el_cross_c:.2f}")
        CROSS_el_area_ratio_h.SetBinContent(bin_x, bin_y, round(area_ratio,2))
        pairs_filled.append([i_op1, i_op2])
# clean up with 1s
for i_all_pair in list(itertools.combinations(all_ops,2)):
    if list(i_all_pair) in pairs_filled or [i_all_pair[1],i_all_pair[0]] in pairs_filled: continue
    bin_x = all_ops.index(i_all_pair[0]) + 1
    bin_y = all_ops.index(i_all_pair[1]) + 1
    CROSS_el_area_ratio_h.SetBinContent(bin_x, bin_y, 1.0)

# tables
mrd.one_plot([CROSS_el_area_ratio_h],
             f"{plot_dir}/CROSS_el_area_ratio_h.pdf", make_out_missing_dirs = True,
             draw_opt = "text colz",
             # canvas_aspect=[5,3],
             )
