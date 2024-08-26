# take all hists in root file, maybe rebin and cut axis, draw with nicer labels
import os
import lib_utils as lu
import myROOTDrawing.utils as mrd
import ROOT
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--topFilesDir", default = "ZZ_llll")
parser.add_option("--tGenProdDec", default = "ZZ_llll")
parser.add_option("--rivetOutDir", default = "routine_ZZ_llll_cut_SR")
opts, _ = parser.parse_args()

base_dir = "/sps/atlas/kurdysh/vbs_cross_terms_study/eft_files/" + opts.topFilesDir
pdf_out_dir = f"/exp/atlas/kurdysh/vbs_cross_terms_study/plots/{opts.tGenProdDec}/kinematics/"
display_params = lu.get_hists_bounds_cuts(opts.tGenProdDec)

for op_dir in [i_obj for i_obj in os.listdir(base_dir) if os.path.isdir(base_dir + "/" + i_obj)]:
    full_dir = os.path.join(base_dir, op_dir, opts.rivetOutDir, "")
    print("#### plotting new conf", op_dir, "full path to rivet", full_dir)
    _, regime, ops_str = lu.get_op_from_dir(op_dir, opts.tGenProdDec)
    draw_dir = pdf_out_dir + f"{regime}_{ops_str}/"
    if regime=="CROSS": continue
    if not os.path.exists(f"{full_dir}/hists.root"): continue
    root_in_file = ROOT.TFile(f"{full_dir}/hists.root", "READ")
    # root_out_file = ROOT.TFile(full_dir + "/hists_nicified.root", "UPDATE")
    for i_hist_name in [ih.GetName() for ih in list(root_in_file.GetListOfKeys())]:
        i_h = root_in_file.Get(i_hist_name).Clone()
        i_h.SetDirectory(0)
        i_display_params = [1, i_h.GetXaxis().GetXmin(), i_h.GetXaxis().GetXmax()]
        if i_hist_name in display_params.keys():
            i_display_params = display_params[i_hist_name]
        mrd.one_plot([i_h],
                     f"{draw_dir}/{i_hist_name}.pdf", make_out_missing_dirs = True,
                     draw_opt = "hist",
                     x_label = lu.get_var_latex(i_hist_name),
                     canvas_aspect=[5,3],
                     rebins_x_arr = [i_display_params[0]],
                     custom_x_range = [i_display_params[1], i_display_params[2]]
                     )
        # i_h.Write("", ROOT.TObject.kOverwrite)
    root_in_file.Close()
    # root_out_file.Close()
