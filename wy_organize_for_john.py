# take hists from my ususal organization and stack them into folders like john did originally
import os
import lib_utils as lu
import ROOT

my_in_baserdir = "/sps/atlas/kurdysh/vbs_cross_terms_study/eft_files/WpyWmy_lvy/"
out_basedir = "/exp/atlas/kurdysh/wy_share_conventional_hists/"

for op_dir in os.listdir(my_in_baserdir):
    _, order, ops_str = lu.get_op_from_dir(op_dir,"WpyWmy_lvy")
    ops_str = ops_str.replace("F","") # as he did't have F in naming
    print("working with ops and order", ops_str, order)
    h_file = ROOT.TFile(f"{my_in_baserdir}/{op_dir}/routine_Wy_John_cut_SR/hists_norm_run2.root","read")
    for clip in ["inf","3000","2000","1500","1000","700"]:
        clip_h = h_file.Get(f"pt_lepton_clip_{clip}").Clone()
        clip_h.SetName("h_truth_lep_pt")
        out_clip_dir = "nominal" if clip=="inf" else f"clip_{clip}"
        out_dir = f"{out_basedir}/{out_clip_dir}/"
        if not os.path.exists(out_dir): os.makedirs(out_dir)
        out_file = ROOT.TFile(f"{out_dir}/{ops_str}_{order}.root", "recreate")
        clip_h.Write("", ROOT.TObject.kOverwrite)
        out_file.Close()