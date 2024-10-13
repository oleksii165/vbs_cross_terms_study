# go from integral 1 hist to the one normalized to xsec kind of the way it appears in analysis
# not a part of overall chain as may want or not include overflow in norm, change binning etc
import os
import ROOT
import lib_utils as lu

# ana_dirs = {"Wmy_lvy": ["Wy_John", "SR", "m_Wy"], "Wpy_lvy": ["Wy_John", "SR", "m_Wy"]}

# ana_dirs = {"WmZ_lllv": ["WZ_lllv", "SR", "m_WZ"], "WpZ_lllv": ["WZ_lllv", "SR", "m_WZ"]}

# ana_dirs = {"Zy_vvy":["Zy_vvy", "SR", "m_Zy"]}
# ana_dirs = {"ZZ_llll":["ZZ_llll", "SR", "m_ZZ"]}
ana_dirs = {"Zy_lly":["ATLAS_2023_I2663725", "SR", "m_Zy"]}
skip_cross = True
output_counts = False # need to multiply by 139 1/fb or not?

base_dir = "/lapp_data/atlas/kurdysh/vbs_eft_files/"
clips = ["inf", "3000", "2000", "1500", "1000", "700"]
for i_ana_dir in ana_dirs.keys():
    routine, cut, m_vv_name = ana_dirs[i_ana_dir][0], ana_dirs[i_ana_dir][1], ana_dirs[i_ana_dir][2]
    fit_var, fit_bins, re_overflow, exclude_underverflow_from_norm = lu.get_fitted_plot(routine, cut)
    full_ana_dir = f"{base_dir}/{i_ana_dir}/"
    op_dirs = [i_dir for i_dir in os.listdir(full_ana_dir) if "_EXT0" in i_dir]
    for i_op_dir in op_dirs:
        print("working on ", i_op_dir)
        hist_xsec_dir = f"{full_ana_dir}/{i_op_dir}/routine_{routine}_cut_{cut}/"
        if skip_cross:
            if "CROSS" in i_op_dir: continue
        # find hist file and if relevant hists are there
        hist_path_in = f"{hist_xsec_dir}/hists.root"
        if not os.path.exists(hist_path_in):
            print("-------- hfile doesn't exist will skip", hist_path_in)
            continue
        hist_root_in = ROOT.TFile(hist_path_in, "read")
        hists_names_in = [h.GetName() for h in hist_root_in.GetListOfKeys() if h.GetName().startswith(fit_var)]
        if len(hists_names_in)==0:
            print("---cant find any fitted", fit_var, "hist in",hist_root_in)
        var_clip_str = f"{fit_var}_clip_"
        clips_hist_found = [ihn[ihn.find(var_clip_str)+len(var_clip_str):] for ihn in hists_names_in]
        # find xsec files
        xsecs = {} # key is clip
        for i_clip in clips:
            i_xsec_file = f"{hist_xsec_dir}/xsec_times_frac_fb_clip_{i_clip}.txt"
            if not os.path.exists(i_xsec_file):
                print("--cant find xsec file", i_xsec_file)
                continue
            xsecs[i_clip] = lu.read_oneline_file(i_xsec_file)
        print(xsecs)
        # combine hist and xsec
        lumi = 139 if output_counts else 1
        o_file_suff = "norm_run2" if output_counts else "norm_xsec"
        clips_found = list(set(clips_hist_found).intersection(xsecs.keys()))
        hist_path_out = hist_path_in.replace("hists.root",f"hists_{fit_var}_{o_file_suff}.root")
        hist_root_out = ROOT.TFile(hist_path_out,"recreate")
        lumi = 139 if output_counts else 1
        for i_clip_found in clips_found:
            i_h_name = f"{var_clip_str}{i_clip_found}"
            if i_h_name in hists_names_in:
                i_h = hist_root_in.Get(i_h_name)
                i_h_normed = lu.dress_hist(i_h, i_h_name, my_norm = xsecs[i_clip_found]*1,
                                           re_bins = fit_bins,re_overflow = re_overflow,
                                           exclude_underverflow_from_norm = exclude_underverflow_from_norm)
                i_h_normed.Write("", ROOT.TObject.kOverwrite)
        # also save m_vv
        # i_h = hist_root_in.Get(m_vv_name)
        i_h_normed = lu.dress_hist(hist_root_in.Get(m_vv_name), m_vv_name, my_norm = xsecs["inf"]*1)
        i_h_normed.Write("", ROOT.TObject.kOverwrite)
        hist_root_out.Close()




