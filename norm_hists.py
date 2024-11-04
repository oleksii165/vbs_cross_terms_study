# go from integral 1 hist to the one normalized to xsec kind of the way it appears in analysis
# not a part of overall chain as may want or not include overflow in norm, change binning etc
import os
import ROOT
import lib_utils as lu
import yoda
import numpy as np

def get_hist_yodaname(varname, clip, uncert=""):
    yodaname = f"/ATLAS_2023_I2663725:cut=SR/{varname}_clip_{clip}[_{uncert}_]"
    if len(uncert)==0:
        yodaname = yodaname.replace("[_","").replace("_]","")
    return yodaname
uncert_RF = ["MUR05_MUF05_PDF260000", "MUR05_MUF10_PDF260000", "MUR05_MUF20_PDF260000",
             "MUR10_MUF05_PDF260000", "MUR10_MUF10_PDF260000", "MUR10_MUF20_PDF260000",
             "MUR20_MUF05_PDF260000", "MUR20_MUF10_PDF260000", "MUR20_MUF20_PDF260000"]
uncert_PDF = [f"MUR10_MUF10_PDF{PDF}" for PDF in
              [
13100,
25200,
260000,
260001,
260002,
260003,
260004,
260005,
260006,
260007,
260008,
260009,
260010,
260011,
260012,
260013,
260014,
260015,
260016,
260017,
260018,
260019,
260020,
260021,
260022,
260023,
260024,
260025,
260026,
260027,
260028,
260029,
260030,
260031,
260032,
260033,
260034,
260035,
260036,
260037,
260038,
260039,
260040,
260041,
260042,
260043,
260044,
260045,
260046,
260047,
260048,
260049,
260050,
260051,
260052,
260053,
260054,
260055,
260056,
260057,
260058,
260059,
260060,
260061,
260062,
260063,
260064,
260065,
260066,
260067,
260068,
260069,
260070,
260071,
260072,
260073,
260074,
260075,
260076,
260077,
260078,
260079,
260080,
260081,
260082,
260083,
260084,
260085,
260086,
260087,
260088,
260089,
260090,
260091,
260092,
260093,
260094,
260095,
260096,
260097,
260098,
260099,
260100,
265000,
266000,
90400,
90401,
90402,
90403,
90404,
90405,
90406,
90407,
90408,
90409,
90410,
90411,
90412,
90413,
90414,
90415,
90416,
90417,
90418,
90419,
90420,
90421,
90422,
90423,
90424,
90425,
90426,
90427,
90428,
90429,
90430,
90431,
90432]
]
uncert_DYNSCALE = [f"MUR10_MUF20_PDF260000_{scale}" for scale in
                  [
                      # "DYNSCALE1",
                   "DYNSCALE2","DYNSCALE3",
                   # "DYNSCALE4"
                  ]
                  ]
def get_efficiency(yoda_f, fit_var, i_clip_found, unc):
    ref_histname = get_hist_yodaname(fit_var, i_clip_found, uncert=unc)
    pos_in = yoda_f[ref_histname.replace(fit_var,"pos_w_initial").replace(f"_clip_{i_clip_found}","")].sumW()
    neg_in = yoda_f[ref_histname.replace(fit_var,"neg_w_initial").replace(f"_clip_{i_clip_found}","")].sumW()
    pos_out = yoda_f[ref_histname.replace(fit_var,"pos_w_final")].sumW()
    neg_out = yoda_f[ref_histname.replace(fit_var,"neg_w_final")].sumW()
    eff = (pos_in+neg_in) / (pos_out+neg_out)
    return eff
def get_envelope_hist(ref_hist, var_hists, unc_name, clip):
    max_env_per_bin = [0.0] * ref_hist.GetNbinsX()
    for ivar_hist in var_hists:
        for ibin in range(1, ref_hist.GetNbinsX()+1):
            ref_count = ref_hist.GetBinContent(ibin)
            var_count =  ivar_hist.GetBinContent(ibin)
            env = abs(ref_count-var_count)
            if abs(ref_count-var_count) > max_env_per_bin[ibin-1]:
                max_env_per_bin[ibin-1] = env
        # print("after ", ivar_hist.GetName(), "now env will be", max_env_per_bin)
    o_h_name = f"{unc_name}_envelope_on_top_of_nom_clip_{clip}"
    env_on_top_of_nom_hist = ROOT.TH1F(o_h_name, o_h_name,
                                       ref_hist.GetNbinsX(),
                                       ref_hist.GetXaxis().GetXmin(), ref_hist.GetXaxis().GetXmax())
    for ind in range(len(max_env_per_bin)):
        env_on_top_of_nom_hist.SetBinContent(ind+1, max_env_per_bin[ind]+ref_hist.GetBinContent(ind+1))
    return env_on_top_of_nom_hist

# ana_dirs = {"Wmy_lvy": ["Wy_John", "SR", "m_Wy"], "Wpy_lvy": ["Wy_John", "SR", "m_Wy"]}

# ana_dirs = {"WmZ_lllv": ["WZ_lllv", "SR", "m_WZ"], "WpZ_lllv": ["WZ_lllv", "SR", "m_WZ"]}

# ana_dirs = {"Zy_vvy":["Zy_vvy", "SR", "m_Zy"]}
# ana_dirs = {"ZZ_llll":["ZZ_llll", "SR", "m_ZZ"]}
ana_dirs = {"Zy_lly":["ATLAS_2023_I2663725", "SR", "m_Zy"]}
skip_cross = False
output_counts = True # need to multiply by 139 1/fb or not?

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
        for i_clip_found in clips_found:
            i_h_name = f"{var_clip_str}{i_clip_found}"
            if i_h_name in hists_names_in:
                i_h = hist_root_in.Get(i_h_name)
                i_h_normed = lu.dress_hist(i_h, i_h_name, my_norm = xsecs[i_clip_found]*lumi,
                                           re_bins = fit_bins,re_overflow = re_overflow,
                                           exclude_underverflow_from_norm = exclude_underverflow_from_norm)
                i_h_normed.Write("", ROOT.TObject.kOverwrite)
                ########
                # save theory uncert for EFTFun of zlly
                #######
                if i_ana_dir=="Zy_lly":
                    yoda_path = f"{hist_xsec_dir}/MyOutput.yoda.gz"
                    yoda_f = yoda.read(yoda_path)
                    for uncerts_arr, uncert_class_name in zip([uncert_RF,uncert_PDF,uncert_DYNSCALE], ["RF","PDF","DYNSCALE"]):
                        hists_names_unc = [get_hist_yodaname(fit_var, i_clip_found, uncert=unc) for unc in uncerts_arr]
                        hists_unc = [lu.yoda_to_root_1d(yoda_f[hist_name_unc], f"{name_unc}")
                                     for hist_name_unc, name_unc in zip(hists_names_unc, uncerts_arr)]
                        eff_nominal = get_efficiency(yoda_f, fit_var, i_clip_found, "")
                        effs_unc = [get_efficiency(yoda_f, fit_var, i_clip_found, unc) for unc in uncerts_arr]
                        effs_unc_rel_to_nom = np.array(effs_unc) / eff_nominal
                        hists_unc_normed = []
                        for i_h_unc, i_w_factor_unc in zip(hists_unc,effs_unc_rel_to_nom):
                            i_h_unc_normed = lu.dress_hist(i_h_unc, f"{i_h_unc.GetName()}_normed_clip_{i_clip_found}",
                                                           my_norm = xsecs[i_clip_found] * lumi * i_w_factor_unc,
                                                           re_bins = fit_bins,re_overflow = re_overflow,
                                                           exclude_underverflow_from_norm = exclude_underverflow_from_norm)
                            hists_unc_normed.append(i_h_unc_normed)
                            i_h_unc_normed.Write("", ROOT.TObject.kOverwrite)
                        i_env_hist = get_envelope_hist(i_h_normed, hists_unc_normed, uncert_class_name, i_clip_found)
                        i_env_hist.Write("", ROOT.TObject.kOverwrite)
        # also save m_vv
        i_h_normed = lu.dress_hist(hist_root_in.Get(m_vv_name), m_vv_name, my_norm = xsecs["inf"]*lumi)
        i_h_normed.Write("", ROOT.TObject.kOverwrite)
        hist_root_out.Close()
        print("hi")




