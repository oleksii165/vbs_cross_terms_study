import glob
import ROOT
import myROOTDrawing.hist_ops as mho
import matplotlib.pyplot as plt
import numpy as np
import os

dp = {} # [rebin, min_x, max_x, log_y]
dp["MET_by_HT"] = [1, -999, -999, 0]
dp["dR_ll"] = [1, 0, 2, 0]
dp["dphi_MET_pt_ll"] = [1, 0, 4, 0]
dp["dphi_ll"] = [1, 0, 2, 0]
dp["dphi_tagjets"] = [1, 0, 4, 0]
dp["dy_ll"] = [1, 0, 2, 0]
dp["dy_tagjets"] = [1, 0, 10, 0]
dp["eta_lepton"] = [1, -999, -999, 0]
dp["eta_tagjets"] = [1, -999, -999, 0]
dp["m_ZZ"] = [1, 0, 7000, 0]
dp["m_ll"] = [1, -999, -999, 0]
dp["m_tagjets"] = [1, 0, 3000, 0]
dp["n_jets"] = [1, -999, -999, 0]
dp["n_lepton_stable"] = [1, -999, -999, 0]
dp["phi_tagjets"] = [1, -999, -999, 0]
dp["pt_MET"] = [1, 0, 2500, 0]
dp["pt_Z_clip_inf"] = [1, 0, 2000, 0]
dp["pt_lepton"] = [1, 0, 1500, 0]
dp["pt_tagjet1"] = [1, 0, 1500, 0]
dp["pt_tagjet2"] = [1, 0, 500, 0]


def get_hist(order, op, hist_name, cut, llvv_dir, sample_suff=""):
    assert llvv_dir in ["ZZ_llvv_original", "ZZ_llvv", "ZZ_llvv_comparisons"]
    my_dir = "/lapp_data/atlas/kurdysh/vbs_eft_files/"
    root_file_names = glob.glob(f"{my_dir}/{llvv_dir}/*{sample_suff}_ZZ_llvv*{op}*{order}*/routine_ZZ_llvv_cut_{cut}/hists_norm_xsec.root")
    hist = -1
    if len(root_file_names)==1:
        root_file = ROOT.TFile(root_file_names[0], "READ")
        if hist_name in [ih.GetName() for ih in list(root_file.GetListOfKeys())]:
            hist = root_file.Get(hist_name).Clone()
            hist.SetDirectory(0)
            hist.RebinX(dp[hist_name][0])
    return hist

hists={} # key is order, op, hist_name, cut, original_sample
for i_order in ["QUAD","INT"]:
    for i_op in ["FM0","FM1","FM2","FM3","FM4","FM5","FM7","FS02","FS1","FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]:
        print("getting", i_op, i_order)
        for i_hist_name in dp.keys():
            for i_cut in ["SR","NO"]:
                for llvv_dir in ["ZZ_llvv","ZZ_llvv_original"]:
                    i_ntuple = (i_order,i_op,i_hist_name,i_cut,llvv_dir,"")
                    hists[i_ntuple] = get_hist(*i_ntuple)
# get FT5 for comparisons
suffixes = [
            # "rep_cuts",
            # "rep_cuts_coupling",
            # "rep_cuts_coupling_pdf",
            # "rep_cuts_coupling_pdf_dynscale",
            # "rep_cuts_coupling_pdf_dynscale_qed6",
            # "rep_cuts_coupling_pdf_dynscale_qed6_autoptjmjj_cutdecays"
            "rep_cuts_coupling1_pdf_dynscale_qed6_autoptjmjj_cutdecays_sde_hardstrat",
            "rep_cuts_coupling1_pdf_dynscale_qed6_autoptjmjj_cutdecays_sde_hardstrat_r21"
            ]
for i_hist_name in dp.keys():
    for i_cut in ["SR","NO"]:
        for i_suff in suffixes:
            i_ntuple = ("QUAD","FT5",i_hist_name,i_cut,"ZZ_llvv_comparisons",i_suff)
            hists[i_ntuple] = get_hist(*i_ntuple)

print("hi")

plot_basedir = "/lapp_data/atlas/kurdysh/vbs_truth/plots/ZZ_llvv/"
# compare original and current
# for i_op in ["FT0","FT2","FT5","FT8"]:
#     odir = f"{plot_basedir}/original_new_comparison/{i_op}/"
#     if not os.path.exists(odir): os.makedirs(odir)
#     for i_order in ["INT","QUAD"]:
#         for i_hist_name in dp.keys():
#             x_min, x_max = dp[i_hist_name][1], dp[i_hist_name][2]
#             for i_cut in ["SR","NO"]:
#                 h_original_counts, h_original_err, h_x = mho.hist_to_arr(hists[(i_order,i_op,i_hist_name,i_cut,"ZZ_llvv_original","")])
#                 h_new_counts, h_new_err, _ = mho.hist_to_arr(hists[(i_order,i_op,i_hist_name,i_cut,"ZZ_llvv","")])
#                 ratio_new_to_old = np.divide(np.array(h_new_counts), np.array(h_original_counts))
#                 #
#                 plt.clf()
#                 fig, ax = plt.subplots(2,1, gridspec_kw={'height_ratios': [2, 1]}, sharex='col')
#                 fig.tight_layout()
#                 ax[0].plot(h_x, h_original_counts, 'o', label="original")
#                 ax[0].plot(h_x, h_new_counts, 'o', label="new")
#                 ax[0].set_ylabel("diff xsec")
#                 ax[0].set_title(f"{i_order} {i_op} {i_hist_name}, cut:{i_cut}")
#                 ax[0].legend()
#                 if x_min!=-999: ax[0].set_xlim([x_min, x_max])
#                 #
#                 ax[1].plot(h_x, ratio_new_to_old, 'o')
#                 ax[1].set_ylabel("new/original")
#                 if x_min!=-999: ax[1].set_xlim([x_min, x_max])
#                 oname = odir + f"{i_order}_cut_{i_cut}_{i_hist_name}.png"
#                 print("saving", oname)
#                 plt.savefig(oname, bbox_inches='tight')

#compare different versions of new FT5
for i_op in ["FT5"]:
    odir = f"{plot_basedir}/mg_settings_comparison/{i_op}/"
    if not os.path.exists(odir): os.makedirs(odir)
    for i_order in ["QUAD"]:
        for i_hist_name in dp.keys():
            x_min, x_max = dp[i_hist_name][1], dp[i_hist_name][2]
            # for i_cut in ["SR","NO"]:
            for i_cut in ["NO"]:
                h_original_counts, h_original_err, h_x = mho.hist_to_arr(hists[(i_order,i_op,i_hist_name,i_cut,"ZZ_llvv_original","")])
                arr_counts_per_setting = {}
                arr_err_per_setting = {}
                arr_ratio_to_original_per_setting = {}
                for i_setting in suffixes:
                    h_counts, h_err, _ = mho.hist_to_arr(hists[(i_order,i_op,i_hist_name,i_cut,"ZZ_llvv_comparisons",i_setting)])
                    arr_counts_per_setting[i_setting] = h_counts
                    arr_err_per_setting[i_setting] = h_err
                    arr_ratio_to_original_per_setting[i_setting] = np.divide(np.array(h_counts), np.array(h_original_counts))
                # draw all together
                plt.clf()
                colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
                fig, ax = plt.subplots(2,1,
                                       # gridspec_kw={'height_ratios': [2, 1]},
                                       sharex='col')
                fig.tight_layout()
                ax[0].plot(h_x, h_original_counts, 'o', label="original",color="black")
                for num_color,i_setting in enumerate(suffixes):
                    ax[0].plot(h_x, arr_counts_per_setting[i_setting],'o',label=i_setting,color=colors[num_color])
                ax[0].set_ylabel("diff xsec")
                ax[0].set_title(f"{i_order} {i_op} {i_hist_name}, cut:{i_cut}")
                ax[0].legend()
                if x_min!=-999: ax[0].set_xlim([x_min, x_max])
                #
                for num_color,i_setting in enumerate(suffixes):
                    ax[1].plot(h_x,arr_ratio_to_original_per_setting[i_setting],'o',label=i_setting,color=colors[num_color])
                ax[1].set_ylabel("new/original")
                if x_min!=-999: ax[1].set_xlim([x_min, x_max])
                # ax[1].set_ylim([0, 1])
                oname = odir + f"{i_order}_cut_{i_cut}_{i_hist_name}.png"
                print("saving", oname)
                plt.savefig(oname, bbox_inches='tight')

# dump kinematics of new plots
for i_op in ["FM0","FM1","FM2","FM3","FM4","FM5","FM7","FS02","FS1","FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]:
    odir = f"{plot_basedir}/kinematics/{i_op}/"
    if not os.path.exists(odir): os.makedirs(odir)
    for i_order in ["INT","QUAD"]:
        for i_hist_name in dp.keys():
            x_min, x_max = dp[i_hist_name][1], dp[i_hist_name][2]
            for i_cut in ["NO"]:
                i_hist = hists[(i_order,i_op,i_hist_name,i_cut,"ZZ_llvv","")]
                if i_hist==-1: continue
                h_new_counts, h_new_err, h_x = mho.hist_to_arr(i_hist)
                #
                plt.clf()
                fig, ax = plt.subplots(1,1)
                fig.tight_layout()
                ax.plot(h_x, h_new_counts, 'o', label="new")
                ax.set_ylabel("diff xsec")
                ax.set_title(f"{i_order} {i_op} {i_hist_name}, cut:{i_cut}")
                ax.legend()
                if x_min!=-999: ax.set_xlim([x_min, x_max])
                oname = odir + f"{i_order}_cut_{i_cut}_{i_hist_name}.png"
                print("saving", oname)
                plt.savefig(oname, bbox_inches='tight')