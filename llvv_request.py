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

ops_arr = ["FM0","FM1","FM2","FM3","FM4","FM5","FM7",
           "FS02",
           # "FS1",
           "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]
cuts_arr = ["NO"]
orders_arr = ["QUAD","INT"]

def get_hist(order, op, hist_name, cut, llvv_dir, sample_suff=""):
    assert llvv_dir in ["ZZ_llvv_original", "ZZ_llvv", "ZZ_llvv_comparisons", "ZZ_llvv_request", "ZZ_llvv_request_comparisons"]
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
for i_order in orders_arr:
    for i_op in ops_arr:
        print("getting", i_op, i_order)
        for i_hist_name in dp.keys():
            for i_cut in cuts_arr:
                for llvv_dir in ["ZZ_llvv_request"]:
                    i_ntuple = (i_order,i_op,i_hist_name,i_cut,llvv_dir,"")
                    hists[i_ntuple] = get_hist(*i_ntuple)
# get FT5 for comparisons
suffixes = [
    "requestJO",
    "100500_pdf",
    "100501_pdf_dynscale",
    "100502_pdf_dynscale_qed6",
    "100503_pdf_dynscale_qed6_cutdecaysautoptjmjj",
    "100504_pdf_dynscale_qed6_cutdecaysautoptjmjj_cuts",
    "legacyJO"
]
suffixes_labels = [
    "requestJO",
    "requestJO+change1",
    "requestJO+change1,2",
    "requestJO+change1,2,3",
    "requestJO+change1,2,3,4",
    "requestJO+change1,2,3,4,5",
    "legacyJO"
]
for i_hist_name in dp.keys():
    for i_cut in cuts_arr:
        for i_suff in suffixes:
            i_ntuple = ("QUAD","FT5",i_hist_name,i_cut,"ZZ_llvv_request_comparisons",i_suff)
            hists[i_ntuple] = get_hist(*i_ntuple)

print("hi")

plot_basedir = "/lapp_data/atlas/kurdysh/vbs_truth/plots/ZZ_llvv/"

#compare different versions of new FT5
for i_op in ["FT5"]:
    odir = f"{plot_basedir}/mg_settings_comparison/{i_op}/"
    if not os.path.exists(odir): os.makedirs(odir)
    for i_order in ["QUAD"]:
        for i_hist_name in dp.keys():
            x_min, x_max = dp[i_hist_name][1], dp[i_hist_name][2]
            for i_cut in cuts_arr:
                h_original_counts, h_original_err, h_x = mho.hist_to_arr(hists[(i_order,i_op,i_hist_name,i_cut,"ZZ_llvv_request_comparisons","legacyJO")])
                arr_counts_per_setting = {}
                arr_err_per_setting = {}
                arr_ratio_to_original_per_setting = {}
                for i_setting in suffixes:
                    if i_setting=="legacyJO":
                        continue
                    h_counts, h_err, _ = mho.hist_to_arr(hists[(i_order,i_op,i_hist_name,i_cut,"ZZ_llvv_request_comparisons",i_setting)])
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
                ax[0].errorbar(h_x, h_original_counts, yerr=h_original_err, marker='o', label="legacyJO",color="black")
                for num_color,i_setting in enumerate(suffixes):
                    if i_setting=="legacyJO":
                        continue
                    ax[0].errorbar(h_x, arr_counts_per_setting[i_setting], yerr=arr_err_per_setting[i_setting], marker='o',
                                   label = suffixes_labels[suffixes.index(i_setting)],
                                   color=colors[num_color])
                ax[0].set_ylabel("diff xsec")
                ax[0].set_title(f"{i_order} {i_op} {i_hist_name}, cut:{i_cut}")
                ax[0].legend()
                if x_min!=-999: ax[0].set_xlim([x_min, x_max])
                #
                for num_color,i_setting in enumerate(suffixes):
                    if i_setting=="legacyJO":
                        continue
                    ax[1].plot(h_x,arr_ratio_to_original_per_setting[i_setting],marker='o',color=colors[num_color])
                ax[1].set_ylim([0.5, 1.5])
                ax[1].set_ylabel("new/legacy")
                if x_min!=-999: ax[1].set_xlim([x_min, x_max])
                # ax[1].set_ylim([0, 1])
                oname = odir + f"{i_order}_cut_{i_cut}_{i_hist_name}.png"
                print("saving", oname)
                plt.savefig(oname, bbox_inches='tight')

# dump kinematics of new plots
for i_op in ops_arr:
    odir = f"{plot_basedir}/kinematics/{i_op}/"
    if not os.path.exists(odir): os.makedirs(odir)
    for i_order in orders_arr:
        for i_hist_name in dp.keys():
            x_min, x_max = dp[i_hist_name][1], dp[i_hist_name][2]
            for i_cut in cuts_arr:
                i_hist = hists[(i_order,i_op,i_hist_name,i_cut,"ZZ_llvv_request","")]
                if i_hist==-1: continue
                h_new_counts, h_new_err, h_x = mho.hist_to_arr(i_hist)
                #
                plt.clf()
                fig, ax = plt.subplots(1,1)
                fig.tight_layout()
                ax.errorbar(h_x, h_new_counts, yerr=h_new_err, marker='o', label="new")
                ax.set_ylabel("diff xsec")
                ax.set_title(f"{i_order} {i_op} {i_hist_name}, cut:{i_cut}")
                ax.legend()
                if x_min!=-999: ax.set_xlim([x_min, x_max])
                oname = odir + f"{i_order}_cut_{i_cut}_{i_hist_name}.png"
                print("saving", oname)
                plt.savefig(oname, bbox_inches='tight')