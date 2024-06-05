# at the moment have clip hists for ssWW Zvvy and ZZllvv
import lib_utils as lu
import os
import ROOT
import pandas as pd
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gROOT.ProcessLine( "gErrorIgnoreLevel = 1001;")

h_name = {"ZZ_llvv": "m_ZZ",
          "ssWW_lvlv": "m_diboson",
          "Zy_vvy": "m_Zy",
          "WZ_lllv":"m_WZ"}
plot_folder = {"ZZ_llvv": "DOCUT_YES", # didnt update it to new format yet
               "ssWW_lvlv": "routine_ssWW_lvlv_cut_SR",
               "Zy_vvy": "routine_Zy_vvy_cut_SR",
               "WZ_lllv": "routine_WZ_lllv_cut_SR"}

def x_val_for_y(cdf_h1, y_val, tolerance=0.01):
    x_val = -1
    for ibin in range(cdf_h1.GetNbinsX()):
        i_y = cdf_h1.GetBinContent(ibin)
        if y_val-tolerance <= i_y < y_val+tolerance:
            x_val = cdf_h1.GetBinCenter(ibin)
    return round(x_val,1)

base_plot_folder = "/exp/atlas/kurdysh/vbs_cross_terms_study/plots/"
dict_bounds = {}
for ana in plot_folder.keys(): dict_bounds[ana]={}

for ana in list(h_name.keys()):
    plots_dir = f"{base_plot_folder}/{ana}/{plot_folder[ana]}/{ana}_clipping_bounds/"
    plots_dir_diboson, plots_dir_integ = plots_dir+"/diboson_hist/", plots_dir+"/diboson_hist_integral/"
    _, top_files_dir = lu.find_prod_dec_and_dir(f"user.okurdysh.MadGraph_{ana}_FM0_SM")
    for op_dir in [i_obj for i_obj in os.listdir(top_files_dir) if os.path.isdir(top_files_dir + "/" + i_obj)]:
        full_op_dir = os.path.join(top_files_dir,op_dir,plot_folder[ana],"")
        _, order, op = lu.get_op_from_dir(op_dir, ana)
        if order != "QUAD": continue
        pdfname = f"{op}_{order}.pdf"
        hfile =  ROOT.TFile.Open(full_op_dir+"hists.root", "READ")
        diboson_hist = hfile.Get(h_name[ana])
        diboson_hist.RebinX(2)
        peak_mass = diboson_hist.GetBinCenter(diboson_hist.GetMaximumBin())
        diboson_hist.SetTitle(diboson_hist.GetTitle()+f"_peak_{peak_mass:.2f}")
        lu.save_plot(diboson_hist,plots_dir_diboson+pdfname, draw_option = "")
        diboson_hist_cumulative = diboson_hist.GetCumulative()
        lu.save_plot(diboson_hist_cumulative,plots_dir_integ+pdfname, draw_option = "")
        mass_5_percent = x_val_for_y(diboson_hist_cumulative, 0.05)
        mass_95_percent = x_val_for_y(diboson_hist_cumulative, 0.95)
        # print("for ana", ana, "op order", op, order, "diboson 5% is is at", mass_5_percent, "95% at", mass_95_percent)
        dict_bounds[ana][op] = [mass_5_percent, mass_95_percent]

all_ops_dublication = []
for i_ana in list(plot_folder.keys()):
    i_ana_ops = list(dict_bounds[i_ana].keys())
    for i_op in i_ana_ops:
        all_ops_dublication.append(i_op)
all_ops = sorted(list(set(all_ops_dublication)))
df_sum = pd.DataFrame(index=all_ops, columns=list(plot_folder.keys()))
for i_ana in plot_folder.keys():
    for i_op in dict_bounds[i_ana].keys():
        b_arr = dict_bounds[i_ana][i_op]
        df_sum.at[i_op, i_ana] = f"{b_arr[0]}, {b_arr[1]}"
lu.save_df(df_sum, f"{base_plot_folder}/clipping_reasonable_bounds.pdf", aspect = (9, 9))
