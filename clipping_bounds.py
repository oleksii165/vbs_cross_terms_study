# at the moment have clip hists for ssWW Zvvy and ZZllvv
import lib_utils as lu
import os
import myROOTDrawing.utils as mrd
import myROOTDrawing.hist_ops as mho
import ROOT
import pandas as pd
ROOT.gROOT.SetBatch(ROOT.kTRUE)

h_name = {"WpyWmy_lvy": "m_Wy",
          "Zy_vvy": "m_Zy",
          "WpZWmZ_lllv": "m_WZ"}
h_name_latex = {"WpyWmy_lvy": "m_{W\gamma}",
                "Zy_vvy": "m_{Z\gamma}",
                "WpZWmZ_lllv": "m_{WZ}"}
plot_folder = {"WpyWmy_lvy": "routine_Wy_John_cut_SR",
               "Zy_vvy": "routine_Zy_vvy_cut_SR",
               "WpZWmZ_lllv": "routine_WZ_lllv_cut_SR"}
rebin_x = {"WpyWmy_lvy": 6,
           "Zy_vvy": 6,
           "WpZWmZ_lllv": 6}

# def x_val_for_y(cdf_h1, y_val, tolerance=0.01):
#     x_val = -1
#     for ibin in range(cdf_h1.GetNbinsX()):
#         i_y = cdf_h1.GetBinContent(ibin)
#         if y_val-tolerance <= i_y < y_val+tolerance:
#             x_val = cdf_h1.GetBinCenter(ibin)
#     return round(x_val,1)

base_plot_folder = "/exp/atlas/kurdysh/vbs_cross_terms_study/plots/"
dict_m_vv_trans = {}
for ana in plot_folder.keys(): dict_m_vv_trans[ana]={}

for ana in list(h_name.keys()):
    plots_dir = f"{base_plot_folder}/{ana}/{plot_folder[ana]}/{ana}_clipping_bounds/"
    plots_dir_diboson, plots_dir_integ = plots_dir+"/diboson_hist/", plots_dir+"/diboson_hist_integral/"
    _, top_files_dir = lu.find_prod_dec_and_dir(f"user.okurdysh.MadGraph_{ana}_FM0_SM")
    rebin_fact= rebin_x[ana]
    hists_QUAD, hists_INT = {}, {} # key is op
    # just save all the INT QUAD hists available
    for op_dir in [i_obj for i_obj in os.listdir(top_files_dir) if os.path.isdir(top_files_dir + "/" + i_obj)]:
        full_op_dir = os.path.join(top_files_dir, op_dir, plot_folder[ana],"")
        _, order, op = lu.get_op_from_dir(op_dir, ana)
        if order not in ["INT", "QUAD"]: continue
        if "FT" not in op: continue
        hfile = ROOT.TFile.Open(full_op_dir+"hists_norm_run2.root", "READ")
        diboson_hist = hfile.Get(h_name[ana]).Clone()
        diboson_hist.SetDirectory(0)
        diboson_hist.SetName(f"{op}_{order}")
        if order=="QUAD":
            hists_QUAD[op] = diboson_hist
        elif order=="INT":
            hists_INT[op] = diboson_hist
    # now make operations on these ops
    ops_QUAD, ops_INT = list(hists_QUAD.keys()), list(hists_INT.keys())
    ops = list(set(ops_QUAD).intersection(set(ops_INT)))
    for i_op in ops:
        ih_QUAD, ih_INT = hists_QUAD[i_op], hists_INT[i_op]
        ih_QUAD.RebinX(rebin_fact)
        ih_INT.RebinX(rebin_fact)
        icdf_QUAD, icdf_INT = ih_QUAD.GetCumulative(), mho.make_abs_hist(ih_INT).GetCumulative()
        # find where int/quad are equal
        m_vv_trans = -1
        prev_ratio = 10000.0
        for ibin in range(icdf_QUAD.GetNbinsX()):
            cdf_QUAD, cdf_INT = icdf_QUAD.GetBinContent(ibin), icdf_INT.GetBinContent(ibin)
            if cdf_QUAD!=0:
                iratio = cdf_INT/cdf_QUAD
                icenter = icdf_QUAD.GetBinCenter(ibin)
                if prev_ratio>1 and iratio<1:
                    prev_ratio = iratio
                    if icenter > icdf_QUAD.GetBinCenter(5): # to avoid being stuck at first bins
                        m_vv_trans = icenter
        print("+++for ", ana, i_op, "got int/quad=1 at m_vv",m_vv_trans)
        # plot together with lines where they are equal, if any
        pdfname = f"{i_op}_INT_QUAD.pdf"
        standalone_text = f"INT/QUAD=1 at {m_vv_trans} GeV" if m_vv_trans>0 else ""
        mrd.one_plot([ih_QUAD,ih_INT],
                     plots_dir_diboson+pdfname, make_out_missing_dirs = True,
                     draw_opt = "",
                     x_label = h_name_latex[ana], y_label = "",
                     legends_arr = [f"{i_op} QUAD", f"{i_op} INT"], legends_type_arr_in = ["p","p"],
                     canvas_aspect=[4,3],
                     use_colours_1d = True,
                     custom_x_range = [0,4000],
                     custom_y_range = [min([ih_QUAD.GetMinimum(),ih_INT.GetMinimum()]),
                                       max([ih_QUAD.GetMaximum(),ih_INT.GetMaximum()])])
        mrd.one_plot([icdf_QUAD, icdf_INT],
                     plots_dir_integ+pdfname, make_out_missing_dirs = True,
                     draw_opt = "",
                     x_label = h_name_latex[ana], y_label = "|Cumulative Distribution Function|",
                     legends_arr = [f"{i_op} QUAD", f"{i_op} INT"], legends_type_arr_in = ["p","p"],
                     canvas_aspect=[4,3],
                     use_colours_1d = True,
                     custom_x_range = [0,4000],
                     standalone_text = standalone_text)
        dict_m_vv_trans[ana][i_op] = m_vv_trans

all_ops_dublication = []
for i_ana in list(plot_folder.keys()):
    i_ana_ops = list(dict_m_vv_trans[i_ana].keys())
    for i_op in i_ana_ops:
        all_ops_dublication.append(i_op)
all_ops = sorted(list(set(all_ops_dublication)))
df_sum = pd.DataFrame(index=all_ops, columns=list(plot_folder.keys()))
for i_ana in plot_folder.keys():
    for i_op in dict_m_vv_trans[i_ana].keys():
        df_sum.at[i_op, i_ana] = dict_m_vv_trans[i_ana][i_op]
lu.save_df(df_sum, f"{base_plot_folder}/clipping_reasonable_bounds.pdf", aspect = (9, 9), save_latex=True)
