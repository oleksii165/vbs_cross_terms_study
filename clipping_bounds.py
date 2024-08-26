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
          "WpZWmZ_lllv": "m_WZ",
          "ZZ_llll": "m_ZZ"}
h_name_latex = {"WpyWmy_lvy": "m_{W\gamma}",
                "Zy_vvy": "m_{Z\gamma}",
                "WpZWmZ_lllv": "m_{WZ}",
                "ZZ_llll": "m_{ZZ}"}
plot_folder = {"WpyWmy_lvy": "routine_Wy_John_cut_SR",
               "Zy_vvy": "routine_Zy_vvy_cut_SR",
               "WpZWmZ_lllv": "routine_WZ_lllv_cut_SR",
               "ZZ_llll": "routine_ZZ_llll_cut_SR"}
rebin_x = {"WpyWmy_lvy": 6,
           "Zy_vvy": 6,
           "WpZWmZ_lllv": 6,
           "ZZ_llll": 6}


base_plot_folder = "/exp/atlas/kurdysh/vbs_cross_terms_study/plots/"
dict_frac_at_low_clip = {}
check_frac_at_clip = 1000
for ana in plot_folder.keys(): dict_frac_at_low_clip[ana]={}

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
        # if "FT" not in op: continue
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
        ih_QUAD_INT = ih_QUAD.Clone()
        ih_QUAD_INT.Add(ih_INT)

        icdf_QUAD_INT = mho.make_abs_hist(ih_QUAD_INT).GetCumulative()
        if icdf_QUAD_INT.Integral() != 0:
            icdf_QUAD_INT.Scale(1/icdf_QUAD_INT.Integral())

        pdfname = f"{i_op}_INT_QUAD.pdf"
        mrd.one_plot([ih_QUAD, ih_INT, ih_QUAD_INT],
                     plots_dir_diboson+pdfname, make_out_missing_dirs = True,
                     draw_opt = "",
                     x_label = h_name_latex[ana], y_label = "",
                     legends_arr = [f"{i_op} QUAD", f"{i_op} INT", "QUAD+INT"], legends_type_arr_in = ["p","p", "p"],
                     canvas_aspect=[4,3],
                     use_colours_1d = True,
                     custom_x_range = [0,5000],
                     custom_y_range = [min([ih_QUAD.GetMinimum(), ih_INT.GetMinimum(), ih_QUAD_INT.GetMinimum()]),
                                       max([ih_QUAD.GetMaximum(), ih_INT.GetMaximum(), ih_QUAD_INT.GetMaximum()])])
        mrd.one_plot([icdf_QUAD_INT],
                     plots_dir_integ+pdfname, make_out_missing_dirs = True,
                     draw_opt = "",
                     x_label = h_name_latex[ana], y_label = "|Cumulative Distribution Function|",
                     legends_arr = [f"{i_op} QUAD+INT"], legends_type_arr_in = ["p"],
                     canvas_aspect=[4,3],
                     use_colours_1d = True,
                     custom_x_range = [0,5000])
        frac_at_clip = icdf_QUAD_INT.Integral(1,icdf_QUAD_INT.FindBin(check_frac_at_clip))
        dict_frac_at_low_clip[ana][i_op] = '{:.2f}'.format(float(frac_at_clip*100), 2)

all_ops_dublication = []
for i_ana in list(plot_folder.keys()):
    i_ana_ops = list(dict_frac_at_low_clip[i_ana].keys())
    for i_op in i_ana_ops:
        all_ops_dublication.append(i_op)
all_ops = sorted(list(set(all_ops_dublication)))
df_sum = pd.DataFrame(index=all_ops, columns=list(plot_folder.keys()))
for i_ana in plot_folder.keys():
    for i_op in dict_frac_at_low_clip[i_ana].keys():
        df_sum.at[i_op, i_ana] = dict_frac_at_low_clip[i_ana][i_op]
lu.save_df(df_sum, f"{base_plot_folder}/frac_INT+QUAD_at_clipping_{check_frac_at_clip}.pdf", aspect = (9, 9), save_latex=True)
