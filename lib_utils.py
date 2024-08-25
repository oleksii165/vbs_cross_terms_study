import glob
import os
from array import array
import ROOT
import itertools
import matplotlib.pyplot as plt
plt.rcParams.update({'text.usetex': True,
                     "backend": 'PDF'})
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import numpy as np
import math
import json
import re
from pandas.plotting import table

def get_fitted_plot(routine, cut):
    mystr,bins = "",[]
    overflow_bin, exclude_underverflow_from_norm = -1, -1
    if routine == "Zy_vvy":
        mystr = "pt_photon"
        bins = array('d', [150,300,450,600,750,900,1050,1200,2000])
        overflow_bin, exclude_underverflow_from_norm = 1, 0
    elif routine in ["WmZ_lllv","WpZ_lllv","WZ_lllv"]:
        mystr="m_T_WZ"
        bins = array('d', [0,400,750,1050,1350,1800]) # last bin is overflow
        overflow_bin, exclude_underverflow_from_norm = 1, 0 # should be opposute
    elif routine == "ZZ_llvv":
        mystr="pt_Z"
        bins = array('d', [50,100,150,200,250,350,1500])
        # bins = array('d', [50,110,130,150,170,200,250,350,1500])
    elif routine in ["WmWm_lvlv","WpWp_lllv","ssWW_lvlv","ATLAS_2023_I2729396"]:
        if cut=="SR":
            mystr, bins = "m_ll", array('d', [0,250,500,750,1000,1500])
        elif cut=="LowmjjCR":
            mystr, bins = "m_ll", array('d', [0,1500])
        elif cut=="WZCR":
            mystr, bins = "m_tagjets", array('d', [0,3000])
    elif routine in ["Wmy_lvy","Wpy_lvy","Wy_lvy", "Wy_John"]:
        mystr = "pt_lepton"
        bins = array('d', [30,43,60,85,130,550])
        overflow_bin, exclude_underverflow_from_norm = 0, 1
    return mystr, bins, overflow_bin, exclude_underverflow_from_norm

def get_rivet_job_name(genJobName,routine,cut):
    return f'{genJobName}_rivet_{routine}_cut_{cut}'

def get_missing_ops(routine):
    if routine=="Zy_vvy":
        ops = ["FS0", "FS1", "FS2", "FS02",
                       "FM7",
                       "FM3", "FM4", "FM5",
                       "FT1", "FT2",
                       "FT6", "FT7"]
    elif routine in ["WpZ_lllv","WmZ_lllv", "WZ_lllv"]:
        ops = ["FM2", "FM3", "FM4", "FM5",
               "FT5", "FT6", "FT7",
               "FT8", "FT9"]
    elif routine in ["ZZ_llvv"]:
        ops = ["FS02", "FS1", 
                "FM0", "FM1", "FM2", "FM3", "FM4", "FM5", "FM7"]
    elif routine in ["WmWm_lvlv","WpWp_lvlv","ssWW_lvlv"]:
        ops = ["FM2", "FM3", "FM4", "FM5", "FT5", "FT6", "FT7", "FT8", "FT9",
            #    "FT2"
               ] # not missing but bad agreement with ws at clip inf for quad
    else:
        ops=[]
    return ops

def get_var_latex(varname): # to be used with ROOT latex
    latexstr=varname
    if varname=="pt_photon": latexstr="p_{T}^{\gamma} [GeV]"
    elif varname=="m_WZ_T": latexstr="m_{T}^{WZ} [GeV]"
    return latexstr

def latex_ana_str(prod_dec):
    mystr=""
    if "Zy_lly" in prod_dec: mystr = "Z(" + r"$\rightarrow$" + f"ll)y"
    return mystr

def get_cutflow_arrays(cutflow_file):
    names = []
    cumu = []
    incr = []
    with open(cutflow_file) as f:
        for i_line in f.readlines():
            line_arr = [i_thing.strip().replace("%", "") for i_thing in i_line.split(" ") if i_thing != ""]
            print(line_arr)
            if len(line_arr) != 5: continue
            names.append(line_arr[1])
            cumu.append(float(line_arr[3]))
            incr.append(float(line_arr[4]))
    return names, cumu, incr

def draw_cutflows(cut_names, y_arrays, labels_arr, outname, prod_dec=""):
    cycle_colors = itertools.cycle('bgrcmk')
    plt.clf()
    fig, ax = plt.subplots(figsize=(6, 4))
    for y_arr,label in zip(y_arrays, labels_arr):
        plt.plot(range(len(cut_names)), y_arr, marker='o', label=label, color=next(cycle_colors))
    plt.xticks(range(len(cut_names)), cut_names, fontsize=7)
    plt.xticks(rotation=90)
    for xc in cut_names: ax.axvline(x=xc, color='0.3', linestyle='--', linewidth=0.3)
    plt.ylim(0, 110)
    plt.legend()
    plt.ylabel('fraction $[\%]$')
    plt.xlabel('cut')
    op = outname[outname.find("_F")+1 : outname.find("_EXT0")]
    plt.title(latex_ana_str(prod_dec)+" "+op)
    plt.savefig(outname, bbox_inches='tight')

# @TODO maybe dublication with find_evnt_dir_and_file
def get_evnt_log_files(base_dir,i_job_name):
    evnt_did, evnt_dir, log_did, log_dir_before_untar = get_envt_log_names_dirs(base_dir,i_job_name)
    evnt_candidates_out = -1
    log_file = -1
    if os.path.exists(evnt_dir):
        print("directories for evnt  exist")
        evnt_candidates = glob.glob(evnt_dir + "/*EVNT.root")
        print("evnt candidates of len",len(evnt_candidates), evnt_candidates)
        if len(evnt_candidates)>0:
            evnt_candidates_out = evnt_candidates
    else:
        print("directories for evnt  DOESN:T exist")

    if os.path.exists(log_dir_before_untar):
        print("directories for log exist")
        log_candidates = glob.glob(log_dir_before_untar + "/tarball_PandaJob*/log.generate")
        print("log candidates of len",len(log_candidates), log_candidates)
        if len(log_candidates)==1:
            log_file =  log_candidates[0]
    else:
        print("directories for log  DOESN:T exist")

    print("returning evnt files", evnt_candidates_out)
    print("returning log file", log_file)
    return evnt_candidates_out, log_file

def get_plotdir(prod_dec, routine, cut, DOCUT_format=False):
    routine_dir = f"routine_{routine}_cut_{cut}" if not DOCUT_format else f"routine_{routine}_cut_{cut}"
    my_dir = f"/exp/atlas/kurdysh/vbs_cross_terms_study/plots/{prod_dec}/{routine_dir}/"
    if not os.path.exists(my_dir): os.makedirs(my_dir)
    return my_dir, routine_dir

def get_pair_str(op1,op2):
    mypair = sorted([op1,op2])
    return f"{mypair[0]}vs{mypair[1]}"

def get_bookletdir(start_path, normalized=""):
    my_dir = start_path + "/booklets/"
    if len(normalized)>0: my_dir = my_dir[:-1] + "_normalized/"
    if not os.path.exists(my_dir): os.makedirs(my_dir)
    if not os.path.exists(my_dir + "/svg/"): os.makedirs(my_dir + "/svg/")
    return my_dir

def get_ops(include_fs0_2):
    all_ops = ["FM0","FM1","FM2","FM3","FM4","FM5","FM7",
                "FS02","FS1",
                "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]
    if include_fs0_2: 
        all_ops += ["FS0","FS2"]
        all_ops.remove("FS02")
    op_pairs = sorted(list(itertools.combinations(all_ops,2)))
    return all_ops, op_pairs

def get_hists_bounds_cuts(prod_dec):
    with open(f"{prod_dec}_hists.json") as fo: total_h = json.load(fo)
    # always load analysis specific hists and jet ones
    with open("jet_hists.json") as fo: jet_h = json.load(fo)
    with open("photon_hists.json") as fo: photon_h = json.load(fo)
    with open("lepton_hists.json") as fo: lepton_h = json.load(fo)
    total_h.update(jet_h)
    # update others depending on ana
    if prod_dec in ["Zy_lly", "Wmy_lvy"]:
        total_h.update(lepton_h)
        total_h.update(photon_h)
    elif prod_dec=="Zy_vvy":
        total_h.update(photon_h)
    if prod_dec in ["WmZ_lllv", "WpZ_lllv"]:
        total_h.update(lepton_h)
    
    # this loop is unncessecary
    return_dict = {} # key is hist name, value[nbin,min,max,cut,cutdir="+,-"]
    for i_hist, i_h_arr in total_h.items():
        if len(i_h_arr)!=3: 
            print("----were not able to  for",i_hist,"will not be there in kin plots")
            continue
        return_dict[i_hist]=i_h_arr
    # replace original bin params with what want to have [rebin_x, x_low, x_up]
    params={}
    params["pt_tagjet1"]  = [10, 0, 2500]
    params["pt_tagjet2"] = [10, 0, 1000]
    params["eta_tagjets"] = [3, -1, -1]
    params["phi_tagjets"] = [5, 0, 6]
    params["m_tagjets"] = [10, -1, -1]
    params["dy_tagjets"]= [5, 0, 9]
    for i_dphi in ["dphi_tagjets","dphi_MET_photon", "dphi_MET_tagjet"]:
        params[i_dphi]= [-1, 0, 3.5]
    params["pt_lepton"]  = [10, 0, 2500]
    params["eta_lepton"] = [3, -3, 3]
    # params["pt_photon"]  = params["pt_MET"] = [10, 1200, 4000]
    params["pt_MET"] = [10, 1200, 4000]
    params["eta_photon"] = [3, -3, 3]
    params["m_ll"] = [10,0,300]
    params["m_lly"] =  params["m_ly"] = [10,-1,1]
    for i_cent in ["centrality_lly","centrality_jjy","centrality_jjly"]:
        params[i_cent] = [-1,0,0.8]
    params["cone_frac_photon"] = [-1, 0, 0.1]
    params["m_W_T"] = [3, 0, 150]
    # params["m_WZ_T"] = [10, 0, 4000]
    params["dR_lepton_photon"] = [2, 0, 6]
    params["dR_tagjets"] = [2,0,10]
    default_params = [-1,-1,-1]
    for i_hist, i_h_arr in return_dict.items():
        if i_hist in params.keys(): replace_params=params[i_hist]
        else: replace_params = default_params
        i_h_arr[0]=replace_params[0]
        i_h_arr[1]=replace_params[1]
        i_h_arr[2]=replace_params[2]

    print("returning plotting dict")
    return return_dict 

def find_prod_dec_and_dir(conf):
    prod_temp = conf[conf.find(".MadGraph_")+len(".MadGraph_") : ]
    print("start from string", prod_temp)
    prod_dec = prod_temp[:prod_temp.find("_F")]
    print("from conf found production dec", prod_dec)
    conf_dir = f"/sps/atlas/kurdysh/vbs_cross_terms_study/eft_files/{prod_dec}/"
    print("dir would be", conf_dir)
    return prod_dec, conf_dir

def find_evnt_dir_and_file(search_com):
    conf_dir_arr = glob.glob(search_com)
    print("found possibilities for dir", conf_dir_arr)
    conf_dir = conf_dir_arr[0] if len(conf_dir_arr)>=1 else -1  
    if conf_dir == -1: raise ValueError("did not find folder for this config ", search_com)

    evnt_file_candidates = glob.glob(conf_dir + "/*EVNT.root")
    if len(evnt_file_candidates)==0: evnt_file_candidates=-1

    return conf_dir, evnt_file_candidates

def get_conf_cut_dir(evnt_dir, routine, cut):
    mydir = evnt_dir + f"/routine_{routine}_cut_{cut}/"
    if not os.path.exists(mydir): os.makedirs(mydir) 
    return mydir

def find_last_match_job(names_arr, str_op_EFTmode):
    matches_temp = []
    for i_name in names_arr:
        if str_op_EFTmode in i_name: matches_temp.append(i_name)
    matches = sorted(matches_temp)
    print("found possible jobs", matches)
    last_try_job = matches[-1] if len(matches)>0 else -1 
    print("select the last job ", last_try_job)
    return last_try_job

def get_envt_log_names_dirs(base_dir,i_job_name):
    evnt_did = i_job_name + "_EXT0"
    evnt_dir = base_dir + "/" + evnt_did 
    log_did = i_job_name + ".log"
    log_dir = evnt_dir + "/" + log_did
    print("returnning names structure for envt did and dir",evnt_did, evnt_dir)
    print("returnning names structure log did and dir",log_did, log_dir)
    return evnt_did, evnt_dir, log_did, log_dir  

def get_xsec(log_file, prod_xsec=False):
    with open(log_file) as textf:
        xsec_val, xsec_unit = -999.0 , "fb" # here pb but later for plots will convert to fb
        for line in textf:
            if prod_xsec:
                if 'Cross-section :' in line:
                    line_vec = line.split(" ")
                    xsec_val = float(line_vec[10])
                    xsec_unit = line_vec[13].strip()
            else:
                if 'MetaData: cross-section' in line:
                    xsec_val = float(line[line.find('=')+1:])
                    xsec_unit = line[line.find('(')+1:line.find(')')]
    conv_fact_to_pb = {"mb":1e9, "um":1e6, "nb":1e3, "pb":1, "fb":1e-3}
    xsec_fb = xsec_val * conv_fact_to_pb[xsec_unit] * 1000
    # print("found xsec value ",xsec_val,"with unit",xsec_unit,"converting to fb get in fb",xsec_fb)
    return xsec_fb

def write_to_f(product_file,product):
    f = open(product_file, "w") #since opening yoda a bit slow for 50 configs save when have it
    f.write(str(product))
    f.close()

def yoda_to_root_1d(h_yoda, yoda_title):
    print("converting to root hist with name", yoda_title)
    root_title =  yoda_title
    mjj_h_root = ROOT.TH1D(root_title, root_title, h_yoda.numBins(), array('d', h_yoda.xEdges()))
    mjj_h_root.Sumw2()
    mjj_h_root.GetXaxis().SetTitle(yoda_title)
    rtErrs = mjj_h_root.GetSumw2()
    for i in range(mjj_h_root.GetNbinsX()):
        mjj_h_root.SetBinContent(i + 1, h_yoda.bin(i).sumW())
        rtErrs.AddAt(h_yoda.bin(i).sumW2(), i+1)
    mjj_h_root.SetDirectory(0)
    return mjj_h_root.Clone()

def get_op_from_dir(mydir,prod_dec):
    temp1 = mydir[len(f"user.okurdysh.MadGraph_{prod_dec}_"):]
    temp2 = temp1[:temp1.find("_EXT0")]
    if "try" in temp2: temp3 = temp2[:temp2.find("_try")]
    else: temp3  = temp2
    arr = temp3.split("_")
    # print("vec for this str", arr)
    ops = arr[0]
    regime = arr[1]
    ops_arr = []
    if "vs" in ops:
        ops_arr.append(ops[:ops.find("vs")])
        ops_arr.append(ops[ops.find("vs")+2:])
    else:
        ops_arr.append(ops)
    return sorted(ops_arr), regime, ops

def save_plot(plot,path_to_save, draw_option = "text45", log_scale = False, legend=False):
    if not os.path.exists(os.path.dirname(path_to_save)):
        os.makedirs(os.path.dirname(path_to_save))
    c=ROOT.TCanvas()
    plot.Draw(draw_option)
    if log_scale: ROOT.gPad.SetLogy()
    if legend: 
        l=ROOT.gPad.BuildLegend()
        l.SetFillColorAlpha(0, 0) # transparent label
    c.Modified()
    c.Update()
    c.Show()
    c.SaveAs(path_to_save)

def read_hists(root_file_name, h_names_arr):
    hists = {}
    if os.path.exists(root_file_name):
        h_file =  ROOT.TFile.Open(root_file_name, "READ")
        hists_in_file = [i_obj.GetName() for i_obj in h_file.GetListOfKeys()]
        # print("list of hists in this file", hists_in_file, "in file ", root_file_name)
        for i_name in h_names_arr:
            if i_name in hists_in_file: 
                i_hist = h_file.Get(i_name)
                i_hist.SetDirectory(0)
                hists[i_name] = i_hist 
            else:
                print("there is no hist", i_name, "in file ", root_file_name)
        h_file.Close()
    return hists

def dress_hist(my_hist_in, my_title, my_color=1, my_norm = 1.0, re_bins=-1, re_overflow=0, exclude_underverflow_from_norm=False):
    if exclude_underverflow_from_norm and re_overflow:
        raise ValueError("doesnt make sense when both exclude_underverflow_from_norm and re_overflow")
    my_hist = my_hist_in.Clone() # otherwise will modify name of original hist etc
    my_hist.SetTitle(my_title)
    my_hist.SetName(my_title)
    my_hist.SetLineColor(my_color)
    my_hist.SetMarkerColor(my_color)
    if re_bins!=-1 and re_overflow==0:
        my_hist = my_hist.Rebin(len(re_bins) - 1, my_hist.GetName(), re_bins)
        print("hi")
    elif re_bins!=-1 and re_overflow:
        orig_nbins = my_hist.GetNbinsX()
        orig_end_x = my_hist.GetBinLowEdge(orig_nbins) + my_hist.GetBinWidth(orig_nbins)
        orig_hist_rebin_with_overflow = my_hist.Clone()
        re_bins_overflow = array('d',re_bins.tolist().copy())
        re_bins_overflow[-1] = orig_end_x
        orig_hist_rebin_with_overflow = orig_hist_rebin_with_overflow.Rebin(len(re_bins_overflow)-1, "over", re_bins_overflow)
        my_hist = my_hist.Rebin(len(re_bins) - 1, my_hist.GetName(), re_bins)
        rebin_nbins = my_hist.GetNbinsX()
        # keep visually same binning but include everything in last bin
        my_hist.SetBinContent(rebin_nbins,orig_hist_rebin_with_overflow.GetBinContent(rebin_nbins))

    hist_integ = my_hist.Integral()
    if exclude_underverflow_from_norm:
        under_c = my_hist.GetBinContent(0)
        over_c = my_hist.GetBinContent(my_hist.GetNbinsX()+1)
        integ_with_under_over = hist_integ + over_c + under_c
        if integ_with_under_over!=0:
            my_norm_corr = (hist_integ / integ_with_under_over) * my_norm
        else:
            my_norm_corr = 0
    else:
        my_norm_corr = my_norm
    if hist_integ!=0:
        my_hist.Scale(my_norm_corr/hist_integ)

    return my_hist.Clone()

def make_stack(hist_arr, title="",norm=-1):
    my_stack = ROOT.THStack(f"{title}/{norm}", f"{title}/{norm}")
    for i_plot in hist_arr:
        plot_copy = i_plot.Clone() # since will use several time with/without normalzietions  
        if norm!=-1 and plot_copy.Integral()!=0: plot_copy.Scale(norm/plot_copy.Integral())
        my_stack.Add(plot_copy)
    return my_stack.Clone()

def get_ratio_plot_tests(hist_1_in, hist_2_in): # _2 is the one with to respect to which hist_1 is taken
    hist_1 = hist_1_in
    hist_2 = hist_2_in
    ratio_hist_for_gr_range = ROOT.TH1F("dummy", "dummy", hist_1.GetNbinsX(),
                                        hist_1.GetXaxis().GetXmin(),hist_1.GetXaxis().GetXmax())
    for binx in range(1, ratio_hist_for_gr_range.GetNbinsX() + 1):
        i_1 = hist_1.GetBinContent(binx)
        i_2 = hist_2.GetBinContent(binx)
        if i_1 <= 0 or i_2 <= 0: continue
        err_1 = hist_1.GetBinError(binx)
        err_2 = hist_2.GetBinError(binx)
        ratio = i_1 / i_2
        rel_err_1 = err_1 / i_1
        rel_err_2 = err_2 / i_2
        uncert = abs(ratio) * math.sqrt(rel_err_1 ** 2 + rel_err_2 ** 2)
        ratio_hist_for_gr_range.SetBinContent(binx, ratio)
        ratio_hist_for_gr_range.SetBinError(binx, uncert)
    rchi2 = hist_1.Chi2Test(hist_2, "WW CHI2/NDF")
    ks = hist_2.KolmogorovTest(hist_1)
    ratio_name = f"{hist_1_in.GetName()}/{hist_2_in.GetName()} rChi2={rchi2:.2f} KS={ks:.2f}"
    ratio_hist_for_gr_range.SetName(ratio_name)
    ratio_hist_for_gr_range.SetTitle(ratio_name)
    ratio_color = hist_1_in.GetLineColor() # notice it's hist_1 where color is taken
    ratio_hist_for_gr_range.SetMarkerColor(ratio_color)
    ratio_hist_for_gr_range.SetLineColor(ratio_color)
    return ratio_hist_for_gr_range.Clone(), rchi2, ks

def draw_stack_with_ratio(my_stack, mg_ratios, xtitle, outname, stack_x_range=[]):
    c = ROOT.TCanvas()
    y_divide = 0.4 if mg_ratios!=-1 else 0.0
    pad_1 = ROOT.TPad("pad1","pad1", 0.0, y_divide, 1.0, 1.0)
    pad_2 = ROOT.TPad("pad1","pad1", 0.0, 0.0, 1.0, y_divide)
    pad_1.Draw()
    pad_2.Draw()

    pad_1.cd()
    my_stack.Draw("nostack")
    if len(stack_x_range)!=0: 
        my_stack.GetXaxis().SetRangeUser(stack_x_range[0], stack_x_range[1])
    ROOT.gPad.SetLogy()
    my_stack.GetXaxis().SetTitle(xtitle)
    l = ROOT.gPad.BuildLegend()
    l.SetFillColorAlpha(0, 0)

    if mg_ratios!=-1:
        pad_2.cd()
        mg_ratios.Draw("")
        if len(stack_x_range) != 0:
            mg_ratios.GetXaxis().SetRangeUser(stack_x_range[0], stack_x_range[1])
        l = ROOT.gPad.BuildLegend()
        l.SetFillColorAlpha(0, 0)
        ROOT.gPad.SetLogy()

    c.Modified()
    c.Update()
    c.Show()
    c.SaveAs(outname)


def get_replacement(df, op_miss, list_to_search_in):
    column = df[op_miss]
    def do_search(sublist_to_search_in):
        best_chi2 = 1000000
        best_rep = ""
        for rep, chi2 in column.items():
            if rep not in sublist_to_search_in: continue
            if rep==op_miss: continue
            if chi2 < best_chi2:
                best_chi2 = chi2
                best_rep = rep
        return best_rep, round(best_chi2,2)

    op_miss_family = op_miss[:re.search(r"\d", op_miss).start()]
    # first seach within same family
    same_family_search_list = [i_op for i_op in  list_to_search_in if op_miss_family in i_op]
    rep_sf, chi2_sf = do_search(same_family_search_list)
    # then inside all of them
    rep_all, chi2_all = do_search(list_to_search_in)
    if abs(chi2_sf-chi2_all) < 1:
        best_rep, best_chi2 = rep_sf, chi2_sf
    else:
        best_rep, best_chi2 = rep_all, chi2_all
    return best_rep, best_chi2


def save_df(df, out_path, save_csv=False, aspect = (16, 9), save_pdf=True, save_latex=False):
    if save_csv: df.to_csv(out_path+".csv", sep=";")
    if save_latex:
        fname_no_dir_ext = os.path.basename(out_path)
        label = fname_no_dir_ext[:fname_no_dir_ext.find(".pdf")]
        df.to_latex(out_path+".tex", label="table:"+label)
    if save_pdf:
        plt.clf()
        fig, ax = plt.subplots(figsize=aspect)
        ax.axis('tight')
        ax.axis('off')
        table_quad_pairs = ax.table(cellText = df.values, rowLabels = df.index, colLabels = df.columns, loc='center')
        with PdfPages(out_path) as pdf:
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()

def arr_to_hist1d(arr_counts, hname, bin_edges):
    my_hist = ROOT.TH1F(hname, hname, len(bin_edges)-1, bin_edges)
    for i_num_bin, i_bin_count in zip(range(1,my_hist.GetNbinsX()+1),arr_counts):
        my_hist.SetBinContent(i_num_bin, i_bin_count)
    return my_hist.Clone()

def get_rivet_resub_name(last_job_name):
    rivet_str_ind = last_job_name.find("rivet")
    rivet_part_name = last_job_name[rivet_str_ind:]
    if "try" not in rivet_part_name:
        job_name = rivet_part_name + "_try2"
    else:
        last_job_name_try_part = rivet_part_name[rivet_part_name.find("try"):]
        last_job_name_try_attempt = int(last_job_name_try_part[3:])
        job_name = rivet_part_name.replace(last_job_name_try_part,"try"+str(last_job_name_try_attempt+1))
    full_job_name = last_job_name[:rivet_str_ind] + job_name
    return full_job_name

def ssWW_get_ws_hist(ws_json_clip,op,order_capital, i_hist_out_name,
                gen_prod_dec, cut):
    substr_cut_rivet_to_ws_dict = {"SR": "CutSR", "LowmjjCR": "CutLowMjj", "WZCR": "CutCRWZ3Lep"}
    eft_substr = "EFTWZ_" if gen_prod_dec=="WZ_lllv" else "EFT_"
    op_substr = op[1:] if "F" in op else op
    arr_counts = -1
    # to select certain region it's "cutSR"
    index_region_dict= {"SR": 3, "WZCR": 0, "LowmjjCR": 2}
    sr_samples = ws_json_clip['distributions'][index_region_dict[cut]]['samples']
    for i_hist_json in sr_samples:
        i_name = i_hist_json['name'] # like 'EFTWZ_M1_28_quad_CutSR_all_overallSyst_x_StatUncert'
        if eft_substr not in i_name: continue
        if op_substr not in i_name: continue
        if order_capital.lower() not in i_name: continue
        if substr_cut_rivet_to_ws_dict[cut] not in i_name: continue
        arr_counts = i_hist_json['data']['contents']
    _, fit_plot_bins, _, _ = get_fitted_plot("ssWW_lvlv", cut)
    i_hist = arr_to_hist1d(arr_counts, i_hist_out_name, fit_plot_bins) if arr_counts!=-1 else -1
    return i_hist, arr_counts

def ssWW_get_wilson_coef(ops_arr, gen_prod_dec, order):
    # final wilson = 8 * their wilson. table gives those numbers. theis is brough up for some reason
    wilson_corr_coef = {}
    wilson_corr_coef["ssWW_lvlv"] = {"FS02":8.0,"FS1":20.0,"FM0":6.0,"FM1":10.0,"FM7":13.0,"FT0":0.6,"FT1":0.3,"FT2":1.0}
    wilson_corr_coef["WZ_lllv"] = {"FS02":82,"FS1":128,"FM0":27,"FM1":28,"FM7":30,"FT0":2.4,"FT1":1.6,"FT2":5.5}
    mult_factor_on_my_xsec = -1
    op_search = ops_arr[0]
    if op_search in wilson_corr_coef[gen_prod_dec].keys():
        coef = wilson_corr_coef[gen_prod_dec][op_search]
        if order == "QUAD":
            mult_factor_on_my_xsec = coef ** 2
        elif order == "INT":
            mult_factor_on_my_xsec = coef
    return mult_factor_on_my_xsec

def get_clip_hist_name(fitvar,clip):
    return f"{fitvar}_clip_{clip}"

def read_oneline_file(xsec_file):
    with open(xsec_file, 'r') as f:
        fid_xsec_fb = float(f.read())
    return fid_xsec_fb

    