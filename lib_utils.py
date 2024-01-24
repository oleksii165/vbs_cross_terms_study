import glob
import os
from array import array
import ROOT
import itertools

def get_plotdir(prod_dec, DOCUT_str):
    my_dir = f"/exp/atlas/kurdysh/vbs_cross_terms_study/plots/{prod_dec}/{DOCUT_str}/" 
    if not os.path.exists(my_dir): os.makedirs(my_dir)
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

def get_hists_to_draw_with_params(prod_dec):
    hists_dict= {}
    hists_dict["WmWm_lvlv"] = ["pt_tagjet1", "m_tagjets", "deta_tagjets", "lepton_pt","lepton_eta","m_ll", "MET", "m_T"] 
    hists_dict["WpWm_lvlv"] =  ["pt_tagjet1", "m_tagjets", "deta_tagjets", "pt_tagjet2",  "jet3_centralty",  "m_ll", "centrality", "MET"]
    return hists_dict[prod_dec]

def find_prod_dec_and_dir(conf):
    prod_temp = conf[conf.find("user.okurdysh.MadGraph_")+len("user.okurdysh.MadGraph_"):]
    print("start from string", prod_temp)
    prod_dec = prod_temp[:prod_temp.find("_F")]
    print("from conf found production dec", prod_dec)
    conf_dir = f"/exp/atlas/kurdysh/vbs_cross_terms_study/eft_files/{prod_dec}/"
    print("dir would be", conf_dir)
    return prod_dec, conf_dir

def find_evnt_dir_and_file(search_com):
    conf_dir_arr = glob.glob(search_com)
    print("found possibilities for dir", conf_dir_arr)
    conf_dir = conf_dir_arr[0] if len(conf_dir_arr)>=1 else -1  
    if conf_dir == -1: raise ValueError("did not find folder for this config ",search_com)

    evnt_file_candidates = glob.glob(conf_dir + "/*EVNT.root")
    evnt_file = evnt_file_candidates[0] if len(evnt_file_candidates)>0 else -1

    return conf_dir, evnt_file

def get_conf_cut_dir(evnt_dir, docut):
    mydir = evnt_dir + f"/DOCUT_{docut}/"
    if not os.path.exists(mydir): os.makedirs(mydir) 
    return mydir

def get_job_name_template(prod,dec,op,EFTmode):
    return f"user.okurdysh.MadGraph_{prod}_{dec}_{op}_{EFTmode}"

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
    print("returnning envt did and dir",evnt_did, evnt_dir)
    print("returnning log did and dir",log_did, log_dir)
    return evnt_did, evnt_dir, log_did, log_dir  

def get_rivet_com(job_name, evtMax=-1, redoRivet=-1, redoPlots=-1, DOCUT=-1):
    mycom = f'python run_rivet.py --conf="{job_name}" '
    if evtMax!=-1: mycom += f' --evtMax {evtMax} '
    if redoRivet!=-1: mycom += f' --redoRivet "{redoRivet}" '
    if redoPlots!=-1: mycom += f' --redoPlots "{redoPlots}" '
    if DOCUT!=-1: mycom += f' --DOCUT "{DOCUT}" '
    return mycom

def get_xsec(log_file):
    with open(log_file) as textf:
        xsec_val, xsec_unit = -999.0 , "fb" # here pb but later for plots will convert to fb
        for line in textf:
            if 'MetaData: cross-section' in line:
                xsec_val = float(line[line.find('=')+1:])
                xsec_unit = line[line.find('(')+1:line.find(')')]
    conv_fact_to_pb = {"mb":1e9, "um":1e6, "nb":1e3, "pb":1, "fb":1e-3}
    xsec_fb = xsec_val * conv_fact_to_pb[xsec_unit] * 1000
    print("found xsec value ",xsec_val,"with unit",xsec_unit,"converting to fb get in fb",xsec_fb)
    return xsec_fb

def save_xsec_frac_prod(savedir,xsec_fb,frac):
    write_to_f(savedir + "xsec_fb.txt",xsec_fb)
    write_to_f(savedir + "frac_after_cuts.txt",frac)
    write_to_f(savedir + "xsec_times_frac_fb.txt",xsec_fb*frac)

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
    print("vec for this str", arr)
    ops = arr[0]
    regime = arr[1]
    ops_arr = []
    if "vs" in ops:
        ops_arr.append(ops[:ops.find("vs")])
        ops_arr.append(ops[ops.find("vs")+2:])
    else:
        ops_arr.append(ops)
    return sorted(ops_arr), regime

def save_plot(plot,path_to_save, draw_option = "text45", log_scale = False):
    c=ROOT.TCanvas()
    plot.Draw(draw_option)
    if log_scale: ROOT.gPad.SetLogy()
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

def dress_hist(my_hist, my_title, my_color, my_norm = 1.0):
    my_hist.SetTitle(my_title)
    my_hist.SetLineColor(my_color)
    my_hist.SetMarkerColor(my_color)
    hist_integ = my_hist.Integral()
    print("normalize", my_hist.GetName(), my_hist.GetTitle(), "to", my_norm, "start from integ", hist_integ)
    my_hist.Scale(my_norm/hist_integ)
    print("get back integ", my_hist.Integral())
    return my_hist.Clone()




