import ROOT
import lib_utils as lu
import os
import glob

# basedir="/sps/atlas/kurdysh/vbs_cross_terms_study/eft_files/Zy_vvy/"
# op_dirs = {"FT0": "user.okurdysh.MadGraph_Zy_vvy_FT0_QUAD_try2_EXT0/DOCUT_YES/",
# "FT2": "user.okurdysh.MadGraph_Zy_vvy_FT2_QUAD_EXT0/DOCUT_YES/", "FT5": "user.okurdysh.MadGraph_Zy_vvy_FT5_QUAD_try2_EXT0/DOCUT_YES/",
# "FM2": "user.okurdysh.MadGraph_Zy_vvy_FM2_QUAD_EXT0/DOCUT_YES/", "FM5": "user.okurdysh.MadGraph_Zy_vvy_FM5_QUAD_EXT0/DOCUT_YES/",
# }
# h_stack = ROOT.THStack()
# for num_hist, i_op in enumerate(op_dirs.keys(),start=1):
#     i_root_file_name  = basedir + op_dirs[i_op] + "hists.root"
#     i_h_file =  ROOT.TFile.Open(i_root_file_name, "READ")
#     i_hist = i_h_file.Get("m_Zy").Clone()
#     i_hist.SetDirectory(0)
#     i_hist.RebinX(4)
#     i_hist.SetName(i_op);i_hist.SetTitle(i_op)
#     i_hist.SetLineColor(num_hist);i_hist.SetMarkerColor(num_hist)
#     h_stack.Add(i_hist)
# plot_path  = "../plots/Zy_vvy/m_Zy_comaprison_no_rwgt.pdf"
# lu.save_plot(h_stack,plot_path, draw_option = "nostack", log_scale = False, legend=True)

# test MG reweighting - kinematics
# dir_rw = "/sps/atlas/kurdysh/vbs_cross_terms_study/eft_files/Zy_vvy/user.okurdysh.MadGraphRw100kNoHelicityChangeLO_Zy_vvy_FT0_QUAD_EXT0/DOCUT_YES/"
# dir_m2 = "/sps/atlas/kurdysh/vbs_cross_terms_study/eft_files/Zy_vvy/user.okurdysh.MadGraph_Zy_vvy_FM2_QUAD_EXT0/DOCUT_YES/"
# dir_t0 = "/sps/atlas/kurdysh/vbs_cross_terms_study/eft_files/Zy_vvy/user.okurdysh.MadGraph_Zy_vvy_FT0_QUAD_EXT0/DOCUT_YES/"
# hists = ["dphi_MET_photon","dphi_MET_tagjet","dphi_tagjets","dy_tagjets","eta_photon","eta_tagjets","m_Zy","m_tagjets","n_jets","n_lepton_stable","n_photons_iso","phi_tagjets","pt_MET","pt_photon_clip_inf","pt_tagjet1","pt_tagjet2"]
# plot_dir  = "../plots/Zy_vvy/plot_rw_m2/"
# if not os.path.exists(plot_dir): os.makedirs(plot_dir)
# rfile_rw = ROOT.TFile(dir_rw + "hists_rwgt_M2_QUAD.root","read") 
# rfile_nom = ROOT.TFile(dir_m2 + "hists.root","read")
# for ihname in hists:
#     name_nom = ihname
#     name_rw = ihname + "_rwgt_M2_QUAD_"
#     h_rw = rfile_rw.Get(name_rw)  
#     h_rw.SetLineColor(2);h_rw.SetMarkerColor(2) 
#     h_nom = rfile_nom.Get(name_nom)
#     istack = ROOT.THStack()
#     istack.Add(h_rw)
#     istack.Add(h_nom)
#     lu.save_plot(istack,plot_dir+f"{ihname}.pdf", draw_option = "nostack", log_scale = False, legend=True)

# check xsec of different ops - can merge T012? --- no can't but maybe good later to identify blocks
# ana_files_dir = "/sps/atlas/kurdysh/vbs_cross_terms_study/eft_files/WmZ_lllv/"
# job_folders = os.listdir(ana_files_dir)
# for i_job in job_folders:
#     i_prod_dec, _ = lu.find_prod_dec_and_dir(i_job[:i_job.find("_EXT0")])
#     _, _, i_op  = lu.get_op_from_dir(i_job,i_prod_dec)
#     i_log_files = glob.glob(f"{ana_files_dir}{i_job}/*.log/tarball*/log.generate")
#     if len(i_log_files)==0: continue
#     xsec_fb =  lu.get_xsec(i_log_files[0])


# check xsec for two dynscaled in zvvy
def get_xsec(my_dir,op_pattern):
    xsec_fb = -1
    i_log_files = glob.glob(f"{my_dir}/*{op_pattern}*/tarball*/log.generate")
    # print(i_log_files)
    if len(i_log_files)==1:
    	xsec_fb =  lu.get_xsec(i_log_files[0])
    return xsec_fb

# ana_files_orig_dyn_3 = "/sps/atlas/kurdysh/vbs_cross_terms_study/eft_files/Zy_vvy_dyntest/"
# ana_files_semi_dyn = "/sps/atlas/kurdysh/vbs_cross_terms_study/eft_files/Zy_vvy/"

# for i_op_pattern in ["FT0_QUAD","FT1_QUAD","FT2_QUAD","FT5_QUAD","FT6_QUAD","FT7_QUAD","FT8_QUAD","FT9_QUAD","FM0_QUAD","FM1_QUAD","FM2_QUAD","FM3_QUAD","FM4_QUAD","FM5_QUAD","FM7_QUAD"]:
#  dyn3 =   get_xsec(ana_files_orig_dyn_3,i_op_pattern)	
#  dynsemi = get_xsec(ana_files_semi_dyn,i_op_pattern)
#  print(f"{i_op_pattern} semi in fb:",  dynsemi, "and with orin dynscale 3",dyn3, "ratio to dyn3", dynsemi/dyn3)


# check if anything changed  from-regenerating ssWW
# dir_base = "/sps/atlas/kurdysh/vbs_cross_terms_study/eft_files/"
# dir_ssww = dir_base + "ssWW_lvlv"
# dir_wmwm, dir_wpwp = dir_base + "WmWm_lvlv", dir_base + "WpWp_lvlv"
# for i_op_pattern in ["FM0_QUAD","FM1_QUAD","FM7_QUAD","FS02_QUAD","FS1_QUAD","FT0_QUAD","FT1_QUAD","FT2_QUAD"]:
#     xsec_orig = get_xsec(dir_ssww,i_op_pattern)
#     xsec_charges_sep = get_xsec(dir_wmwm,i_op_pattern)+get_xsec(dir_wpwp,i_op_pattern)
#     print(f"{i_op_pattern} re-gen orig in fb:", xsec_orig , "and for my prev wmwm+wpwp", xsec_charges_sep, "ratip orig/prev", xsec_orig/xsec_charges_sep)

# def read_xsec_not_log(my_path):
#     file = open(my_path, "r")
#     content = file.read()
#     file.close()
#     return float(content)

# for i_op_pattern in ["FM0_QUAD","FM1_QUAD","FM7_QUAD","FS02_QUAD","FS1_QUAD","FT0_QUAD","FT1_QUAD","FT2_QUAD"]:
#     path_xsec_orig = glob.glob(f"{dir_ssww}/*{i_op_pattern}*/routine_ATLAS_2023_I2729396_cut_SR/xsec_times_frac_fb_clip_inf.txt")
#     path_xsec_wmwm = glob.glob(f"{dir_wmwm}/*{i_op_pattern}*/routine_ssWW_lvlv_cut_SR/xsec_times_frac_fb_clip_inf.txt")
#     path_xsec_wpwp = glob.glob(f"{dir_wmwm}/*{i_op_pattern}*/routine_ssWW_lvlv_cut_SR/xsec_times_frac_fb_clip_inf.txt")
#     if len(path_xsec_orig)!=1 or len(path_xsec_wmwm)!=1 or len(path_xsec_wpwp)!=1: continue 
#     xsec_orig = read_xsec_not_log(path_xsec_orig[0])
#     xsec_charges_sep = read_xsec_not_log(path_xsec_wmwm[0]) + read_xsec_not_log(path_xsec_wpwp[0])
#     print(f"{i_op_pattern} re-gen orig in fb:", xsec_orig , "and for my prev wmwm+wpwp", xsec_charges_sep, "ratip orig/prev", xsec_orig/xsec_charges_sep)

base_dir = "/sps/atlas/kurdysh/vbs_cross_terms_study/eft_files/test/zzllll_closure_dynscale/"
for op in ["FM0","FT2"]:
    for dyn in [2,3]:
        xsecs = {}
        for order in ["FULL","QUAD","INT","SM"]:
            xsecs[order] = get_xsec(f"{base_dir}/{op}/", f"dyn_{dyn}*{op}_{order}*.log")
        sum_components = xsecs["SM"]+xsecs["INT"]+xsecs["QUAD"]
        print("for op dyn", op, dyn, "SM+INT+QUAD in fb", sum_components , "and FULL", xsecs["FULL"], f"ratio full/sum_components {xsecs['FULL']/sum_components:.2f}", )