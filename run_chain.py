# combine panda downloading, running rivet, making plots for many combinations
from pandaclient import panda_api
c = panda_api.get_api()
import subprocess
import os
import glob
import math
import sympy
from sympy import abc
import spb
from array import array
import matplotlib
import yoda
import ROOT
import lib_utils as lu
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)
from optparse import OptionParser

def prepare_grid_files(i_job_name):
    print("will download+untar files evnt and log for", i_job_name)
    evnt_did, evnt_dir, log_did, log_dir = lu.get_envt_log_names_dirs(base_dir,i_job_name)
    if not os.path.exists(evnt_dir):
        print ("will download", evnt_did, "since dir doesn't exist", evnt_dir)
        subprocess.call(f"rucio download {evnt_did}", shell=True, cwd=base_dir)
    else:
        print("have ", evnt_dir ,"already so dont download")
    if not os.path.exists(log_dir):
        print ("will download", log_did, "since dir doesn't exist", log_dir)
        subprocess.call(f"rucio download {log_did}", shell=True, cwd=evnt_dir)
    else:
        print("have ", evnt_dir ,"already so dont download")
    # if job were reran tehre will be several logs and last one comes from succesfull
    untared_dir_cand = glob.glob(f"{evnt_dir}/*.log/tarball_PandaJob*")
    if len(untared_dir_cand)==0:
        tar_file_candidates = sorted(glob.glob(f"{log_dir}/*log.tgz"))
        tar_file = os.path.basename(tar_file_candidates[-1])
        print("have candidates for tar files", tar_file_candidates)
        print("from them select the file to untar", tar_file)
        subprocess.call(f"tar -xvf {tar_file}", shell=True, cwd=log_dir)
    else:
        print("dont untar since did it before res in ", untared_dir_cand[0])

def get_evnt_log_files(i_job_name):
    evnt_did, evnt_dir, log_did, log_dir_before_untar = lu.get_envt_log_names_dirs(base_dir,i_job_name)
    evnt_file = -1
    log_file = -1
    if os.path.exists(evnt_dir) and os.path.exists(log_dir_before_untar):
        print("directories for evnt and log exist")
        evnt_candidates = glob.glob(evnt_dir + "/*EVNT.root")
        log_candidates = glob.glob(log_dir_before_untar + "/tarball_PandaJob*/log.generate")
        print("evnt candidates of len",len(evnt_candidates), evnt_candidates)
        print("log candidates of len",len(log_candidates), log_candidates)
        if len(evnt_candidates)==1 and len(log_candidates)==1:
            evnt_file = evnt_candidates[0]
            log_file =  log_candidates[0]
    else:
        print("directories for evnt and log DOESN:T exist")

    print("returning evnt file", evnt_file)
    print("returning log file", log_file)
    return evnt_file, log_file

def save_fid_xsec_root_hists(DOCUT_str, mydir, xsec_fb):
    yoda_f_str = mydir + "MyOutput.yoda.gz"
    if os.path.exists(yoda_f_str):
        yoda_f = yoda.read(yoda_f_str)
        print("reading from yoda file ", yoda_f_str)

        root_file = mydir + "/hists.root"
        if not os.path.exists(root_file) or not os.path.exists(mydir + "xsec_fb.txt"): proceed = 1
        elif (os.path.exists(root_file) or os.path.exists(mydir + "xsec_fb.txt")) and opts.runAgain=="yes": proceed = 1
        else: proceed = 0
        if proceed:
            # save fid xsec
            integral = yoda_f[f"/VBS_CROSS_TERMS:{DOCUT_str}/m_tagjets"].integral()
            lu.save_xsec_frac_prod(mydir,xsec_fb,integral)
            # save hists in root for further plotting
            root_file = ROOT.TFile(root_file,"UPDATE")
            for i_hist in hists_to_root:
                h_yoda =  yoda_f[f"/VBS_CROSS_TERMS:{DOCUT_str}/" + i_hist]
                h_root = lu.yoda_to_root_1d(h_yoda, i_hist)
                h_root.Write("", ROOT.TObject.kOverwrite)
            root_file.Close()
        else:
            print("dont do xsec and hisst to root")
    else:
        print("dont see yoda file in dir ", mydir)

def save_hists_log_get_xsec_after_cuts(job_name):
    prepare_grid_files(job_name)
    evnt_file, log_file = get_evnt_log_files(job_name)
    if evnt_file!=-1 and log_file!=-1:
        com_run_rivet_no_cut = lu.get_rivet_com(job_name, evtMax = 20000, DOCUT = "NO", redoRivet=opts.runAgain, redoPlots=opts.runAgain)
        com_run_rivet_with_cut = lu.get_rivet_com(job_name, evtMax = 20000, DOCUT = "YES", redoRivet=opts.runAgain, redoPlots=opts.runAgain)
        print("will run rivet+untar in two ways in sequence with/without cut \n", com_run_rivet_no_cut, "\n",com_run_rivet_with_cut)
        if opts.runNoCuts=="yes":subprocess.call(com_run_rivet_no_cut, shell=True)
        if opts.runWithCuts=="yes": subprocess.call(com_run_rivet_with_cut, shell=True)
        # get xsec after cuts
        evnt_dir = os.path.dirname(evnt_file)
        xsec_fb = lu.get_xsec(log_file)
        if opts.runNoCuts=="yes": save_fid_xsec_root_hists("DOCUT=NO",evnt_dir + "/DOCUT_NO/", xsec_fb)
        if opts.runWithCuts=="yes": save_fid_xsec_root_hists("DOCUT=YES",evnt_dir + "/DOCUT_YES/", xsec_fb)

def main():
    parser = OptionParser()
    parser.add_option("--runNoCuts", default = "no")
    parser.add_option("--runWithCuts", default = "yes")
    parser.add_option("--runAgain", default = "no")
    parser.add_option("--jobName", default = "")
    global opts
    opts, _ = parser.parse_args()

    global base_dir
    global hists_to_root    
    if len(opts.jobName)>0:
        print("##################### \n ############# will work on job", opts.jobName)
        prod_dec, base_dir = lu.find_prod_dec_and_dir(opts.jobName) # dir where all files are stored
        hists_to_root = lu.get_hists_arr(prod_dec)
        save_hists_log_get_xsec_after_cuts(opts.jobName)

    return 0

if __name__ == "__main__":
    main()

########
# download files and logs(+untar) - from last succesfol try per config
############
# if opts.doDownloadAndRivet:
#     tasks = c.get_tasks(limit=1000, days=14, username="Oleksii Kurdysh", status="done") # get already last try since only retry if it failed
#     task_names_temp = [i_task['taskname'].replace("/","") for i_task in tasks if "MadGraph" in i_task['taskname'] and opts.tProd in i_task['taskname']]
#     if ignore_fs02:
#         task_names = [i_task for i_task in task_names_temp if "FS02" not in i_task]
#     else:
#         task_names = task_names_temp

#     conf_dirs_dict = {}
#     for i_name in task_names:
#         print("downloading  and untaring for" ,i_name)
#         evnt_did = i_name + "_EXT0"
#         log_did = i_name + ".log"
#         if not os.path.exists(f"{top_files_dir}/{evnt_did}"):
#             subprocess.call(f"rucio download {evnt_did}", shell=True, cwd=top_files_dir)
#             subprocess.call(f"rucio download {log_did}", shell=True, cwd=f'{top_files_dir}/{evnt_did}')
#             tar_file_candidates = sorted(glob.glob(f"{top_files_dir}/{evnt_did}/{log_did}/*log.tgz"))
#             if len(tar_file_candidates)>=1:
#                 tar_file = os.path.basename(tar_file_candidates[-1]) #-1 because sometimes there are several attemps and only last if succesflyy
#                 subprocess.call(f"tar -xvf {tar_file}", shell=True, cwd=f'{top_files_dir}/{evnt_did}/{log_did}')
#             print("done")
#         else:
#             print("EXT0 folder for this exsits already, do nothing")

#         # check that dir contains both evnt and log files
#         conf_dir_cont = os.listdir(f"{top_files_dir}/{evnt_did}")
#         if len(conf_dir_cont)>=2:
#             evnt_candidates = [iobj for iobj in conf_dir_cont if "EVNT.root" in iobj]
#             log_candidates = [iobj for iobj in conf_dir_cont if ".log" in iobj]
#             if len(evnt_candidates)==1 and len(log_candidates)==1:  valid_dir = True
#             else: valid_dir=False
#         else:
#             valid_dir = False
#         if valid_dir: conf_dirs_dict[i_name] =  os.path.join(os.getcwd(),top_files_dir,evnt_did,'')
#     print("left with dirs", conf_dirs_dict)

#     ##########
#     # run rivet and rivet-mkhtml
#     ###########
#     def rivet_run_conf(conf):
#         conf_dir = f"{top_files_dir}/{conf}_EXT0/"
#         if not os.path.exists(conf_dir + "rivet-plots/") or not os.path.exists(conf_dir + "/MyOutput.yoda.gz"):
#             run_com = "athena rivet_job.py -c 'conf=" + f'"{conf}"' + "'"
#             print("#### will run rivet with", run_com)
#             subprocess.call(run_com, shell=True)

#             print("#### will run mkhtml in dir", conf_dir)
#             subprocess.call("rivet-mkhtml MyOutput.yoda.gz", shell=True, cwd = conf_dir)
#         else:
#             print("dont run rivet here since have tgz html plots already")

#     if len(opts.conf)==0: # run all
#         for num_conf,i_name in enumerate(conf_dirs_dict.keys()):
#             print("##############")
#             print("############# will do config num", num_conf, "out of", len(conf_dirs_dict.keys()), "this time it is", i_name)
#             rivet_run_conf(i_name)
#     else:
#         rivet_run_conf(opts.conf)

# ##########
# # summary table of all the xsec*filt
# ###########
# def get_op_from_dir(mydir):
#     temp1 = mydir[len(f"user.okurdysh.MadGraph_{opts.tProd}_"):]
#     temp2 = temp1[:temp1.find("_EXT0")]
#     if "try" in temp2: temp3 = temp2[:temp2.find("_try")]
#     else: temp3  = temp2
#     arr = temp3.split("_")
#     print("vec for this str", arr)
#     ops = arr[0]
#     regime = arr[1]
#     ops_arr = []
#     if "vs" in ops:
#         ops_arr.append(ops[:ops.find("vs")])
#         ops_arr.append(ops[ops.find("vs")+2:])
#     else:
#         ops_arr.append(ops)
#     return sorted(ops_arr), regime

# xsec_dict = {'SM':{}, 'FULL':{}, 'QUAD':{}, 'CROSS':{}, 'INT':{}}
# for op_dir in os.listdir(top_files_dir):
#     print("#### new conf", op_dir)
#     product_file = os.path.join(top_files_dir,op_dir,"xsec_times_frac_pb.txt") # here pb but later for plots will convert to fb
#     if not os.path.exists(product_file):
#         yoda_path = os.path.join(top_files_dir,op_dir,"MyOutput.yoda.gz")
#         print("for yoda use path", yoda_path)
#         integral = yoda.read(yoda_path)["/VBS_CROSS_TERMS/mjj"].integral() # do one line not to have problems in loop 
#         #
#         log_path = glob.glob(f"{top_files_dir}/{op_dir}/*.log/tarball_PandaJob*/log.generate")[0]
#         print("for log use path", log_path)
#         xsec = get_xsec(log_path)
#         #
#         product = integral*xsec
#         print("for ",op_dir,"decoded",get_op_from_dir(op_dir), "integral", integral,"xsec in pb",xsec, "product",product)
#         f = open(product_file, "w") #since opening yoda a bit slow for 50 configs save when have it
#         f.write(str(product))
#         f.close()
#     else:
#         print("have xsec*frac for this file already")
#     f = open(product_file, "r")
#     fid_xsec_fb = float(f.read()) * 1000
#     print("reading back fid_xsec and converting from pb to fb", fid_xsec_fb)
#     ops_arr, regime = get_op_from_dir(op_dir)
#     print("this is for ops", ops_arr, "in regime", regime)
#     if regime=="CROSS":
#         my_op1 = ops_arr[0]
#         my_op2 = ops_arr[1]
#         if my_op1 not in xsec_dict[regime].keys():
#             xsec_dict[regime][my_op1] = {}
#         xsec_dict[regime][my_op1][my_op2] = fid_xsec_fb 
#     else:
#         my_op = ops_arr[0]
#         xsec_dict[regime][my_op] = fid_xsec_fb

# print("SM", xsec_dict["SM"])
# print("FULL", xsec_dict["FULL"])
# print("QUAD", xsec_dict["QUAD"])
# for i_key1 in xsec_dict["CROSS"].keys():
#     print("CROSS",i_key1,xsec_dict["CROSS"][i_key1])

# #######
# #### convert dict into plot
# ############
# def save_plot(plot,path_to_save, draw_option = "text45", log_scale = False):
#     c=ROOT.TCanvas()
#     plot.Draw(draw_option)
#     if log_scale: ROOT.gPad.SetLogy()
#     c.Modified()
#     c.Update()
#     c.Show()
#     c.SaveAs(path_to_save)

# all_ops =  sorted(list(xsec_dict["QUAD"].keys()))
# print("all ops", all_ops)
# nbins = len(all_ops)
# if len(xsec_dict["SM"])>=1 and len(xsec_dict["FULL"])>=1:
#     SM_ref =  xsec_dict["SM"][list(xsec_dict["SM"].keys())[0]] # anyway they are all the same as it should be
#     print("SM xsec in fb", SM_ref)
#     ############# FULL and in comparsion with SM
#     FULL_h = ROOT.TH2F("FULL_h","FULL_h", nbins,0,nbins,1,0,1)
#     FULL_ratio_SM_h = ROOT.TH2F("FULL_ratio_SM_h","FULL_ratio_SM_h", nbins,0,nbins,1,0,1)
#     for num_bin,i_op in enumerate(all_ops, start=1):
#         FULL_h.GetXaxis().SetBinLabel(num_bin,i_op)
#         FULL_ratio_SM_h.GetXaxis().SetBinLabel(num_bin,i_op)
#         if i_op in xsec_dict["FULL"].keys(): 
#             FULL_h.SetBinContent(num_bin, 1, xsec_dict["FULL"][i_op])
#             FULL_ratio_SM_h.SetBinContent(num_bin, 1, xsec_dict["FULL"][i_op] / SM_ref)
#     save_plot(FULL_h,plotdir + "FULL.pdf")
#     save_plot(FULL_ratio_SM_h,plotdir + "FULL_ratio_SM_h.pdf")

#     # draw on one plot jet_pt1 for FULL and SM 
#     if opts.makeJetPt1Plot:


#         FULL_arr = []
#         SM_arr = []
#         stop_num = 115
#         counter_num = 0
#         for op_dir in os.listdir(top_files_dir):
#             ops_arr, regime = get_op_from_dir(op_dir)
#             my_op = ops_arr[0]
#             yoda_path = os.path.join(top_files_dir,op_dir,"MyOutput.yoda.gz")
#             if regime in ["FULL","SM"]:
#                 if counter_num > stop_num: break
#                 print("for yoda use path", yoda_path)
#                 yoda_f = yoda.read(yoda_path)
#                 h_yoda = yoda_f["/VBS_CROSS_TERMS/pt_jet1"]
#                 h_name = my_op + "_" + regime
#                 h_root = yoda_to_root_1d(h_yoda,h_name)
#                 if regime=="FULL":
#                     counter_num +=1
#                     FULL_arr.append(h_root)
#                 elif regime=="SM":
#                     SM_arr.append(h_root)
#                 save_plot(h_root,plotdir +"/per_op/pt_jet1_" + h_name +".pdf", draw_option = "", log_scale = True)
#                 # print("integrals yoda and root", mjj_h_yoda.integral(), mjj_h_root.Integral())        
            
#         stack = ROOT.THStack("pt_jet1", "pt_jet1")
#         stack.Add(SM_arr[0])
#         SM_arr[0].SetLineColor(1)
#         SM_arr[0].SetMarkerColor(1)
#         SM_arr[0].SetMarkerStyle(59)
#         SM_arr[0].SetMarkerSize(1.0)
#         for num_color,ih in enumerate(FULL_arr,start=2):
#             ih.SetLineColor(num_color); ih.SetMarkerColor(num_color)
#             stack.Add(ih)
#         c = ROOT.TCanvas()
#         stack.Draw("nostack")
#         ROOT.gPad.SetLogy()
#         c.BuildLegend()
#         c.Modified()
#         c.Update()
#         c.Show()
#         c.SaveAs(plotdir + "pt_jet1_hist_FULL_SM.pdf")

# ################ QUAD
# QUAD_h = ROOT.TH2F("QUAD_h","QUAD_h", nbins,0,nbins,1,0,1)
# QUAD_ratio_FULL_h = ROOT.TH2F("QUAD_ratio_FULL_h","QUAD_ratio_FULL_h", nbins,0,nbins,1,0,1)
# for num_bin,i_op in enumerate(all_ops, start=1):
#     QUAD_h.GetXaxis().SetBinLabel(num_bin,i_op)
#     QUAD_ratio_FULL_h.GetXaxis().SetBinLabel(num_bin,i_op)
#     if i_op in xsec_dict["QUAD"].keys():
#         quad =  xsec_dict["QUAD"][i_op]
#         QUAD_h.SetBinContent(num_bin, 1, quad)
#         if i_op in xsec_dict["FULL"].keys(): QUAD_ratio_FULL_h.SetBinContent(num_bin, 1, quad / xsec_dict["FULL"][i_op]) 
# save_plot(QUAD_h,plotdir + "QUAD.pdf")
# save_plot(QUAD_ratio_FULL_h,plotdir + "QUAD_ratio_FULL.pdf")
# ############## CROSS + geometric average
# CROSS_h = ROOT.TH2F("CROSS_h","CROSS_h", nbins,0,nbins,nbins,0,nbins)
# CROSS_geom_QUAD_h = ROOT.TH2F("CROSS_geom_QUAD_h","CROSS_geom_QUAD_h", nbins,0,nbins,nbins,0,nbins)
# CROSS_el_area_ratio_h = ROOT.TH2F("CROSS_el_area_ratio_h","CROSS_el_area_ratio_h", nbins,0,nbins,nbins,0,nbins)
# CROSS_el_area_max_h = ROOT.TH2F("CROSS_el_area_max_h","CROSS_el_area_max_h", nbins,0,nbins,nbins,0,nbins)
# for i_op1 in all_ops:
#     if i_op1 in xsec_dict["CROSS"].keys(): 
#         bin_x = all_ops.index(i_op1) + 1
#         CROSS_h.GetXaxis().SetBinLabel(bin_x,i_op1)
#         CROSS_geom_QUAD_h.GetXaxis().SetBinLabel(bin_x,i_op1)
#         CROSS_el_area_ratio_h.GetXaxis().SetBinLabel(bin_x,i_op1)
#         CROSS_el_area_max_h.GetXaxis().SetBinLabel(bin_x,i_op1)
#         for i_op2 in xsec_dict["CROSS"][i_op1]:
#             bin_y = all_ops.index(i_op2) + 1
#             CROSS_h.GetYaxis().SetBinLabel(bin_y,i_op2)
#             CROSS_geom_QUAD_h.GetYaxis().SetBinLabel(bin_y,i_op2)
#             CROSS_el_area_ratio_h.GetYaxis().SetBinLabel(bin_y,i_op2)
#             CROSS_el_area_max_h.GetYaxis().SetBinLabel(bin_y,i_op2)
#             cross = xsec_dict["CROSS"][i_op1][i_op2]
#             CROSS_h.SetBinContent(bin_x, bin_y, cross)
#             # fill geometric average
#             if i_op1 in xsec_dict["QUAD"].keys() and i_op2 in xsec_dict["QUAD"].keys():
#                 quad1 = xsec_dict["QUAD"][i_op1]
#                 quad2 = xsec_dict["QUAD"][i_op2]
#                 geom_average = cross / math.sqrt(quad1*quad2)
#                 # print("for i_op1 i_op2", i_op1, i_op2, "using quads", quad1, quad2, "and cross",xsec_fid ,"with geom ave", geom_average)
#                 CROSS_geom_QUAD_h.SetBinContent(bin_x, bin_y, geom_average)
#                 ##### ellipses
#                 lumi = 140 # lumi of run2 in 1/fb, expect xsec to be in fb
#                 el_q1_c = lumi * quad1 / 3 # dont square xsec since already squared from madgraph
#                 el_q2_c = lumi * quad2 / 3 # divide by 3 since formula for are expect factor to be 1
#                 el_cross_c = lumi * cross / 3
#                 den_with_cross2 = 4*el_q1_c*el_q2_c - el_cross_c**2
#                 if den_with_cross2 > 0:
#                     area_no_cross = 2 * math.pi / math.sqrt(4*el_q1_c*el_q2_c)
#                     area_with_cross = 2 * math.pi / math.sqrt(den_with_cross2)
#                     max_area = max([area_no_cross, area_with_cross])
#                     area_ratio = area_with_cross/area_no_cross
#                     print(f"for i_op1 i_op2 {i_op1} {i_op2} got areas no cross {area_no_cross:.2f} with cross {area_with_cross:.2f} ratio {area_ratio:.2f} used ellipse q1 q2 cross {el_q1_c:.2f} {el_q2_c:.2f} {el_cross_c:.2f}")
#                 else:
#                     print("for i_op1 i_op2", i_op1, i_op2, "cannot do area sqrt in this case",den_with_cross2)
#                     max_area, area_ratio = -99, -99
#                 CROSS_el_area_ratio_h.SetBinContent(bin_x, bin_y, area_ratio)
#                 CROSS_el_area_max_h.SetBinContent(bin_x, bin_y, max_area)
#                 #plot
#                 # plot_b = 5*max_area #0.15
#                 # eq_no_cross = sympy.Eq(el_q1_c*abc.x**2 + el_q2_c*abc.y**2, 1)
#                 # eq_with_cross = sympy.Eq(el_q1_c*abc.x**2 + el_q2_c*abc.y**2 + cross*abc.x*abc.y, 1)
#                 # plot_no_cross = sympy.plot_implicit(eq_no_cross,(abc.x,-1*plot_b,plot_b),(abc.y,-1*plot_b,plot_b),
#                 #                                 show=False,line_color='blue')
#                 # plot_with_cross = sympy.plot_implicit(eq_with_cross,(abc.x,-1*plot_b,plot_b),(abc.y,-1*plot_b,plot_b),
#                 #                                 show=False,line_color='red')
#                 # plot_no_cross.append(plot_with_cross[0])
#                 # plot_no_cross.save(plotdir + f"/ellipses/el_{i_op1}_{i_op2}.png")
#                 # print("done for this pair of ops")
#                 # # spb.backends.matplotlib.MatplotlibBackend.close(plot_no_cross)
#                 # matplotlib.pyplot.close()


# save_plot(CROSS_h, plotdir + "CROSS.pdf")
# save_plot(CROSS_geom_QUAD_h, plotdir + "CROSS_geom_QUAD.pdf")
# save_plot(CROSS_el_area_ratio_h, plotdir + "CROSS_el_area_ratio_h.pdf")
# save_plot(CROSS_el_area_max_h, plotdir + "CROSS_el_area_max_h.pdf")



