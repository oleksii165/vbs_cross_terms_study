import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)
import os
import lib_utils as lu
import math
import subprocess
from pandaclient import panda_api
c = panda_api.get_api()
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--tProd", default = "WmWm")
parser.add_option("--tDec", default = "lvlv")
parser.add_option("--runSMAndFULL", default = "no")
parser.add_option("--runQUADAndCROSS", default = "no")
parser.add_option("--numJobsParallel", default = 15)
parser.add_option("--runNoCuts", default = "no")
parser.add_option("--runWithCuts", default = "yes")
parser.add_option("--runAgain", default = "no")
parser.add_option("--sumPlotsOnly", default = "yes")
opts, _ = parser.parse_args()
assert [opts.runNoCuts, opts.runWithCuts]!=["yes","yes"], "once at the time"

prod_dec = f"{opts.tProd}_{opts.tDec}"
_, top_files_dir = lu.find_prod_dec_and_dir(f"user.okurdysh.MadGraph_{prod_dec}_FM0_SM")
docut_dir = "DOCUT_YES" if opts.runWithCuts=="yes" else "DOCUT_NO"
# plots_to_save = ["pt_tagjet1","eta_tagjets","m_tagjets", "dy_tagjets", "lepton_pt","lepton_eta", "m_ll", "MET"]
plots_to_save = ["m_tagjets", "m_ll"]
plot_dir = lu.get_plotdir(prod_dec, docut_dir)

tasks = c.get_tasks(limit=1000, days=30, username="Oleksii Kurdysh", status="done") # get already last try since only retry if it failed
task_names = [i_task['taskname'].replace("/","") for i_task in tasks if "MadGraph" in i_task['taskname'] and opts.tProd in i_task['taskname'] and opts.tDec in i_task['taskname']]
all_ops, op_pairs = lu.get_ops(include_fs0_2=True)
splitedSize = int(opts.numJobsParallel)
blocks_of_single_ops = [all_ops[x:x + splitedSize] for x in range(0, len(all_ops), splitedSize)]
blocks_of_op_pairs = [op_pairs[x:x + splitedSize] for x in range(0, len(op_pairs), splitedSize)]
############
# get xsec*frac and save hists to root
##############
def get_com(jobname):
    return f'python run_chain.py --jobName "{jobname}" --runAgain "{opts.runAgain}" --runNoCuts "{opts.runNoCuts}" --runWithCuts "{opts.runWithCuts}"'

def call_bloc_proc(op_blocks, eft_config):
    for i_num_block,i_ops_block in enumerate(op_blocks,start=1): # loops over blocks of 5
        print ("############## \n ####calling processing of ops", i_ops_block, f"block {i_num_block} out of {len(op_blocks)}")
        i_block_bashcoms = []
        for i_op in i_ops_block:
            if eft_config!="CROSS":
                i_job = lu.find_last_match_job(task_names, f"{i_op}_{eft_config}")
                if i_job!=-1: i_block_bashcoms.append(get_com(i_job))
                else: print("-------------- did not find done job for", i_job)
            else:
                i_job_order_1 = lu.find_last_match_job(task_names, f"{i_op[0]}vs{i_op[1]}_{eft_config}")
                i_job_order_2 = lu.find_last_match_job(task_names, f"{i_op[1]}vs{i_op[0]}_{eft_config}")
                if i_job_order_1!=-1: i_block_bashcoms.append(get_com(i_job_order_1))
                elif i_job_order_2!=-1: i_block_bashcoms.append(get_com(i_job_order_2))
                else: print("-------------- did not find done in both orders for", i_job_order_1, i_job_order_2)
        procs = [subprocess.Popen(i_bash, shell=True) for i_bash in i_block_bashcoms]
        # wait until all processes are finished
        if len(procs)>0: 
            for p in procs: p.wait()
        print (f"############## \n ####finished with block ops num {i_num_block} out of {len(op_blocks)}")

if opts.sumPlotsOnly!="yes":
    if opts.runSMAndFULL=="yes": 
        call_bloc_proc([["FM0"]], "SM")
        call_bloc_proc(blocks_of_single_ops, "FULL")
        
    if opts.runQUADAndCROSS=="yes": 
        call_bloc_proc(blocks_of_single_ops, "QUAD")
        call_bloc_proc(blocks_of_op_pairs, "CROSS")

###########
# build summary plots out of files created above
################


plots_dict = {'SM':{}, 'FULL':{}, 'QUAD':{}, 'CROSS':{}, 'INT':{}}
for i_eft in plots_dict.keys():
    for i_plot in plots_to_save:
        plots_dict[i_eft][i_plot] = {}
xsec_dict = {'SM':{}, 'FULL':{}, 'QUAD':{}, 'CROSS':{}, 'INT':{}}

for op_dir in os.listdir(top_files_dir):
    full_op_dir = os.path.join(top_files_dir,op_dir,docut_dir,"")
    print("#### plotting new conf", op_dir, "full path", full_op_dir)
    ops_arr, regime = lu.get_op_from_dir(op_dir, prod_dec)
    print("this is for ops", ops_arr, "in regime", regime)
    product_file = full_op_dir + "xsec_times_frac_fb.txt"
    hists_file = full_op_dir + "hists.root"
    if os.path.exists(product_file):
        f = open(product_file, "r")
        fid_xsec_fb = float(f.read())
        print("reading back fid_xsec in fb", fid_xsec_fb)    
        if regime=="CROSS":
            my_op1 = ops_arr[0]
            my_op2 = ops_arr[1]
            if my_op1 not in xsec_dict[regime].keys(): xsec_dict[regime][my_op1] = {}
            xsec_dict[regime][my_op1][my_op2] = fid_xsec_fb 
        else:
            my_op = ops_arr[0]
            xsec_dict[regime][my_op] = fid_xsec_fb
    else:
        print("didnt find the txt xsec fid file for", op_dir)
print("xsec SM", xsec_dict["SM"])
print("xsec FULL", xsec_dict["FULL"])
print("xsec QUAD", xsec_dict["QUAD"])
for i_key1 in xsec_dict["CROSS"].keys():
    print("xsec CROSS",i_key1,xsec_dict["CROSS"][i_key1])

all_ops =  sorted(list(xsec_dict["QUAD"].keys()))
print("all ops", all_ops)
nbins = len(all_ops)
if opts.runSMAndFULL=="yes" and len(xsec_dict["SM"])>=1 and len(xsec_dict["FULL"])>=1:
    SM_ref =  xsec_dict["SM"][list(xsec_dict["SM"].keys())[0]] 
    print("SM xsec in fb", SM_ref)
    lu.write_to_f(top_files_dir+"/SM_xsec_times_frac_fb.txt", SM_ref)
    ############# FULL and in comparsion with SM
    FULL_h = ROOT.TH2F("FULL_h","FULL_h", nbins,0,nbins,1,0,1)
    FULL_ratio_SM_h = ROOT.TH2F("FULL_ratio_SM_h","FULL_ratio_SM_h", nbins,0,nbins,1,0,1)
    for num_bin,i_op in enumerate(all_ops, start=1):
        FULL_h.GetXaxis().SetBinLabel(num_bin,i_op)
        FULL_ratio_SM_h.GetXaxis().SetBinLabel(num_bin,i_op)
        if i_op in xsec_dict["FULL"].keys(): 
            FULL_h.SetBinContent(num_bin, 1, xsec_dict["FULL"][i_op])
            FULL_ratio_SM_h.SetBinContent(num_bin, 1, xsec_dict["FULL"][i_op] / SM_ref)
    lu.save_plot(FULL_h,plot_dir + "FULL.pdf")
    lu.save_plot(FULL_ratio_SM_h,plot_dir + "FULL_ratio_SM_h.pdf")

if opts.runQUADAndCROSS=="yes":
    ##########
    # make tables
    ###########
    QUAD_h = ROOT.TH2F("QUAD_h","QUAD_h", nbins,0,nbins,1,0,1)
    for num_bin,i_op in enumerate(all_ops, start=1):
        QUAD_h.GetXaxis().SetBinLabel(num_bin,i_op)
        if i_op in xsec_dict["QUAD"].keys(): QUAD_h.SetBinContent(num_bin, 1, xsec_dict["QUAD"][i_op])
    lu.save_plot(QUAD_h,plot_dir + "QUAD.pdf")

    CROSS_h = ROOT.TH2F("CROSS_h","CROSS_h", nbins,0,nbins,nbins,0,nbins)
    CROSS_geom_QUAD_h = ROOT.TH2F("CROSS_geom_QUAD_h","CROSS_geom_QUAD_h", nbins,0,nbins,nbins,0,nbins)
    CROSS_el_area_ratio_h = ROOT.TH2F("CROSS_el_area_ratio_h","CROSS_el_area_ratio_h", nbins,0,nbins,nbins,0,nbins)
    for i_op1 in all_ops:
        if i_op1 in xsec_dict["CROSS"].keys(): 
            bin_x = all_ops.index(i_op1) + 1
            CROSS_h.GetXaxis().SetBinLabel(bin_x,i_op1)
            CROSS_geom_QUAD_h.GetXaxis().SetBinLabel(bin_x,i_op1)
            CROSS_el_area_ratio_h.GetXaxis().SetBinLabel(bin_x,i_op1)
            for i_op2 in xsec_dict["CROSS"][i_op1]:
                bin_y = all_ops.index(i_op2) + 1
                CROSS_h.GetYaxis().SetBinLabel(bin_y,i_op2)
                CROSS_geom_QUAD_h.GetYaxis().SetBinLabel(bin_y,i_op2)
                CROSS_el_area_ratio_h.GetYaxis().SetBinLabel(bin_y,i_op2)
                cross = xsec_dict["CROSS"][i_op1][i_op2]
                CROSS_h.SetBinContent(bin_x, bin_y, cross)
                # fill geometric average
                if i_op1 in xsec_dict["QUAD"].keys() and i_op2 in xsec_dict["QUAD"].keys():
                    quad1 = xsec_dict["QUAD"][i_op1]
                    quad2 = xsec_dict["QUAD"][i_op2]
                    geom_average = cross / math.sqrt(quad1*quad2)
                    # print("for i_op1 i_op2", i_op1, i_op2, "using quads", quad1, quad2, "and cross",xsec_fid ,"with geom ave", geom_average)
                    CROSS_geom_QUAD_h.SetBinContent(bin_x, bin_y, geom_average)
                    ##### ellipses
                    lumi = 140 # lumi of run2 in 1/fb, expect xsec to be in fb
                    el_q1_c = lumi * quad1 / 3 # dont square xsec since already squared from madgraph
                    el_q2_c = lumi * quad2 / 3 # divide by 3 since formula for are expect factor to be 1
                    el_cross_c = lumi * cross / 3
                    den_with_cross2 = 4*el_q1_c*el_q2_c - el_cross_c**2
                    if den_with_cross2 > 0:
                        area_no_cross = 2 * math.pi / math.sqrt(4*el_q1_c*el_q2_c)
                        area_with_cross = 2 * math.pi / math.sqrt(den_with_cross2)
                        area_ratio = area_with_cross/area_no_cross
                        print(f"for i_op1 i_op2 {i_op1} {i_op2} got areas no cross {area_no_cross:.2f} with cross {area_with_cross:.2f} ratio {area_ratio:.2f} used ellipse q1 q2 cross {el_q1_c:.2f} {el_q2_c:.2f} {el_cross_c:.2f}")
                    else:
                        print("for i_op1 i_op2", i_op1, i_op2, "cannot do area sqrt in this case",den_with_cross2)
                        area_ratio = -99
                    CROSS_el_area_ratio_h.SetBinContent(bin_x, bin_y, area_ratio)
    lu.save_plot(CROSS_h, plot_dir + "CROSS.pdf")
    lu.save_plot(CROSS_geom_QUAD_h, plot_dir + "CROSS_geom_QUAD.pdf")
    lu.save_plot(CROSS_el_area_ratio_h, plot_dir + "CROSS_el_area_ratio_h.pdf")

    ##########
    # make plots
    ###########
