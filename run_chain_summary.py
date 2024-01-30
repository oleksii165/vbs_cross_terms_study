import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)
# ROOT.gROOT.gErrorIgnoreLevel(ROOT.kWarning)
import os
import lib_utils as lu
import matplotlib.pyplot as plt
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
parser.add_option("--SMOnSumPlots", default = "no")
opts, _ = parser.parse_args()
assert [opts.runNoCuts, opts.runWithCuts]!=["yes","yes"], "once at the time"

prod_dec = f"{opts.tProd}_{opts.tDec}"
_, top_files_dir = lu.find_prod_dec_and_dir(f"user.okurdysh.MadGraph_{prod_dec}_FM0_SM")
docut_dir = "DOCUT_YES" if opts.runWithCuts=="yes" else "DOCUT_NO"
plots_to_save = lu.get_hists_to_draw(prod_dec)
plot_dir = lu.get_plotdir(prod_dec, docut_dir)
start_numev = 10000
big_pairs_cutoff = 1.05

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
frac_dict = {'SM':{}, 'FULL':{}, 'QUAD':{}, 'CROSS':{}, 'INT':{}}
sumw_initial_dict = {'SM':{}, 'FULL':{}, 'QUAD':{}, 'CROSS':{}, 'INT':{}}
for op_dir in [i_obj for i_obj in os.listdir(top_files_dir) if os.path.isdir(top_files_dir + "/" + i_obj)]:
    full_op_dir = os.path.join(top_files_dir,op_dir,docut_dir,"")
    print("#### plotting new conf", op_dir, "full path", full_op_dir)
    _, log_file = lu.get_evnt_log_files(top_files_dir,op_dir[:op_dir.find("_EXT0")])
    ops_arr, regime = lu.get_op_from_dir(op_dir, prod_dec)
    print("this is for ops", ops_arr, "in regime", regime)
    product_file = full_op_dir + "xsec_times_frac_fb.txt"
    frac_file = full_op_dir + "frac_after_cuts.txt"
    hists_file = full_op_dir + "hists.root"

    if not os.path.exists(product_file) or not os.path.exists(frac_file) or not os.path.exists(hists_file) or not os.path.exists(log_file):
        print("didnt find the txt xsec fid and/or root file or frac file or log file for", op_dir)
        continue 

    sumw_in = lu.get_sumw_initial(log_file)
    print("reading back sumw_initial", sumw_in)
    # xsec organization
    with open(product_file, 'r') as f:
        fid_xsec_fb = float(f.read())
    print("reading back fid_xsec in fb", fid_xsec_fb)
    with open(frac_file, 'r') as f:
        frac = float(f.read())    
    print("reading back efficiency", frac)
    if regime=="CROSS":
        my_op1, my_op2  = ops_arr[0], ops_arr[1]
        if my_op1 not in xsec_dict[regime].keys(): xsec_dict[regime][my_op1] = {}
        xsec_dict[regime][my_op1][my_op2] = fid_xsec_fb
        #
        if my_op1 not in frac_dict[regime].keys(): frac_dict[regime][my_op1] = {}
        frac_dict[regime][my_op1][my_op2] = frac
        #
        if my_op1 not in sumw_initial_dict[regime].keys(): sumw_initial_dict[regime][my_op1] = {}
        sumw_initial_dict[regime][my_op1][my_op2] = sumw_in
    else:
        my_op = ops_arr[0]
        xsec_dict[regime][my_op] = fid_xsec_fb
        frac_dict[regime][my_op] = frac
        sumw_initial_dict[regime][my_op] = sumw_in
    # hists organization
    op_hists = lu.read_hists(hists_file, plots_to_save)
    for i_hist_name in plots_to_save:
        if i_hist_name not in op_hists.keys(): continue

        if regime=="CROSS":
            my_op1,my_op2 = ops_arr[0], ops_arr[1]
            if my_op1 not in plots_dict[regime][i_hist_name].keys(): 
                plots_dict[regime][i_hist_name][my_op1] = {}
            plots_dict[regime][i_hist_name][my_op1][my_op2] = op_hists[i_hist_name] 
        else:
            my_op = ops_arr[0]
            plots_dict[regime][i_hist_name][my_op] = op_hists[i_hist_name]
    
print("xsec SM", xsec_dict["SM"])
print("xsec FULL", xsec_dict["FULL"])
print("xsec QUAD", xsec_dict["QUAD"])
for i_key1 in xsec_dict["CROSS"].keys():
    print("xsec CROSS",i_key1,xsec_dict["CROSS"][i_key1])
#
print("saved hists SM", [[i_plot, plots_dict["SM"][i_plot].keys()] for i_plot in plots_to_save])
print("saved hists FULL", [[i_plot, plots_dict["FULL"][i_plot].keys()] for i_plot in plots_to_save])
print("saved hists QUAD", [[i_plot, plots_dict["QUAD"][i_plot].keys()] for i_plot in plots_to_save])
print("saved hists CROSS")
for i_plot in plots_to_save:
    print("for plot", i_plot)
    for i_key1 in plots_dict["CROSS"][i_plot].keys():
        print("have hists", plots_dict["CROSS"][i_plot][i_key1].keys())

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
    # make tables and efficiency plot
    ###########
    QUAD_h = ROOT.TH2F("QUAD_h","QUAD_h", nbins,0,nbins,1,0,1)
    for num_bin,i_op in enumerate(all_ops, start=1):
        QUAD_h.GetXaxis().SetBinLabel(num_bin,i_op)
        if i_op in xsec_dict["QUAD"].keys(): QUAD_h.SetBinContent(num_bin, 1, xsec_dict["QUAD"][i_op])
    lu.save_plot(QUAD_h,plot_dir + "QUAD.pdf")

    CROSS_h = ROOT.TH2F("CROSS_h","CROSS_h", nbins,0,nbins,nbins,0,nbins)
    CROSS_geom_QUAD_h = ROOT.TH2F("CROSS_geom_QUAD_h","CROSS_geom_QUAD_h", nbins,0,nbins,nbins,0,nbins)
    CROSS_el_area_ratio_h = ROOT.TH2F("CROSS_el_area_ratio_h","CROSS_el_area_ratio_h", nbins,0,nbins,nbins,0,nbins)
    pairs_str_avalaible = []
    fracs_quad1 = []
    fracs_quad2 = []
    fracs_cross = []
    fracs_errors_quad1 = []
    fracs_errors_quad2 = []
    fracs_errors_cross = []
    #
    pairs_str_big_c = []
    fracs_quad1_big_c = []
    fracs_quad2_big_c = []
    fracs_cross_big_c = []
    fracs_errors_quad1_big_c = []
    fracs_errors_quad2_big_c = []
    fracs_errors_cross_big_c = []
    for i_op1 in all_ops:
        if i_op1 not in xsec_dict["CROSS"].keys(): continue 

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

            if i_op1 not in xsec_dict["QUAD"].keys() or i_op2 not in xsec_dict["QUAD"].keys(): continue

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
            ##### fractions
            pair_str = lu.get_pair_str(i_op1,i_op2)
            pairs_str_avalaible.append(pair_str)
            i_frac_quad1 = frac_dict["QUAD"][i_op1]
            i_frac_quad2 = frac_dict["QUAD"][i_op2]
            i_frac_cross = frac_dict["CROSS"][i_op1][i_op2]
            print("for pair", pair_str, "founds fracs q1 q2 cross", i_frac_quad1, i_frac_quad2, i_frac_cross)
            fracs_quad1.append(i_frac_quad1)
            fracs_quad2.append(i_frac_quad2)
            fracs_cross.append(i_frac_cross)
            i_frac_quad1_unc = math.sqrt(i_frac_quad1*(1-i_frac_quad1)/start_numev) if 1>i_frac_quad1>0 else 0 
            i_frac_quad2_unc = math.sqrt(i_frac_quad2*(1-i_frac_quad2)/start_numev) if 1>i_frac_quad2>0 else 0
            i_frac_cross_unc = math.sqrt(i_frac_cross*(1-i_frac_cross)/start_numev) if 1>i_frac_cross>0 else 0 
            fracs_errors_quad1.append(i_frac_quad1_unc)
            fracs_errors_quad2.append(i_frac_quad2_unc)
            fracs_errors_cross.append(i_frac_cross_unc)
            # fractions but for big pairs
            if area_ratio>big_pairs_cutoff:
                pairs_str_big_c.append(pair_str)
                fracs_quad1_big_c.append(i_frac_quad1)
                fracs_quad2_big_c.append(i_frac_quad2)
                fracs_cross_big_c.append(i_frac_cross)
                fracs_errors_quad1_big_c.append(i_frac_quad1_unc)
                fracs_errors_quad2_big_c.append(i_frac_quad2_unc)
                fracs_errors_cross_big_c.append(i_frac_cross_unc)
            
    # tables
    lu.save_plot(CROSS_h, plot_dir + "CROSS.pdf")
    lu.save_plot(CROSS_geom_QUAD_h, plot_dir + "CROSS_geom_QUAD.pdf")
    lu.save_plot(CROSS_el_area_ratio_h, plot_dir + "CROSS_el_area_ratio_h.pdf")
    # fractions 
    for i_pair, i_frac_quad1, i_frac_quad2, i_frac_cross in zip(pairs_str_avalaible,fracs_quad1,fracs_quad2,fracs_cross):
        print(f"for pair {i_pair} q1 q2 c fractions are {i_frac_quad1}, {i_frac_quad2}, {i_frac_cross}")
    plt.clf()
    ops_num = range(len(pairs_str_avalaible))
    fig = plt.figure(figsize=(16,6))
    plt.errorbar(ops_num, fracs_quad1, yerr=fracs_errors_quad1, fmt='o', mfc="blue", mec="blue", ecolor="blue", label="QUAD1")
    plt.errorbar(ops_num, fracs_quad2, yerr=fracs_errors_quad2, fmt='o', mfc="green", mec="green", ecolor="green", label="QUAD2")
    plt.errorbar(ops_num, fracs_cross, yerr=fracs_errors_cross, fmt='o', mfc="black", mec="black", ecolor="black", label="CROSS")
    plt.xticks(ops_num, pairs_str_avalaible, fontsize = 7)
    plt.xticks(rotation=90)
    plt.ylabel('fraction of weights')
    plt.xlabel('operator pair')
    for xc in ops_num: plt.axvline(x=xc, color='0.3', linestyle='--', linewidth=0.3)
    plt.legend()
    plt.savefig(plot_dir + "/frac_w_after_cuts.pdf", bbox_inches='tight')
    plt.savefig(plot_dir + "/frac_w_after_cuts.svg", bbox_inches='tight')
    # fractions for big pairs
    plt.clf()
    ops_num_big_c = range(len(pairs_str_big_c))
    fig = plt.figure(figsize=(16,6))
    plt.errorbar(ops_num_big_c, fracs_quad1_big_c, yerr=fracs_errors_quad1_big_c, fmt='o', mfc="blue", mec="blue", ecolor="blue", label="QUAD1")
    plt.errorbar(ops_num_big_c, fracs_quad2_big_c, yerr=fracs_errors_quad2_big_c, fmt='o', mfc="green", mec="green", ecolor="green", label="QUAD2")
    plt.errorbar(ops_num_big_c, fracs_cross_big_c, yerr=fracs_errors_cross_big_c, fmt='o', mfc="black", mec="black", ecolor="black", label="CROSS")
    plt.xticks(ops_num_big_c, pairs_str_big_c, fontsize = 7)
    plt.xticks(rotation=90)
    plt.ylabel('fraction of weights')
    plt.xlabel('operator pair')
    for xc in ops_num_big_c: plt.axvline(x=xc, color='0.3', linestyle='--', linewidth=0.3)
    plt.legend()
    plt.savefig(plot_dir + f"/frac_w_after_cuts_ratio_above_{big_pairs_cutoff}.pdf", bbox_inches='tight')
    plt.savefig(plot_dir + f"/frac_w_after_cuts_ratio_above_{big_pairs_cutoff}.svg", bbox_inches='tight')
    

    ##########
    # make plots
    ###########
    def get_multihist_all_pairs(i_plot_name):
        plot_hist_cross_quads_sm={}
        quad_plot_ops = plots_dict["QUAD"][i_plot_name].keys()
        for i_key1 in plots_dict["CROSS"][i_plot_name]:
            i_key1_keys2 = plots_dict["CROSS"][i_plot_name][i_key1].keys()
            if len(i_key1_keys2)>0:
                for i_key2 in i_key1_keys2:
                    print("have hist", i_plot_name, "for CROSS pair", [i_key1,i_key2])
                    if i_key1 in quad_plot_ops and i_key2 in quad_plot_ops:
                        print("also have CROSS plot", i_plot_name, "for both of these ops")
                        i_plot_cross = plots_dict["CROSS"][i_plot_name][i_key1][i_key2]
                        i_plot_cross = lu.dress_hist(i_plot_cross, f"CROSS_{i_key1}_{i_key2}", 1, xsec_dict["CROSS"][i_key1][i_key2])
                        #
                        i_plot_quad1 = plots_dict["QUAD"][i_plot_name][i_key1]
                        i_plot_quad1 = lu.dress_hist(i_plot_quad1, f"QUAD_{i_key1}", 2, xsec_dict["QUAD"][i_key1])
                        #
                        i_plot_quad2 = plots_dict["QUAD"][i_plot_name][i_key2]
                        i_plot_quad2 = lu.dress_hist(i_plot_quad2, f"QUAD_{i_key2}", 3, xsec_dict["QUAD"][i_key2])
                        #
                        i_op_pair = lu.get_pair_str(i_key1,i_key2)
                        plot_hist_cross_quads_sm[i_op_pair] = [i_plot_cross.Clone(), i_plot_quad1.Clone(), i_plot_quad2.Clone()]
                        if "FM0" in plots_dict["SM"][i_plot_name].keys() and opts.SMOnSumPlots=="yes":
                            i_plot_sm = plots_dict["SM"][i_plot_name]["FM0"]
                            i_plot_sm = lu.dress_hist(i_plot_sm, f"SM", 4, xsec_dict["SM"]["FM0"])
                            plot_hist_cross_quads_sm[i_op_pair].append(i_plot_sm)
        print("int the end for plot", i_plot_name, "found hists for pairs", plot_hist_cross_quads_sm.keys())
        return plot_hist_cross_quads_sm

    # for each distrib draw on same plot QUAD1,QUAD2,CROSS - normalied to same area and with each xsec*filt
    stacks_arr_per_pair = {} # key is pair
    stacks_arr_per_pair_normalized = {}
    for i_plot_name in plots_to_save:
        i_plot_hist_dict = get_multihist_all_pairs(i_plot_name)
        print("received for plot", i_plot_name, "arrays of c,q,q,?sm for pairs", i_plot_hist_dict.keys())
        for i_pair in i_plot_hist_dict.keys():
            i_pair_hists = i_plot_hist_dict[i_pair]
            if len(i_pair_hists):
                # rebin and maybe change bounds
                display_params = lu.get_root_hist_param(i_plot_name)
                for i_hist in i_pair_hists:
                    if display_params[2]!=-1: i_hist.RebinX(display_params[2]) 
                # save normalized
                i_pair_plot_stack = lu.make_stack(i_pair_hists, title = i_plot_name,norm = 1)
                if i_pair not in stacks_arr_per_pair_normalized.keys(): stacks_arr_per_pair_normalized[i_pair] = [i_pair_plot_stack]
                else: stacks_arr_per_pair_normalized[i_pair].append(i_pair_plot_stack)
                # same without normalization
                i_pair_plot_stack_xsec = lu.make_stack(i_pair_hists, title = i_plot_name)
                if i_pair not in stacks_arr_per_pair.keys(): stacks_arr_per_pair[i_pair] = [i_pair_plot_stack_xsec]
                else: stacks_arr_per_pair[i_pair].append(i_pair_plot_stack_xsec)
                
    def draw_booklets(stacks_arr, outdir):
        for i_pair in stacks_arr.keys():
            stacks_pair = stacks_arr[i_pair]
            c=ROOT.TCanvas()
            c.Divide(4,2)
            for num_canvas, i_stack in enumerate(stacks_pair,start=1):
                display_params = lu.get_root_hist_param(i_stack.GetName().split("/")[0])
                c.cd(num_canvas)
                i_stack.Draw("nostack")
                if display_params[0]!=-1: i_stack.GetXaxis().SetRangeUser(display_params[0], display_params[1]) 
                ROOT.gPad.BuildLegend()
                ROOT.gPad.SetLogy()
            c.Modified()
            c.Update()
            c.Show()
            c.SaveAs(f"{outdir}/{i_pair}.pdf")
            c.SaveAs(f"{outdir}/svg/{i_pair}.svg")

    draw_booklets(stacks_arr_per_pair_normalized, lu.get_bookletdir(plot_dir, normalized="yes"))
    draw_booklets(stacks_arr_per_pair, lu.get_bookletdir(plot_dir))

    stacks_arr_per_pair_normalized_big_pairs = {}
    big_pairs = lu.get_big_pairs()
    print("big pairs", big_pairs)
    print("all pairs that have", stacks_arr_per_pair_normalized.keys())
    for i_pair in stacks_arr_per_pair_normalized.keys():
        print("print checking if pair", i_pair, "is big")
        if i_pair in big_pairs: 
            print("re-saved")
            stacks_arr_per_pair_normalized_big_pairs[i_pair] = stacks_arr_per_pair_normalized[i_pair] 
    print("big pairs re-saved", stacks_arr_per_pair_normalized_big_pairs.keys())
    draw_booklets(stacks_arr_per_pair_normalized_big_pairs, lu.get_bookletdir(plot_dir, normalized="yes",big_pairs=True))

    