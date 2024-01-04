# combine panda downloading, running rivet, making plots for many combinations
from pandaclient import panda_api
c = panda_api.get_api()
import subprocess
import os
import glob
import math
import yoda
import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--conf", default = "")
parser.add_option("--sumPlotsOnly", default = 0)
opts, _ = parser.parse_args()
top_files_dir = "../eft_files/"
plotdir = "../plots/"


########
# download files and logs(+untar) - from last succesfol try per config
############
tasks = c.get_tasks(limit=1000, days=14, username="Oleksii Kurdysh", status="done") # get already last try since only retry if it failed
task_names = [i_task['taskname'].replace("/","") for i_task in tasks if "MadGraph" in i_task['taskname']]

conf_dirs_dict = {}
for i_name in task_names:
    print("downloading  and untaring for" ,i_name)
    evnt_did = i_name + "_EXT0"
    log_did = i_name + ".log"
    if not os.path.exists(f"{top_files_dir}/{evnt_did}"):
        subprocess.call(f"rucio download {evnt_did}", shell=True, cwd=top_files_dir)
        subprocess.call(f"rucio download {log_did}", shell=True, cwd=f'{top_files_dir}/{evnt_did}')
        tar_file = os.path.basename(glob.glob(f"{top_files_dir}/{evnt_did}/{log_did}/*log.tgz")[-1]) #-1 because sometimes there are several attemps and only last if succesflyy
        subprocess.call(f"tar -xvf {tar_file}", shell=True, cwd=f'{top_files_dir}/{evnt_did}/{log_did}')
        print("done")
    else:
        print("EXT0 folder for this exsits already, do nothing")
    conf_dirs_dict[i_name] =  os.path.join(os.getcwd(),top_files_dir,evnt_did,'')
print("left with dirs", conf_dirs_dict)

##########
# run rivet and rivet-mkhtml
###########
def rivet_run_conf(conf):
    run_com = "athena rivet_job.py -c 'conf=" + f'"{conf}"' + "'"
    print("#### will run rivet with", run_com)
    subprocess.call(run_com, shell=True)

    conf_dir = glob.glob(top_files_dir+f"*{conf}*")[0]
    print("#### will run mkhtml in dir", conf_dir)
    subprocess.call("rivet-mkhtml MyOutput.yoda.gz", shell=True, cwd = conf_dir)

if not opts.sumPlotsOnly:
    if len(opts.conf)==0: # run all
        for num_conf,i_name in enumerate(conf_dirs_dict.keys()):
            print("##############")
            print("############# will do config num", num_conf, "out of", len(conf_dirs_dict.keys()), "this time it is", i_name)
            rivet_run_conf(i_name)
    else:
        rivet_run_conf(opts.conf)

##########
# summary table of all the xsec*filt
###########
def get_op_from_dir(mydir):
    temp1 = mydir[len("user.okurdysh.MadGraph_WmWm_"):]
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

def get_xsec(log_file):
    with open(log_file) as textf:
        xsec_val, xsec_unit = 1.0 , "pb"
        for line in textf:
            if 'MetaData: cross-section' in line:
                xsec_val = float(line[line.find('=')+1:])
                xsec_unit = line[line.find('(')+1:line.find(')')]
    conv_fact_to_pb = {"mb":1e9, "um":1e6, "nb":1e3, "pb":1, "fb":1e-3}
    xsec_pb = xsec_val * conv_fact_to_pb[xsec_unit]
    print("found xsec value ",xsec_val,"with unit",xsec_unit,"converting to pb using factor",conv_fact_to_pb[xsec_unit],"get in pb",xsec_pb)
    return xsec_pb


xsec_dict = {'SM':{}, 'FULL':{}, 'QUAD':{}, 'CROSS':{}}
for op_dir in os.listdir(top_files_dir):
    print("#### new conf", op_dir)
    product_file = os.path.join(top_files_dir,op_dir,"xsec_times_frac_pb.txt") 
    if not os.path.exists(product_file):
        yoda_path = os.path.join(top_files_dir,op_dir,"MyOutput.yoda.gz")
        print("for yoda use path", yoda_path)
        integral = yoda.read(yoda_path)["/VBS_CROSS_TERMS/mjj"].integral() # do one line not to have problems in loop 
        #
        log_path = glob.glob(f"{top_files_dir}/{op_dir}/*.log/tarball_PandaJob*/log.generate")[0]
        print("for log use path", log_path)
        xsec = get_xsec(log_path)
        #
        product = integral*xsec
        print("for ",op_dir,"decoded",get_op_from_dir(op_dir), "integral", integral,"xsec in pb",xsec, "product",product)
        f = open(product_file, "w") #since opening yoda a bit slow for 50 configs save when have it
        f.write(str(product))
        f.close()
    else:
        print("have xsec*frac for this file already")
    f = open(product_file, "r")
    fid_xsec = f.read()
    print("reading back fid_xsec",fid_xsec)
    ops_arr, regime = get_op_from_dir(op_dir)
    print("this is for ops", ops_arr, "in regime", regime)
    if regime=="CROSS":
        my_op1 = ops_arr[0]
        my_op2 = ops_arr[1]
        if my_op1 not in xsec_dict[regime].keys():
            xsec_dict[regime][my_op1] = {}
        xsec_dict[regime][my_op1][my_op2] = float(fid_xsec) 
    else:
        my_op = ops_arr[0]
        xsec_dict[regime][my_op] = float(fid_xsec)

print("SM", xsec_dict["SM"])
print("FULL", xsec_dict["FULL"])
print("QUAD", xsec_dict["QUAD"])
for i_key1 in xsec_dict["CROSS"].keys():
    print("CROSS",i_key1,xsec_dict["CROSS"][i_key1])

#######
#### convert dict into plot
############
def save_plot(plot,path_to_save):
    c=ROOT.TCanvas()
    plot.Draw("text")
    c.Modified()
    c.Update()
    c.Show()
    c.SaveAs(path_to_save)

SM_ref =  xsec_dict["SM"][list(xsec_dict["SM"].keys())[0]] # anyway they are all the same as it should be
print("SM xsec in pb", SM_ref)

all_ops =  sorted(list(xsec_dict["QUAD"].keys()))
print("all ops", all_ops)
print("pos of FM0", all_ops.index("FM0"))
print("pos of FS02", all_ops.index("FS02"))
nbins = len(all_ops)

############# FULL and in comparsion with SM
FULL_h = ROOT.TH2F("FULL_h","FULL_h", nbins,0,nbins,1,0,1)
FULL_minus_SM_h = ROOT.TH2F("FULL_minus_SM_h","FULL_minus_SM_h", nbins,0,nbins,1,0,1)
for num_bin,i_op in enumerate(all_ops, start=1):
    FULL_h.GetXaxis().SetBinLabel(num_bin,i_op)
    FULL_minus_SM_h.GetXaxis().SetBinLabel(num_bin,i_op)
    if i_op in xsec_dict["FULL"].keys(): 
        FULL_h.SetBinContent(num_bin, 1, xsec_dict["FULL"][i_op])
        FULL_minus_SM_h.SetBinContent(num_bin, 1, xsec_dict["FULL"][i_op] - SM_ref)
save_plot(FULL_h,plotdir + "FULL.pdf")
save_plot(FULL_minus_SM_h,plotdir + "FULL_minus_SM_h.pdf")
################ QUAD
QUAD_h = ROOT.TH2F("QUAD_h","QUAD_h", nbins,0,nbins,1,0,1)
QUAD_ratio_FULL_h = ROOT.TH2F("QUAD_ratio_FULL_h","QUAD_ratio_FULL_h", nbins,0,nbins,1,0,1)
for num_bin,i_op in enumerate(all_ops, start=1):
    QUAD_h.GetXaxis().SetBinLabel(num_bin,i_op)
    QUAD_ratio_FULL_h.GetXaxis().SetBinLabel(num_bin,i_op)
    if i_op in xsec_dict["QUAD"].keys():
        quad =  xsec_dict["QUAD"][i_op]
        QUAD_h.SetBinContent(num_bin, 1, quad)
        if i_op in xsec_dict["FULL"].keys(): QUAD_ratio_FULL_h.SetBinContent(num_bin, 1, quad / xsec_dict["FULL"][i_op]) 
save_plot(QUAD_h,plotdir + "QUAD.pdf")
save_plot(QUAD_ratio_FULL_h,plotdir + "QUAD_ratio_FULL.pdf")
############## CROSS + geometric average
CROSS_h = ROOT.TH2F("CROSS_h","CROSS_h", nbins,0,nbins,nbins,0,nbins)
CROSS_geom_QUAD_h = ROOT.TH2F("CROSS_geom_QUAD_h","CROSS_geom_QUAD_h", nbins,0,nbins,nbins,0,nbins)
for i_op1 in all_ops:
    if i_op1 in xsec_dict["CROSS"].keys(): 
        bin_x = all_ops.index(i_op1) + 1
        CROSS_h.GetXaxis().SetBinLabel(bin_x,i_op1)
        CROSS_geom_QUAD_h.GetXaxis().SetBinLabel(bin_x,i_op1)
        for i_op2 in xsec_dict["CROSS"][i_op1]:
            bin_y = all_ops.index(i_op2) + 1
            CROSS_h.GetYaxis().SetBinLabel(bin_y,i_op2)
            CROSS_geom_QUAD_h.GetYaxis().SetBinLabel(bin_y,i_op2)
            xsec_fid = xsec_dict["CROSS"][i_op1][i_op2]
            CROSS_h.SetBinContent(bin_x, bin_y, xsec_fid)
            # fill geometric average
            if i_op1 in xsec_dict["QUAD"].keys() and i_op2 in xsec_dict["QUAD"].keys():
                quad1 = xsec_dict["QUAD"][i_op1]
                quad2 = xsec_dict["QUAD"][i_op2]
                geom_average = xsec_fid / math.sqrt(quad1*quad2)
                # print("for i_op1 i_op2", i_op1, i_op2, "using quads", quad1, quad2, "and cross",xsec_fid ,"with geom ave", geom_average)
                CROSS_geom_QUAD_h.SetBinContent(bin_x, bin_y, geom_average)
save_plot(CROSS_h,plotdir + "CROSS.pdf")
save_plot(CROSS_geom_QUAD_h,plotdir + "CROSS_geom_QUAD.pdf")



