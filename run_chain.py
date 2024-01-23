# combine panda downloading, running rivet, making plots for many combinations
from pandaclient import panda_api
c = panda_api.get_api()
import subprocess
import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)
import os
import glob
from array import array
import matplotlib
import yoda
import lib_utils as lu
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