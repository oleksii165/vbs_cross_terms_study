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
import math

def prepare_grid_files(i_job_name):
    print("will download+untar files evnt and log for", i_job_name)
    evnt_did, evnt_dir, log_did, log_dir = lu.get_envt_log_names_dirs(base_dir,i_job_name)
    if opts.runLocally:
        if not os.path.exists(evnt_dir):
            print ("will download", evnt_did, "since dir doesn't exist", evnt_dir)
            subprocess.call(f"rucio download {evnt_did}", shell=True, cwd=base_dir)
        else:
            print("have ", evnt_dir ,"already so dont download")

    if not os.path.exists(log_dir):
        if not opts.runLocally: os.makedirs(evnt_dir)
        print ("will download", log_did, "since dir doesn't exist", log_dir)
        subprocess.call(f"rucio download {log_did}", shell=True, cwd=evnt_dir)
    else:
        print("have ", log_dir ,"already so dont download")
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

def prepare_grid_rivet_files(rivet_job_name):
    rivet_job_files = rivet_job_name.replace("/","") + "_EXT0"
    print("will download yoda files for", rivet_job_name)
    _, evnt_dir, _, _ = lu.get_envt_log_names_dirs(base_dir,opts.genJobName)
    conf_cut_dir = lu.get_conf_cut_dir(evnt_dir, opts.routine,opts.cut)
    if len(glob.glob(f"{conf_cut_dir}/*yoda*"))==0:
        print ("will download", rivet_job_files, "since no yoda in ", conf_cut_dir)
        subprocess.call(f"rucio download {rivet_job_files}", shell=True, cwd=conf_cut_dir)
        # dowloading will create another folder but shoter path names keep it in current
        subprocess.call(f"mv {rivet_job_files}/* .", shell=True, cwd=conf_cut_dir)
        subprocess.call(f"rmdir {rivet_job_files}", shell=True, cwd=conf_cut_dir)
    else:
        print("have yoda in", conf_cut_dir ,"already so dont download")

def save_job_infos(mydir, xsec_no_cuts_fb, runLocally):
    merged_yoda_name = mydir + "MyOutput.yoda.gz"
    if not runLocally:
        yoda_unmerged_names = glob.glob(mydir + "/user*MyOutput.yoda.gz")
        if len(yoda_unmerged_names)==0:
            print("dont see yoda  files for merging in dir ", mydir, ",return")
            return
        if not os.path.exists(merged_yoda_name):
            subprocess.call(f"yodamerge -o {merged_yoda_name} {' '.join(yoda_unmerged_names)}", shell=True, cwd=mydir)

    if not os.path.exists(merged_yoda_name):
        print("dont see yoda  files for merging in dir ", mydir, ",return")
        return

    yoda_f = yoda.read(merged_yoda_name)
    print("reading from yoda file ", merged_yoda_name, "wherehist naes will save those to root")
    all_hists_in_yoda = [iname  for iname in yoda_f.keys() if "[" not in iname and "RAW" not in iname]
    hists_1h_in_yoda = []
    for i_name in all_hists_in_yoda:
        if yoda_f[i_name].type()=="Histo1D": hists_1h_in_yoda.append(i_name)  
    print("have 1d hists to be saved in root:", hists_1h_in_yoda, "in yoda file", merged_yoda_name)
    root_file = ROOT.TFile(mydir + "/hists.root","UPDATE")
    for i_hist in hists_1h_in_yoda:
        h_yoda =  yoda_f[i_hist]
        h_root = lu.yoda_to_root_1d(h_yoda, i_hist.split("/")[-1])
        h_root.Write("", ROOT.TObject.kOverwrite)
    root_file.Close()

    lu.write_to_f(mydir + "xsec_fb.txt", xsec_no_cuts_fb) # before cuts
    rivet_internal_conf_name = f"/{opts.routine}:cut={opts.cut}/"
    print("looking for prefix in counter",rivet_internal_conf_name)
    pos_n_in, pos_w_in = yoda_f[f"{rivet_internal_conf_name}pos_w_initial"].numEntries(), yoda_f[f"{rivet_internal_conf_name}pos_w_initial"].sumW()
    neg_n_in, neg_w_in= yoda_f[f"{rivet_internal_conf_name}neg_w_initial"].numEntries(), yoda_f[f"{rivet_internal_conf_name}neg_w_initial"].sumW()
    for i_clip in ["700", "1000", "1500", "2000", "3000", "inf"]:
        i_counter_pos = yoda_f[f"{rivet_internal_conf_name}pos_w_final_clip_{i_clip}"]
        i_pos_n_f, i_pos_w_f = i_counter_pos.numEntries(), i_counter_pos.sumW()
        i_counter_neg = yoda_f[f"{rivet_internal_conf_name}neg_w_final_clip_{i_clip}"]
        i_neg_n_f, i_neg_w_f = i_counter_neg.numEntries(), i_counter_neg.sumW()
        print("for clip", i_clip, "num neg and pos w after cuts", i_pos_w_f, i_neg_w_f, "sum", i_pos_w_f+i_neg_w_f)
        i_frac_cut = (i_pos_w_f+i_neg_w_f) / (pos_w_in+neg_w_in)
        i_frac_pos = i_pos_n_f / pos_n_in if pos_n_in!=0 else 0
        i_frac_neg = i_neg_n_f / neg_n_in if neg_n_in!=0 else 0
        i_frac_cut_er_bar = 1/(pos_n_in+neg_n_in) * math.sqrt(pos_n_in*i_frac_pos*(1-i_frac_pos) + neg_n_in*i_frac_neg*(1-i_frac_neg))    
        print("after cuts have pos and neg events",i_pos_n_f, i_neg_n_f,"which gives eff", i_frac_cut, "with error", i_frac_cut_er_bar)
        #
        clip_out_suff = f"_clip_{i_clip}" 
        lu.write_to_f(mydir + f"xsec_times_frac_fb{clip_out_suff}.txt", xsec_no_cuts_fb * i_frac_cut)
        lu.write_to_f(mydir + f"frac_after_cuts{clip_out_suff}.txt", i_frac_cut)
        lu.write_to_f(mydir + f"frac_after_cuts_error_bar{clip_out_suff}.txt", i_frac_cut_er_bar)

def get_ext_in_files(routine):
    standard_pack = f"Rivet{routine}.so,{routine}_cuts.json,{routine}_hists.json,"
    standard_pack += "jet_hists.json,"
    if routine=="Zy_vvy":
        files = standard_pack + "photon_hists.json"
    elif routine=="WZ_lllv":
        files = standard_pack + "lepton_hists.json"
    elif routine=="Wy_lvy":
        files = standard_pack + "photon_hists.json,lepton_hists.json"
    elif routine=="Wy_John":
        files = f"Rivet{routine}.so"
    elif routine=="ssWW_lvlv":
        files = standard_pack + "lepton_hists.json"
    else:
        raise Exception("dont know files for this analysis", routine)
    return files

def main():
    parser = OptionParser()
    parser.add_option("--cut", default = "SR")
    parser.add_option("--genJobName", default = "")
    parser.add_option("--evtMax", default = 100000000)
    parser.add_option("--runRivet", default = 0, type='int')
    parser.add_option("--routine", default = "ssWW_lvlv")
    parser.add_option("--genDoDownload", default = 0, type='int')
    parser.add_option("--runLocally", default = 0, type='int')
    parser.add_option("--saveInfoAfterRivet", default = 0, type='int')
    parser.add_option("--doMakeHtml", default = 0, type='int')
    global opts
    opts, _ = parser.parse_args()
    print("got opts", opts)

    global base_dir

    if len(opts.genJobName)==0: return
    if "EXT0" in opts.genJobName or ".log" in opts.genJobName:
        print("provide job name without any EXT or .log")
        return

    print("##################### \n ############# will work on gen job", opts.genJobName)
    prod_dec, base_dir = lu.find_prod_dec_and_dir(opts.genJobName) # dir where all files are stored
    if opts.genDoDownload:
        prepare_grid_files(opts.genJobName)
        evnt_file, log_file = lu.get_evnt_log_files(base_dir,opts.genJobName)
        if (evnt_file==-1 or log_file==-1) and opts.runLocally: return
        if log_file==-1 and not opts.runLocally: return

    tasks = c.get_tasks(limit=10000000, days=13000, username="Oleksii Kurdysh")
    def check_job(job_name):
        matched_tasks=[] # for same taskname but with try2, try3 last try should be first
        for i_task in tasks:
            i_name = i_task['taskname']
            if job_name in i_name: matched_tasks.append(i_task)
        print("found tasks with this name pattern", job_name, "are")
        for i_task in matched_tasks:
            print(i_task['status'], " ", i_task['taskname'].replace("/",""))

        if len(matched_tasks)>0:
            last_task = matched_tasks[0]
            return last_task['status'], last_task['taskname'].replace("/",""), last_task['jeditaskid']
        else:
            return -1,-1,-1

    rivet_job_name = lu.get_rivet_job_name(opts.genJobName,opts.routine,opts.cut)
    if opts.runRivet:
        if opts.runLocally:
            run_com = f'athena rivet_job.py '
            run_com += f'''-c 'runLocally=1;conf="{opts.genJobName}";routine="{opts.routine}";cut="{opts.cut}"' '''
            run_com += f'--evtMax {opts.evtMax}'
        else:
            last_job_status, last_job_name, _ = check_job(rivet_job_name)
            new_job_name = lu.get_rivet_resub_name(last_job_name) if last_job_status != -1 else rivet_job_name
            ext_files_str = get_ext_in_files(opts.routine)
            run_com = 'pathena rivet_job.py '
            run_com += f'''-c 'runLocally=0;conf="{opts.genJobName}";routine="{opts.routine}";cut="{opts.cut}"' '''
            run_com += f'--extOutFile=MyOutput.yoda.gz --extFile={ext_files_str} '
            run_com += f'--inDS={opts.genJobName}_EXT0 --outDS={new_job_name}'
        print("#### will run rivet with", run_com)
        subprocess.call(run_com, shell=True)

    if not opts.runLocally and opts.saveInfoAfterRivet:
        job_status, job_name, _ = check_job(rivet_job_name)
        if job_status!='done':
            print("********** rivet task for this job is not finished", rivet_job_name)
            return
        prepare_grid_rivet_files(job_name)

    if opts.saveInfoAfterRivet:
        _, log_file = lu.get_evnt_log_files(base_dir,opts.genJobName)
        _, evnt_dir, _, _ = lu.get_envt_log_names_dirs(base_dir,opts.genJobName)
        xsec_no_cuts_fb = lu.get_xsec(log_file)
        rivet_out_dir = lu.get_conf_cut_dir(evnt_dir, opts.routine, opts.cut)
        save_job_infos(rivet_out_dir, xsec_no_cuts_fb, opts.runLocally)
        if opts.doMakeHtml:
            plot_com = f"rivet-mkhtml MyOutput.yoda.gz:'Title={prod_dec}' --no-ratio"
            print("#### will run mkhtml in dir", rivet_out_dir, "with com", plot_com)
            subprocess.call(plot_com, shell=True, cwd = rivet_out_dir)

    return 0

if __name__ == "__main__":
    main()
