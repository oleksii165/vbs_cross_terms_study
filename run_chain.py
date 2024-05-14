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
    if not os.path.exists(evnt_dir):
        print ("will download", evnt_did, "since dir doesn't exist", evnt_dir)
        subprocess.call(f"rucio download {evnt_did}", shell=True, cwd=base_dir)
    else:
        print("have ", evnt_dir ,"already so dont download")
    if not os.path.exists(log_dir):
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

def save_job_infos(mydir, xsec_fb):
    yoda_f_str = mydir + "MyOutput.yoda.gz"

    if not os.path.exists(yoda_f_str): 
        print("dont see yoda file in dir ", mydir, ",return")
        return
    
    yoda_f = yoda.read(yoda_f_str)
    print("reading from yoda file ", yoda_f_str)
    all_hists_in_yoda = [iname  for iname in yoda_f.keys() if "[" not in iname and "RAW" not in iname]
    hists_1h_in_yoda = []
    for i_name in all_hists_in_yoda:
        if yoda_f[i_name].type()=="Histo1D": hists_1h_in_yoda.append(i_name)  
    print("have 1d hists to be saved in root:", hists_1h_in_yoda, "in yoda file", yoda_f_str)

    root_file_name = mydir + "/hists.root"
    # save fid xsec
    lu.write_to_f(mydir + "xsec_fb.txt", xsec_fb) # before cuts

    rivet_dir_name = f"/{opts.routine}:cut={opts.cut}/"
    print("looking for prefix in counter",rivet_dir_name)
    pos_n_in, pos_w_in = yoda_f[f"{rivet_dir_name}pos_w_initial"].numEntries(), yoda_f[f"{rivet_dir_name}pos_w_initial"].sumW()
    neg_n_in, neg_w_in= yoda_f[f"{rivet_dir_name}neg_w_initial"].numEntries(), yoda_f[f"{rivet_dir_name}neg_w_initial"].sumW()
    # save res to txt
    # similarly for clipping
    for i_clip in ["700", "1000", "1500", "2000", "3000", "inf"]:
        i_counter_pos = yoda_f[f"{rivet_dir_name}pos_w_final_clip_{i_clip}"]
        i_pos_n_f, i_pos_w_f = i_counter_pos.numEntries(), i_counter_pos.sumW()
        i_counter_neg = yoda_f[f"{rivet_dir_name}neg_w_final_clip_{i_clip}"]
        i_neg_n_f, i_neg_w_f = i_counter_neg.numEntries(), i_counter_neg.sumW()
        print("for clip", i_clip, "num neg and pos w after cuts", i_pos_w_f, i_neg_w_f, "sum", i_pos_w_f+i_neg_w_f)
        i_frac_cut = (i_pos_w_f+i_neg_w_f) / (pos_w_in+neg_w_in)
        i_frac_pos = i_pos_n_f / pos_n_in if pos_n_in!=0 else 0
        i_frac_neg = i_neg_n_f / neg_n_in if neg_n_in!=0 else 0
        i_frac_cut_er_bar = 1/(pos_n_in+neg_n_in) * math.sqrt(pos_n_in*i_frac_pos*(1-i_frac_pos) + neg_n_in*i_frac_neg*(1-i_frac_neg))    
        print("after cuts have pos and neg events",i_pos_n_f, i_neg_n_f,"which gives eff", i_frac_cut, "with error", i_frac_cut_er_bar)
        #
        clip_out_suff = f"_clip_{i_clip}" 
        lu.write_to_f(mydir + f"xsec_times_frac_fb{clip_out_suff}.txt", xsec_fb*i_frac_cut)
        lu.write_to_f(mydir + f"frac_after_cuts{clip_out_suff}.txt", i_frac_cut)
        lu.write_to_f(mydir + f"frac_after_cuts_error_bar{clip_out_suff}.txt", i_frac_cut_er_bar)

    # save hists in root for further plotting
    root_file = ROOT.TFile(root_file_name,"UPDATE")
    for i_hist in hists_1h_in_yoda: # they are in format '/WpWm_lvlv:DOCUT=YES/leptons_pids'
        h_yoda =  yoda_f[i_hist]
        h_root = lu.yoda_to_root_1d(h_yoda, i_hist.split("/")[-1])
        h_root.Write("", ROOT.TObject.kOverwrite)
    root_file.Close()

def get_ext_in_files(routine):
    if routine=="Zy_vvy":
        files = "RivetZy_vvy.so,Zy_vvy_cuts.json,Zy_vvy_hists.json,jet_hists.json,photon_hists.json"
    return files


def main():
    parser = OptionParser()
    parser.add_option("--cut", default = "SR")
    parser.add_option("--genJobName", default = "")
    parser.add_option("--evtMax", default = 100000000)
    parser.add_option("--runRivet", default = 0, type='int')
    parser.add_option("--routine", default = "ssWW_lvlv")
    parser.add_option("--doDownload", default = 0, type='int')
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

    print("##################### \n ############# will work on job", opts.genJobName)
    prod_dec, base_dir = lu.find_prod_dec_and_dir(opts.genJobName) # dir where all files are stored
    if opts.doDownload and opts.runLocally:
        prepare_grid_files(opts.genJobName)
        evnt_file, log_file = lu.get_evnt_log_files(base_dir,opts.genJobName)
        if evnt_file==-1 or log_file==-1: return 
    
    if opts.runRivet:
        if opts.runLocally:
            run_com = f'athena rivet_job.py '
            run_com += f'''-c 'runLocally=1;conf="{opts.genJobName}";routine="{opts.routine}";cut="{opts.cut}"' '''
            run_com += f'--evtMax {opts.evtMax}'
        else:
            ext_files_str = get_ext_in_files(opts.routine)
            run_com = 'pathena rivet_job.py '
            run_com += f'''-c 'runLocally=0;conf="{opts.genJobName}";routine="{opts.routine}";cut="{opts.cut}"' '''
            run_com += f'--extOutFile=MyOutput.yoda.gz --extFile={ext_files_str} '
            run_com += f'--inDS={opts.genJobName}_EXT0 --outDS={lu.get_rivet_job_name(opts.genJobName,opts.routine,opts.cut)}'
        print("#### will run rivet with", run_com)
        subprocess.call(run_com, shell=True)

    # todo implement saving files from running on GRID
    if opts.runLocally and opts.saveInfoAfterRivet:
        # get xsec after cuts
        evnt_files, log_file = lu.get_evnt_log_files(base_dir,opts.genJobName)
        evnt_dir = os.path.dirname(evnt_files[0])
        xsec_fb = lu.get_xsec(log_file)
        rivet_out_dir = lu.get_conf_cut_dir(evnt_dir, opts.routine, opts.cut)
        save_job_infos(rivet_out_dir, xsec_fb)
        if opts.doMakeHtml:
            plot_com = f"rivet-mkhtml MyOutput.yoda.gz:'Title={prod_dec}' --no-ratio"
            print("#### will run mkhtml in dir", rivet_out_dir, "with com", plot_com)
            subprocess.call(plot_com, shell=True, cwd = rivet_out_dir)

    return 0

if __name__ == "__main__":
    main()