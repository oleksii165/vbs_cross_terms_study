# run rivet together with rivet plotting afterwards, if files were already downloded
import subprocess
import lib_utils
import os
import run_chain
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--conf", default = "user.okurdysh.MadGraph_WmWm_lvlv_FT0_FULL")
parser.add_option("--DOCUT", default = "NO")
parser.add_option("--evtMax", default = 100)
parser.add_option("--redoRivet", default = "no")
parser.add_option("--redoPlots", default = "no")
parser.add_option("--plotConfRivet", default = "/exp/atlas/kurdysh/vbs_cross_terms_study/plotting/VBS_CROSS_TERMS.plot")
opts, _ = parser.parse_args()

prod_dec, base_dir = lib_utils.find_prod_dec_and_dir(opts.conf)
conf_dir, _ = lib_utils.find_evnt_dir_and_file(base_dir + f"/*{opts.conf}*EXT0")
conf_cut_dir = lib_utils.get_conf_cut_dir(conf_dir, opts.DOCUT)
rivet_out_name = conf_cut_dir + f'/MyOutput.yoda.gz'
if not os.path.exists(rivet_out_name): do_rivet = 1
elif os.path.exists(rivet_out_name) and opts.redoRivet=="yes": do_rivet = 1
else: do_rivet = 0
if do_rivet:
    run_com = "athena rivet_job.py -c 'conf=" + f'"{opts.conf}";DOCUT=' + f'"{opts.DOCUT}"' + f"' --evtMax {opts.evtMax}"
    print("#### will run rivet with", run_com)
    subprocess.call(run_com, shell=True)
else:
    print("dont do evnt conversion since redoRivet=", opts.redoRivet, f"and file {rivet_out_name} exists ", os.path.exists(rivet_out_name))

plots_dir = conf_cut_dir + "/rivet-plots/"
if not os.path.exists(plots_dir): do_plots = 1
elif os.path.exists(plots_dir) and opts.redoPlots=="yes": do_plots = 1
else: do_plots = 0
if do_plots:
    plot_com = f"rivet-mkhtml MyOutput.yoda.gz:'Title={prod_dec}' -c {opts.plotConfRivet} --no-ratio"
    print("#### will run mkhtml in dir", conf_cut_dir, "with com", plot_com)
    subprocess.call(plot_com, shell=True, cwd = conf_cut_dir)
else:
    print("dont do plot from yoda since redoPlots=", opts.redoPlots, f"and dir {plots_dir} exists ", os.path.exists(plots_dir))