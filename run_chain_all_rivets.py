from pandaclient import panda_api
import subprocess
import lib_utils as lu
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--tGenProd", default = "WmWm")
parser.add_option("--tGenDec", default = "lvlv")
parser.add_option("--runSM", default = 0, type="int")
parser.add_option("--runFULL", default = 0, type="int")
parser.add_option("--runINT", default = 0, type="int")
parser.add_option("--runQUAD", default = 0, type="int")
parser.add_option("--runCROSS", default = 0, type="int")
parser.add_option("--skipOuterCROSS", default = 1, type="int")
parser.add_option("--numLocJobsParallel", default = 15, type="int")
parser.add_option("--runLocally", default = 0, type='int')
parser.add_option("--runRivet", default = 0, type='int')
parser.add_option("--routine", default = "WmWm_lvlv")
parser.add_option("--cut", default = "SR")
parser.add_option("--doMakeHtml", default = 0, type="int")

opts, _ = parser.parse_args()
c = panda_api.get_api()
tasks = c.get_tasks(limit=100000000, days=13000, username="Oleksii Kurdysh", status="done") # get already last try since only retry if it failed
task_pref = "MadGraph"
if opts.tGenDec=="llll": task_pref+="Fixed"
task_names = [i_task['taskname'].replace("/","") for i_task in tasks
              if task_pref in i_task['taskname']
              and f"{opts.tGenProd}_" in i_task['taskname']
              and f"{opts.tGenDec}_" in i_task['taskname']
              and "rivet" not in i_task['taskname']]
all_ops, op_pairs = lu.get_ops(include_fs0_2=False)
splitedSize = opts.numLocJobsParallel if opts.runLocally else 1
blocks_of_single_ops = [all_ops[x:x + splitedSize] for x in range(0, len(all_ops), splitedSize)]
blocks_of_op_pairs = [op_pairs[x:x + splitedSize] for x in range(0, len(op_pairs), splitedSize)]
###########
# get xsec*frac and save hists to root
#############
def get_com(jobname):
    if opts.runLocally:
        com = (f'python run_chain.py '
               f'--runLocally 1 --runRivet {opts.runRivet} '
               f'--genJobName "{jobname}" --routine "{opts.routine}" --cut "{opts.cut}" '
               f'--genDoDownload 1 --saveInfoAfterRivet 1 --doMakeHtml {opts.doMakeHtml}')
    else:
        param_download_proc = 0 if opts.runRivet else 1 
        com = (f'python run_chain.py '
               f'--runLocally 0 --runRivet {opts.runRivet} '
               f'--genDoDownload {param_download_proc} --saveInfoAfterRivet {param_download_proc} --doMakeHtml {opts.doMakeHtml} '
               f'--genJobName "{jobname}" --routine "{opts.routine}" --cut "{opts.cut}" ')
    return com

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
                family1, family2 = i_op[0][1], i_op[1][1] 
                if opts.skipOuterCROSS and family1!=family2: continue
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

if opts.runSM: call_bloc_proc([["FM0"]], "SM")
if opts.runFULL: call_bloc_proc(blocks_of_single_ops, "FULL")
if opts.runINT: call_bloc_proc(blocks_of_single_ops, "INT")
if opts.runQUAD: call_bloc_proc(blocks_of_single_ops, "QUAD")
if opts.runCROSS: call_bloc_proc(blocks_of_op_pairs, "CROSS")
