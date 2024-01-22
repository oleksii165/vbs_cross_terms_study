import lib_utils as lu
import subprocess

from pandaclient import panda_api
c = panda_api.get_api()
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--tProd", default = "WmWm")
parser.add_option("--tDec", default = "lvlv")
parser.add_option("--runSMAndFULL", default = "no")
parser.add_option("--runQUADAndCROSS", default = "no")
parser.add_option("--numJobsParallel", default = 10)
parser.add_option("--runNoCuts", default = "no")
parser.add_option("--runWithCuts", default = "yes")
parser.add_option("--runAgain", default = "no")
opts, _ = parser.parse_args()

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
        for p in procs: p.wait()
        print (f"############## \n ####finished with block ops num {i_num_block} out of {len(op_blocks)}")

if opts.runSMAndFULL=="yes": 
    call_bloc_proc([["FM0"]], "SM")
    call_bloc_proc(blocks_of_single_ops, "FULL")
    
if opts.runQUADAndCROSS=="yes": 
    call_bloc_proc(blocks_of_single_ops, "QUAD")
    call_bloc_proc(blocks_of_op_pairs, "CROSS")



###########
# build summary plots out of files created above
################