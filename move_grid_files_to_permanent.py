# when generating MC resulting EVNTS are stored for 2 weeks - instead move them to permanent storage with
# rucio add-rule user.okurdysh.MadGraph_Zy_vvy_FT2_QUAD_try2.log 1 IN2P3-CC_LOCALGROUPDISK
# somehow cannot be called from clion but ok from normal terminal

from pandaclient import panda_api
from rucio.client import Client
import subprocess
import lib_utils as lu
import os
c = panda_api.get_api()
rucio_client = Client()

prod_dec_to_move = ["Zy_vvy"]
ignore_outer_cross_terms = True
debug_print = False
days_to_search_jobs = 30
def get_com_to_move(container_name):
    return f"rucio add-rule {container_name} 1 IN2P3-CC_LOCALGROUPDISK"

tasks = c.get_tasks(limit=100000000, days=days_to_search_jobs, username="Oleksii Kurdysh", status="done") # get already last try since only retry if it failed
task_names = [i_task['taskname'].replace("/","") for i_task in tasks
              if "rivet" not in i_task['taskname']]

for i_prod_dec in prod_dec_to_move:
    prod_dec_tasks = [i_name for i_name in task_names if i_prod_dec in i_name]
    for i_job_name in prod_dec_tasks:
        ops_arr, order, _ = lu.get_op_from_dir(i_job_name+"_EXT0", i_prod_dec)
        if order not in ["INT","QUAD","CROSS"]:
            if debug_print: print("****drop job with non-standard name which not sure how to work with here", i_job_name)
            continue
        if ignore_outer_cross_terms and order=="CROSS" and ops_arr[0][1]!=ops_arr[1][1]:
            if debug_print: print("****drop outer cross", i_job_name)
            continue
        # drop old job if .log or _EXT0 are empty as useless to have only one
        EXT0_name, log_name = i_job_name+"_EXT0", i_job_name+".log"
        grid_cont_EXT0 = rucio_client.list_content("user.okurdysh", EXT0_name)
        grid_cont_log = rucio_client.list_content("user.okurdysh", log_name)
        if len(list(grid_cont_EXT0))==0 or len(list(grid_cont_log))==0:
            if debug_print: print("****job in which log or EXT0 is not there", i_job_name)
            continue
        print("passed sel then working with job", i_job_name, "that have", ops_arr, order)
        com_move_EXT0, com_move_log = get_com_to_move(EXT0_name), get_com_to_move(log_name)
        print("will call \n",com_move_EXT0, "\n", com_move_log)
        subprocess.call(com_move_EXT0, shell=True)
        subprocess.call(com_move_log, shell=True)





