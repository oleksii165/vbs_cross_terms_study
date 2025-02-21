import shutil
import os

request_dir = "/lapp_data/atlas/kurdysh/mcjoboptions_branches/comb_VBS_ZZ_llvv_request/100xxx/"
logs_dir = "/lapp_data/atlas/kurdysh/vbs_eft_files/ZZ_llvv_request/"
out_logs_dir = "/lapp_data/atlas/kurdysh/mcjoboptions_branches/comb_VBS_ZZ_llvv_request/logs/"

for i_dsid in sorted(os.listdir(request_dir)):
    print("+++++++++++++++++++++++++++ new dsid", i_dsid)
    i_dsid_dir = f"{request_dir}/{i_dsid}/"
    jo = [i_file for i_file in os.listdir(i_dsid_dir) if i_file.endswith(".py")][0]
    op_order_str = jo.replace("mc.MGPy8EG_aQGC","").replace("_ZZ_llvv.py","")
    print("found jo", jo, "giving op and order", op_order_str)
    log_dir = ""
    for i_logdir in os.listdir(logs_dir):
        if op_order_str in i_logdir:
            log_dir = i_logdir
            break
    print("found correspodning log dir", log_dir)
    log_path = f"{logs_dir}/{log_dir}/log.generate"
    print("will copy", log_path, "to", i_dsid_dir)
    shutil.copy(log_path, i_dsid_dir)





