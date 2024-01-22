import glob
import os

def find_prod_dec_and_dir(conf):
    prod_temp = conf[conf.find("user.okurdysh.MadGraph_")+len("user.okurdysh.MadGraph_"):]
    print("start from string", prod_temp)
    prod_dec = prod_temp[:prod_temp.find("_F")]
    print("from conf found production dec", prod_dec)
    conf_dir = f"/exp/atlas/kurdysh/vbs_cross_terms_study/eft_files/{prod_dec}/"
    print("dir would be", conf_dir)
    return prod_dec, conf_dir

def find_evnt_dir_and_file(search_com):
    conf_dir_arr = glob.glob(search_com)
    print("found possibilities for dir", conf_dir_arr)
    conf_dir = conf_dir_arr[0] if len(conf_dir_arr)>=1 else -1  
    if conf_dir == -1: raise ValueError("did not find folder for this config ",search_com)

    evnt_file_candidates = glob.glob(conf_dir + "/*EVNT.root")
    evnt_file = evnt_file_candidates[0] if len(evnt_file_candidates)>0 else -1

    return conf_dir, evnt_file

def get_conf_cut_dir(evnt_dir, docut):
    mydir = evnt_dir + f"/DOCUT_{docut}/"
    if not os.path.exists(mydir): os.makedirs(mydir) 
    return mydir


