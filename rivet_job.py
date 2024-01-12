# example of running:   athena rivet_example.py -c 'conf="WmWm_FT0_FULL"'
theApp.EvtMax = -1
print("#####received conf from cmd through -c 'conf=X':  ",conf)
prod_temp = conf[conf.find("user.okurdysh.MadGraph_")+len("user.okurdysh.MadGraph_"):]
prod = prod_temp[:prod_temp.find("_")]
print("from conf found production ", prod)

import glob
def find_dir(search_com):
    conf_dir_arr = glob.glob(search_com)
    print("found possibilities for dir", conf_dir_arr)
    conf_dir = conf_dir_arr[0] if len(conf_dir_arr)>=1 else -1  
    if conf_dir == -1: raise ValueError("did not find folder for this config ",search_com)
    return conf_dir
base_dir = f"/exp/atlas/kurdysh/vbs_cross_terms_study/eft_files/{prod}/"
evnt_conf_dir = find_dir(base_dir + f"/*{conf}*EXT0")
evnt_file = glob.glob(evnt_conf_dir + "/*EVNT.root")[0]

import AthenaPoolCnvSvc.ReadAthenaPool
svcMgr.EventSelector.InputCollections = [ evnt_file ]

from AthenaCommon.AlgSequence import AlgSequence
job = AlgSequence()

from xAODEventInfoCnv.xAODEventInfoCnvConf import xAODMaker__EventInfoCnvAlg
job += xAODMaker__EventInfoCnvAlg()

from Rivet_i.Rivet_iConf import Rivet_i
rivet = Rivet_i()
import os
rivet.AnalysisPath = os.environ['PWD']

rivet.Analyses += ['VBS_CROSS_TERMS']
rivet.RunName = ''
rivet.HistoFile = evnt_conf_dir + f'/MyOutput.yoda.gz'
rivet.CrossSection = 1.0 #xsec_pb
#rivet.IgnoreBeamCheck = True
job += rivet

