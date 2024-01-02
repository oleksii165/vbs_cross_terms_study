# example of running:   athena rivet_example.py -c 'conf="WmWm_FT0_FULL"'
theApp.EvtMax = -1
print("#####received conf from cmd through -c 'conf=X':  ",conf)

import glob
def find_dir(search_com):
    conf_dir_arr = glob.glob(search_com)
    print("found possibilities for dir", conf_dir_arr)
    conf_dir = conf_dir_arr[0] if len(conf_dir_arr)>=1 else -1  
    if conf_dir == -1: raise ValueError("did not find folder for this config ",search_com)
    return conf_dir
base_dir = "/exp/atlas/kurdysh/vbs_cross_terms_study/"
evnt_conf_dir = find_dir(base_dir + f"/*{conf}*EXT0")
log_conf_dir = find_dir(f"{evnt_conf_dir}/*.log/tarball_PandaJob*")
evnt_file = glob.glob(evnt_conf_dir + "/*EVNT.root")[0]
log_file = log_conf_dir + "/log.generate"
print("will use evnt file and log file ", evnt_file, log_file)

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

with open(log_file) as textf:
    xsec_val, xsec_unit = 1.0 , "pb"
    for line in textf:
        if 'MetaData: cross-section' in line:
            xsec_val = float(line[line.find('=')+1:])
            xsec_unit = line[line.find('(')+1:line.find(')')]
conv_fact_to_pb = {"mb":1e9, "um":1e6, "nb":1e3, "pb":1, "fb":1e-3}
xsec_pb = xsec_val * conv_fact_to_pb[xsec_unit]
print("####found xsec value ",xsec_val,"with unit",xsec_unit,"converting to pb using factor",conv_fact_to_pb[xsec_unit],"get in pb",xsec_pb)
rivet.CrossSection = xsec_pb

#rivet.IgnoreBeamCheck = True
job += rivet

