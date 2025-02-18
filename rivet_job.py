# example of running:   athena rivet_example.py -c 'runLocally=1/0;conf="user.okurdysh.MadGraph_WmWm_lvlv_FT0_FULL";routine="WmWm_lvlv";cut="SR;rivetXsecSet1=1"'

theApp.EvtMax = -1
print("#####received conf from cmd through -c 'runLocally=Z;conf=X;routine=R;cut=Y;rivetXsecSet1=S;ext_files_top_dir=E':  ",runLocally,conf,routine,cut,ext_files_top_dir)
import lib_utils
_, base_dir = lib_utils.find_prod_dec_and_dir(conf,ext_files_top_dir)
if runLocally:
    evnt_conf_dir, evnt_files = lib_utils.find_evnt_dir_and_file(base_dir + f"/*{conf}*EXT0")
    conf_cut_dir = lib_utils.get_conf_cut_dir(evnt_conf_dir, routine, cut)

import AthenaPoolCnvSvc.ReadAthenaPool
if runLocally:
    svcMgr.EventSelector.InputCollections = evnt_files
    print("will run on files", evnt_files)

from AthenaCommon.AlgSequence import AlgSequence
job = AlgSequence()

from xAODEventInfoCnv.xAODEventInfoCnvConf import xAODMaker__EventInfoCnvAlg
job += xAODMaker__EventInfoCnvAlg()

from Rivet_i.Rivet_iConf import Rivet_i
rivet = Rivet_i()
import os
rivet.AnalysisPath = os.environ['PWD']

rivet.Analyses += [f'{routine}:cut={cut}']
rivet.RunName = ''
rivet.HistoFile = f'{conf_cut_dir}/MyOutput.yoda.gz' if runLocally else 'MyOutput.yoda.gz'
if rivetXsecSet1:
    rivet.CrossSection = 1.0 #xsec_pb
job += rivet

