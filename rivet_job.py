# example of running:   athena rivet_example.py -c 'conf="user.okurdysh.MadGraph_WmWm_lvlv_FT0_FULL";DOCUT="YES"'

theApp.EvtMax = -1
print("#####received conf from cmd through -c 'conf=X;DOCUT=Y':  ",conf, DOCUT)
import lib_utils
prod_dec, base_dir = lib_utils.find_prod_dec_and_dir(conf)
evnt_conf_dir,evnt_files  = lib_utils.find_evnt_dir_and_file(base_dir + f"/*{conf}*EXT0")
conf_cut_dir = lib_utils.get_conf_cut_dir(evnt_conf_dir, DOCUT)

import AthenaPoolCnvSvc.ReadAthenaPool
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

if prod_dec in ["Zy_lly", "Zy_vvy", "ZZ_llll", "ZZ_llvv"]:
    routine = prod_dec
elif prod_dec in ["WmWm_lvlv", "WpWp_lvlv"]:
    routine = "ssWW_lvlv"
elif prod_dec in ["WmZ_lllv", "WpZ_lllv"]:
    routine = "WZ_lllv"
elif prod_dec in ["Wmy_lvy", "Wpy_lvy"]:
    routine = "Wy_lvy"
else:
    raise Exception("dont know which routine to use for this prod_dec", prod_dec)
rivet.Analyses += [f'{routine}:OUTDIR={conf_cut_dir}']
rivet.RunName = ''
rivet.HistoFile = conf_cut_dir + f'/MyOutput.yoda.gz'
rivet.CrossSection = 1.0 #xsec_pb
#rivet.IgnoreBeamCheck = True
job += rivet

