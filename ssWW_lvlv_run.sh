###### ssWW routine submission
run_1_download_0=1
####SR
#ssWW sample
# python run_chain_all_rivets.py --tGenProd "WmWm" --tGenDec "lvlv" --runQUAD 1 --runLocally 0 --routine "ssWW_Original" --cut "SR" --runRivet $run_1_download_0
# python run_chain_all_rivets.py --tGenProd "WpWp" --tGenDec "lvlv" --runQUAD 1 --runLocally 0 --routine "ssWW_Original" --cut "SR" --runRivet $run_1_download_0
##WZ sample
# python run_chain_all_rivets.py --tGenProd "WmZ" --tGenDec "lllv" --runINT 1 --runQUAD 1 --runCROSS 1 --runLocally 0 --routine "ssWW_lvlv" --cut "SR" --runRivet $run_1_download_0
# python run_chain_all_rivets.py --tGenProd "WpZ" --tGenDec "lllv" --runINT 1 --runQUAD 1 --runCROSS 1 --runLocally 0 --routine "ssWW_lvlv" --cut "SR" --runRivet $run_1_download_0
# #####LowmjjCR
# ##ssWW sample
# python run_chain_all_rivets.py --tGenProd "WmWm" --tGenDec "lvlv"  --runINT 1 --runQUAD 1 --runCROSS 1 --runLocally 0 --routine "ssWW_lvlv" --cut "LowmjjCR" --runRivet $run_1_download_0
# python run_chain_all_rivets.py --tGenProd "WpWp" --tGenDec "lvlv"  --runINT 1 --runQUAD 1 --runCROSS 1 --runLocally 0 --routine "ssWW_lvlv" --cut "LowmjjCR" --runRivet $run_1_download_0
# ##WZ sample
# python run_chain_all_rivets.py --tGenProd "WmZ" --tGenDec "lllv" --runINT 1 --runQUAD 1 --runCROSS 1  --runLocally 0 --routine "ssWW_lvlv" --cut "LowmjjCR" --runRivet $run_1_download_0
# python run_chain_all_rivets.py --tGenProd "WpZ" --tGenDec "lllv"  --runINT 1 --runQUAD 1 --runCROSS 1 --runLocally 0 --routine "ssWW_lvlv" --cut "LowmjjCR" --runRivet $run_1_download_0
# #####WZCR - WZ only
# python run_chain_all_rivets.py --tGenProd "WmZ" --tGenDec "lllv" --runINT 1 --runQUAD 1 --runCROSS 1 --runLocally 0 --routine "ssWW_lvlv" --cut "WZCR" --runRivet $run_1_download_0
# python run_chain_all_rivets.py --tGenProd "WpZ" --tGenDec "lllv" --runINT 1 --runQUAD 1 --runCROSS 1 --runLocally 0 --routine "ssWW_lvlv" --cut "WZCR" --runRivet $run_1_download_0

# python run_chain.py --cut "SR" --genJobName "user.okurdysh.MadGraph_WpWp_lvlv_FM1_QUAD_try2" --evtMax 15000 --routine "ATLAS_2023_I2729396" --genDoDownload 0 --runRivet 1 --runLocally 1 --saveInfoAfterRivet 1 --doMakeHtml 0 &
# python run_chain.py --cut "SR" --genJobName "user.okurdysh.MadGraph_WmWm_lvlv_FM1_QUAD_try3" --evtMax 15000 --routine "ATLAS_2023_I2729396" --genDoDownload 0 --runRivet 1 --runLocally 1 --saveInfoAfterRivet 1 --doMakeHtml 0 &
# python run_chain.py --cut "SR" --genJobName "user.okurdysh.MadGraph_WpWp_lvlv_FT2_QUAD_try2" --evtMax 15000 --routine "ATLAS_2023_I2729396" --genDoDownload 0 --runRivet 1 --runLocally 1 --saveInfoAfterRivet 1 --doMakeHtml 0 &
# python run_chain.py --cut "SR" --genJobName "user.okurdysh.MadGraph_WmWm_lvlv_FT2_QUAD_try3" --evtMax 15000 --routine "ATLAS_2023_I2729396" --genDoDownload 0 --runRivet 1 --runLocally 1 --saveInfoAfterRivet 1 --doMakeHtml 0 &

# python run_chain.py --cut "SR" --genJobName "originalJO_user.okurdysh.MadGraph_ssWW_lvlv_FM1_QUAD_try2" --evtMax 15000 --routine "ATLAS_2023_I2729396" --genDoDownload 0 --runRivet 0 --runLocally 1 --saveInfoAfterRivet 1 --doMakeHtml 0 &
# python run_chain.py --cut "SR" --genJobName "originalJO_user.okurdysh.MadGraph_ssWW_lvlv_FT1_QUAD_try2" --evtMax 15000 --routine "ATLAS_2023_I2729396" --genDoDownload 0 --runRivet 0 --runLocally 1 --saveInfoAfterRivet 1 --doMakeHtml 0 &
