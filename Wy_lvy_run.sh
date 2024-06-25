run_1_download_0=0
############## Wy
python run_chain_all_rivets.py --tGenProd "Wmy" --tGenDec "lvy"  --runINT 1 --runQUAD 1 --runCROSS 1 --runLocally 0 --routine "Wy_John" --cut "SR" --runRivet $run_1_download_0
python run_chain_all_rivets.py --tGenProd "Wpy" --tGenDec "lvy"  --runINT 1 --runQUAD 1 --runCROSS 1 --runLocally 0 --routine "Wy_John" --cut "SR" --runRivet $run_1_download_0

#  python run_chain.py --cut "SR" --genJobName "original_sample_user.okurdysh.MadGraph_Wy_lvy_FT2_QUAD_try2" --evtMax 30000 --routine "Wy_John" --genDoDownload 0 --runRivet 1 --runLocally 1 --saveInfoAfterRivet 1 --doMakeHtml 0 &
#  python run_chain.py --cut "SR" --genJobName "original_sample_user.okurdysh.MadGraph_Wy_lvy_FM2_QUAD_try2" --evtMax 30000 --routine "Wy_John" --genDoDownload 0 --runRivet 1 --runLocally 1 --saveInfoAfterRivet 1 --doMakeHtml 0 &
