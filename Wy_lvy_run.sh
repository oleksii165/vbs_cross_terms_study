run_1_download_0=0
############## Wy
#python run_chain_all_rivets.py --tGenProd "Wmy" --tGenDec "lvy"  --runINT 1 --runQUAD 1 --runCROSS 1 --runLocally 0 --routine "Wy_John" --cut "SR" --runRivet $run_1_download_0
#python run_chain_all_rivets.py --tGenProd "Wpy" --tGenDec "lvy"  --runINT 1 --runQUAD 1 --runCROSS 1 --runLocally 0 --routine "Wy_John" --cut "SR" --runRivet $run_1_download_0

python run_chain.py --genJobName "user.okurdysh.MadGraph_Wpy_lvy_FM5vsFM7_CROSS_try2" --routine "Wy_John" --cut "SR" --genDoDownload 0 --runLocally 1 --runRivet 1 --saveInfoAfterRivet 1 &
python run_chain.py --genJobName "user.okurdysh.MadGraph_Wpy_lvy_FT5vsFT7_CROSS_try2" --routine "Wy_John" --cut "SR" --genDoDownload 0 --runLocally 1 --runRivet 1 --saveInfoAfterRivet 1 &
python run_chain.py --genJobName "user.okurdysh.MadGraph_Wpy_lvy_FT5vsFT6_CROSS_try2" --routine "Wy_John" --cut "SR" --genDoDownload 0 --runLocally 1 --runRivet 1 --saveInfoAfterRivet 1 &
python run_chain.py --genJobName "user.okurdysh.MadGraph_Wpy_lvy_FT6vsFT7_CROSS_try2" --routine "Wy_John" --cut "SR" --genDoDownload 0 --runLocally 1 --runRivet 1 --saveInfoAfterRivet 1 &
