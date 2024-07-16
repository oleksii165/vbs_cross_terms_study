run_1_download_0=0
############## WZ
#python run_chain_all_rivets.py --tGenProd "WmZ" --tGenDec "lllv" --runINT 1 --runCROSS 1 --runLocally 0 --routine "WZ_lllv" --cut "SR" --runRivet $run_1_download_0
#python run_chain_all_rivets.py --tGenProd "WpZ" --tGenDec "lllv" --runINT 1 --runCROSS 1 --runLocally 0 --routine "WZ_lllv" --cut "SR" --runRivet $run_1_download_0

#python run_chain.py --genJobName "user.okurdysh.MadGraph_WZ_lllv_FT0_QUAD_try2" --routine "WZ_lllv" --cut "SR" --genDoDownload 0 --runLocally 1 --runRivet 1 --evtMax 40000 --saveInfoAfterRivet 1 &
#python run_chain.py --genJobName "user.okurdysh.MadGraph_WZ_lllv_FT2_QUAD_try2" --routine "WZ_lllv" --cut "SR" --genDoDownload 0 --runLocally 1 --runRivet 1 --evtMax 40000 --saveInfoAfterRivet 1 &
#python run_chain.py --genJobName "user.okurdysh.MadGraph_WZ_lllv_FM0_QUAD_try2" --routine "WZ_lllv" --cut "SR" --genDoDownload 1 --runLocally 1 --runRivet 1 --evtMax 40000 --saveInfoAfterRivet 1 &

#for op in "FT7" "FT6" "FS1" "FM7" "FM5" "FM4" "FM3" "FM2" "FM1"
#do
#  python run_chain.py --genJobName "user.okurdysh.MadGraph_WZ_lllv_${op}_QUAD_try2" --routine "WZ_lllv" --cut "SR" --genDoDownload 1 --runLocally 1 --runRivet 0 &
#done

python search_replacements.py --tGenProd "WpZWmZ" --tGenDec "lllv" --routine "WZ_lllv" --cut "SR" --doINT 1 --doCROSS 0 --doReplacement 0 --doReshuffling 1 --savePdf 0
