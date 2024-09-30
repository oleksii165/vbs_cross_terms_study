# for i_op in "FT9" "FT8" "FT5" "FT0" "FM2" "FM1" "FM0"
# do
# python run_chain.py --jobName "user.okurdysh.MadGraph_Zy_vvy_${i_op}_QUAD" --runAgain "yes" --runWithCuts "yes" &
# done

#python run_chain.py --jobName "user.okurdysh.MadGraph_WmWm_lvlv_FS1_QUAD_try3" --runAgain "yes" --runWithCuts "yes" --doDownload "no" --runRivet "yes" --evtMax 30000 &

#for job in "user.okurdysh.MadGraph_ssWW_lvlv_FT1_QUAD_try3" "user.okurdysh.MadGraph_ssWW_lvlv_FM7_QUAD_try3" "user.okurdysh.MadGraph_ssWW_lvlv_FM1_QUAD_try3" "user.okurdysh.MadGraph_ssWW_lvlv_FT2_QUAD_try2" "user.okurdysh.MadGraph_ssWW_lvlv_FT0_QUAD_try2" "user.okurdysh.MadGraph_ssWW_lvlv_FS1_QUAD_try2" "user.okurdysh.MadGraph_ssWW_lvlv_FS02_QUAD_try2" "user.okurdysh.MadGraph_ssWW_lvlv_FM0_QUAD_try2" "user.okurdysh.MadGraph_ssWW_lvlv_FM0vsFM7_CROSS_try3" "user.okurdysh.MadGraph_ssWW_lvlv_FS02vsFS1_CROSS_try3" "user.okurdysh.MadGraph_ssWW_lvlv_FT0vsFT1_CROSS_try3" "user.okurdysh.MadGraph_ssWW_lvlv_FT1vsFT2_CROSS_try3" "user.okurdysh.MadGraph_ssWW_lvlv_FM0vsFM1_CROSS_try4" "user.okurdysh.MadGraph_ssWW_lvlv_FM1vsFM7_CROSS_try4" "user.okurdysh.MadGraph_ssWW_lvlv_FT0vsFT2_CROSS_try4"
#do
# python run_chain.py --cut "SR" --genJobName ${job} --evtMax 200000 --runRivet 1 --routine "ATLAS_2023_I2729396" --genDoDownload 0 --runLocally 1 --saveInfoAfterRivet 1 --rivetXsecSet1 1 &
#done

python run_chain.py --cut "SR" --genJobName "user.okurdysh.MadGraph_ZZ_llll_FM0_SM" --evtMax 2000000 --runRivet 1 --routine "ZZ_llll" --genDoDownload 0 --runLocally 1 --saveInfoAfterRivet 1 &
wait


 for i_op in "FM0" "FM1" "FM2" "FM3" "FM4" "FM5" "FM7" "FS02" "FS1" "FT0" "FT1" "FT2" "FT5" "FT6" "FT7" "FT8" "FT9"
 do
python run_chain.py --cut "SR" --genJobName "user.okurdysh.MadGraph_ZZ_llll_${i_op}_QUAD_try2" --evtMax 2000000 --runRivet 1 --routine "ZZ_llll" --genDoDownload 0 --runLocally 1 --saveInfoAfterRivet 1 &
 done
wait

  for i_op in "FM0" "FM2" "FM3" "FM4" "FM5" "FM7" "FS02" "FS1" "FT2" "FT5" "FT7" "FT8" "FT9"
  do
 python run_chain.py --cut "SR" --genJobName "user.okurdysh.MadGraph_ZZ_llll_${i_op}_INT_try2" --evtMax 2000000 --runRivet 1 --routine "ZZ_llll" --genDoDownload 0 --runLocally 1 --saveInfoAfterRivet 1 &
  done
wait

 python run_chain.py --cut "SR" --genJobName "user.okurdysh.MadGraph_ZZ_llll_FM1_INT_try15" --evtMax 2000000 --runRivet 1 --routine "ZZ_llll" --genDoDownload 0 --runLocally 1 --saveInfoAfterRivet 1 &
 python run_chain.py --cut "SR" --genJobName "user.okurdysh.MadGraph_ZZ_llll_FT0_INT_try3" --evtMax 2000000 --runRivet 1 --routine "ZZ_llll" --genDoDownload 0 --runLocally 1 --saveInfoAfterRivet 1 &
 python run_chain.py --cut "SR" --genJobName "user.okurdysh.MadGraph_ZZ_llll_FT1_INT_try3" --evtMax 2000000 --runRivet 1 --routine "ZZ_llll" --genDoDownload 0 --runLocally 1 --saveInfoAfterRivet 1 &
 python run_chain.py --cut "SR" --genJobName "user.okurdysh.MadGraph_ZZ_llll_FT6_INT_try6" --evtMax 2000000 --runRivet 1 --routine "ZZ_llll" --genDoDownload 0 --runLocally 1 --saveInfoAfterRivet 1 &
