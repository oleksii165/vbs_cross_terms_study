# for i_op in "FT9" "FT8" "FT5" "FT0" "FM2" "FM1" "FM0"
# do
# python run_chain.py --jobName "user.okurdysh.MadGraph_Zy_vvy_${i_op}_QUAD" --runAgain "yes" --runWithCuts "yes" &
# done

#python run_chain.py --jobName "user.okurdysh.MadGraph_WmWm_lvlv_FS1_QUAD_try3" --runAgain "yes" --runWithCuts "yes" --doDownload "no" --runRivet "yes" --evtMax 30000 &

#for job in "user.okurdysh.MadGraph_ssWW_lvlv_FT1_QUAD_try3" "user.okurdysh.MadGraph_ssWW_lvlv_FM7_QUAD_try3" "user.okurdysh.MadGraph_ssWW_lvlv_FM1_QUAD_try3" "user.okurdysh.MadGraph_ssWW_lvlv_FT2_QUAD_try2" "user.okurdysh.MadGraph_ssWW_lvlv_FT0_QUAD_try2" "user.okurdysh.MadGraph_ssWW_lvlv_FS1_QUAD_try2" "user.okurdysh.MadGraph_ssWW_lvlv_FS02_QUAD_try2" "user.okurdysh.MadGraph_ssWW_lvlv_FM0_QUAD_try2" "user.okurdysh.MadGraph_ssWW_lvlv_FM0vsFM7_CROSS_try3" "user.okurdysh.MadGraph_ssWW_lvlv_FS02vsFS1_CROSS_try3" "user.okurdysh.MadGraph_ssWW_lvlv_FT0vsFT1_CROSS_try3" "user.okurdysh.MadGraph_ssWW_lvlv_FT1vsFT2_CROSS_try3" "user.okurdysh.MadGraph_ssWW_lvlv_FM0vsFM1_CROSS_try4" "user.okurdysh.MadGraph_ssWW_lvlv_FM1vsFM7_CROSS_try4" "user.okurdysh.MadGraph_ssWW_lvlv_FT0vsFT2_CROSS_try4"
#do
# python run_chain.py --cut "SR" --genJobName ${job} --evtMax 200000 --runRivet 1 --routine "ATLAS_2023_I2729396" --genDoDownload 0 --runLocally 1 --saveInfoAfterRivet 1 --rivetXsecSet1 1 &
#done

#python run_chain.py --cut "SR" --genJobName "user.okurdysh.MadGraph_ZZ_llll_FM0_SM" --evtMax 2000000 --runRivet 1 --routine "ZZ_llll" --genDoDownload 0 --runLocally 1 --saveInfoAfterRivet 1 &
#wait
#
#
# for i_op in "FM0" "FM1" "FM2" "FM3" "FM4" "FM5" "FM7" "FS02" "FS1" "FT0" "FT1" "FT2" "FT5" "FT6" "FT7" "FT8" "FT9"
# do
#python run_chain.py --cut "SR" --genJobName "user.okurdysh.MadGraph_ZZ_llll_${i_op}_QUAD_try2" --evtMax 2000000 --runRivet 1 --routine "ZZ_llll" --genDoDownload 0 --runLocally 1 --saveInfoAfterRivet 1 &
# done
#wait
#
#  for i_op in "FM0" "FM2" "FM3" "FM4" "FM5" "FM7" "FS02" "FS1" "FT2" "FT5" "FT7" "FT8" "FT9"
#  do
# python run_chain.py --cut "SR" --genJobName "user.okurdysh.MadGraph_ZZ_llll_${i_op}_INT_try2" --evtMax 2000000 --runRivet 1 --routine "ZZ_llll" --genDoDownload 0 --runLocally 1 --saveInfoAfterRivet 1 &
#  done
#wait
#
# python run_chain.py --cut "SR" --genJobName "user.okurdysh.MadGraph_ZZ_llll_FM1_INT_try15" --evtMax 2000000 --runRivet 1 --routine "ZZ_llll" --genDoDownload 0 --runLocally 1 --saveInfoAfterRivet 1 &
# python run_chain.py --cut "SR" --genJobName "user.okurdysh.MadGraph_ZZ_llll_FT0_INT_try3" --evtMax 2000000 --runRivet 1 --routine "ZZ_llll" --genDoDownload 0 --runLocally 1 --saveInfoAfterRivet 1 &
# python run_chain.py --cut "SR" --genJobName "user.okurdysh.MadGraph_ZZ_llll_FT1_INT_try3" --evtMax 2000000 --runRivet 1 --routine "ZZ_llll" --genDoDownload 0 --runLocally 1 --saveInfoAfterRivet 1 &
# python run_chain.py --cut "SR" --genJobName "user.okurdysh.MadGraph_ZZ_llll_FT6_INT_try6" --evtMax 2000000 --runRivet 1 --routine "ZZ_llll" --genDoDownload 0 --runLocally 1 --saveInfoAfterRivet 1 &

# zlly
#  for i_op in "FT0" "FT1" "FT2" "FT5" "FT7" "FT8" "FT9" "FM0" "FM1" "FM2" "FM3" "FM4" "FM5" "FM7"
#  for i_op in "FT6"
#  do
#python run_chain.py --genJobName "user.okurdysh.MadGraph_Zy_lly_${i_op}_INT" --evtMax 100000 --runRivet 1 --routine "ATLAS_2023_I2663725" --runLocally 1 --saveInfoAfterRivet 1 &
#python run_chain.py --genJobName "user.okurdysh.MadGraph_Zy_lly_${i_op}_QUAD_try2" --evtMax 100000 --runRivet 1 --routine "ATLAS_2023_I2663725" --runLocally 1 --saveInfoAfterRivet 1 &
#  done

#  for i_op in "FM0vsFM7" "FM0vsFM3" "FM0vsFM1" "FM0vsFM2" "FM0vsFM5" "FM0vsFM4" "FM1vsFM5"
#  do
#python run_chain.py --genJobName "user.okurdysh.MadGraph_Zy_lly_${i_op}_CROSS_try2" --evtMax 100000 --runRivet 1 --routine "ATLAS_2023_I2663725" --runLocally 1 --saveInfoAfterRivet 1 &
#  done

#  for i_op in "FT5vsFT7" "FT1vsFT5" "FT0vsFT7" "FT5vsFT6" "FT1vsFT8" "FT2vsFT7" "FT7vsFT9" "FM3vsFM7" "FT1vsFT2" "FT7vsFT8" "FT1vsFT9" "FM5vsFM7" "FM1vsFM2" "FT0vsFT8" "FT0vsFT6" "FT2vsFT5" "FM2vsFM4" "FM1vsFM4" "FT5vsFT9" "FT0vsFT2" "FT1vsFT7"
#  do
#python run_chain.py --genJobName "user.okurdysh.MadGraph_Zy_lly_${i_op}_CROSS" --evtMax 100000 --runRivet 1 --routine "ATLAS_2023_I2663725" --runLocally 1 --saveInfoAfterRivet 1 &
#  done

#  for i_op in "FT2vsFT9" "FT2vsFT8" "FM1vsFM3" "FM2vsFM5" "FT5vsFT8" "FT8vsFT9" "FM4vsFM5" "FM4vsFM7" "FT0vsFT5" "FM1vsFM7" "FM2vsFM7" "FT2vsFT6" "FT6vsFT7" "FT6vsFT8" "FT0vsFT1" "FM3vsFM4" "FT1vsFT6" "FT6vsFT9" "FM2vsFM3" "FT0vsFT9"
#  do
#python run_chain.py --genJobName "user.okurdysh.MadGraph_Zy_lly_${i_op}_CROSS" --evtMax 100000 --runRivet 1 --routine "ATLAS_2023_I2663725" --runLocally 1 --saveInfoAfterRivet 1 &
#  done


#python run_chain.py --genJobName "user.okurdysh.MadGraph_Zy_lly_FM3vsFM5_CROSS_try3" --evtMax 100000 --runRivet 1 --routine "ATLAS_2023_I2663725" --runLocally 1 --saveInfoAfterRivet 1 &


# zz llvv theirs
#for i_job in  "user.okurdysh.MadGraph_ZZ_llvv_FT0_INT" "user.okurdysh.MadGraph_ZZ_llvv_FT2_INT" "user.okurdysh.MadGraph_ZZ_llvv_FT5_INT" "user.okurdysh.MadGraph_ZZ_llvv_FT8_INT" "user.okurdysh.MadGraph_ZZ_llvv_FT0_QUAD" "user.okurdysh.MadGraph_ZZ_llvv_FT2_QUAD" "user.okurdysh.MadGraph_ZZ_llvv_FT5_QUAD" "user.okurdysh.MadGraph_ZZ_llvv_FT8_QUAD"
#do
#  python run_chain.py --genJobName "${i_job}" --evtMax 15000 --runRivet 1 --routine "ZZ_llvv" --cut "SR" --genDoDownload 0 --runLocally 1 --extFilesDir "ZZ_llvv_original" --saveInfoAfterRivet 1 &
#done
#
##wait
## zz llvv mine QUAD
#for i_job in "user.okurdysh.MadGraph_ZZ_llvv_FT9_QUAD_try3" "user.okurdysh.MadGraph_ZZ_llvv_FT8_QUAD_try3" "user.okurdysh.MadGraph_ZZ_llvv_FT7_QUAD_try3" "user.okurdysh.MadGraph_ZZ_llvv_FT6_QUAD_try3" "user.okurdysh.MadGraph_ZZ_llvv_FT5_QUAD_try3" "user.okurdysh.MadGraph_ZZ_llvv_FT2_QUAD_try3" "user.okurdysh.MadGraph_ZZ_llvv_FT1_QUAD_try3" "user.okurdysh.MadGraph_ZZ_llvv_FT0_QUAD_try3" "user.okurdysh.MadGraph_ZZ_llvv_FS1_QUAD_try3" "user.okurdysh.MadGraph_ZZ_llvv_FS02_QUAD_try2" "user.okurdysh.MadGraph_ZZ_llvv_FM7_QUAD_try3" "user.okurdysh.MadGraph_ZZ_llvv_FM5_QUAD_try3" "user.okurdysh.MadGraph_ZZ_llvv_FM4_QUAD_try3" "user.okurdysh.MadGraph_ZZ_llvv_FM3_QUAD_try3" "user.okurdysh.MadGraph_ZZ_llvv_FM2_QUAD_try3" "user.okurdysh.MadGraph_ZZ_llvv_FM1_QUAD_try3" "user.okurdysh.MadGraph_ZZ_llvv_FM0_QUAD_try3"
#do
#  python run_chain.py --genJobName "${i_job}" --evtMax 15000 --runRivet 1 --routine "ZZ_llvv" --cut "SR" --genDoDownload 0 --runLocally 1 --extFilesDir "" --saveInfoAfterRivet 1 &
#done
#
#wait
## zz llvv mine INT except FS02
#for i_job in "user.okurdysh.MadGraph_ZZ_llvv_FM7_INT_try3" "user.okurdysh.MadGraph_ZZ_llvv_FM5_INT_try2" "user.okurdysh.MadGraph_ZZ_llvv_FM4_INT_try3" "user.okurdysh.MadGraph_ZZ_llvv_FM3_INT_try2" "user.okurdysh.MadGraph_ZZ_llvv_FM2_INT_try2" "user.okurdysh.MadGraph_ZZ_llvv_FM1_INT_try4" "user.okurdysh.MadGraph_ZZ_llvv_FM0_INT_try2" "user.okurdysh.MadGraph_ZZ_llvv_FT9_INT_try3" "user.okurdysh.MadGraph_ZZ_llvv_FT8_INT_try2" "user.okurdysh.MadGraph_ZZ_llvv_FT7_INT_try2" "user.okurdysh.MadGraph_ZZ_llvv_FT6_INT_try4" "user.okurdysh.MadGraph_ZZ_llvv_FT5_INT_try2" "user.okurdysh.MadGraph_ZZ_llvv_FT2_INT_try5" "user.okurdysh.MadGraph_ZZ_llvv_FT1_INT_try4" "user.okurdysh.MadGraph_ZZ_llvv_FT0_INT_try2" "user.okurdysh.MadGraph_ZZ_llvv_FS1_INT_try2"
#do
#  python run_chain.py --genJobName "${i_job}" --evtMax 15000 --runRivet 1 --routine "ZZ_llvv" --cut "SR" --genDoDownload 0 --runLocally 1 --extFilesDir "" --saveInfoAfterRivet 1 &
#done

#for i_job in  "user.okurdysh.MadGraph_rep_cuts_coupling_pdf_dynscale_qed6_autoptjmjj_cutdecays_ZZ_llvv_FT5_QUAD" "user.okurdysh.MadGraph_rep_cuts_coupling_pdf_dynscale_qed6_ZZ_llvv_FT5_QUAD" "user.okurdysh.MadGraph_rep_cuts_coupling_pdf_dynscale_ZZ_llvv_FT5_QUAD" "user.okurdysh.MadGraph_rep_cuts_coupling_pdf_ZZ_llvv_FT5_QUAD" "user.okurdysh.MadGraph_rep_cuts_coupling_ZZ_llvv_FT5_QUAD" "user.okurdysh.MadGraph_rep_cuts_ZZ_llvv_FT5_QUAD"
for i_job in  "user.okurdysh.MadGraph_rep_cuts_coupling1_pdf_dynscale_qed6_autoptjmjj_cutdecays_sde_hardstrat_ZZ_llvv_FT5_QUAD" "user.okurdysh.MadGraph_rep_cuts_coupling1_pdf_dynscale_qed6_autoptjmjj_cutdecays_sde_hardstrat_r21_ZZ_llvv_FT5_QUAD"
do
#python run_chain.py --genJobName "${i_job}" --evtMax 15000 --runRivet 1 --routine "ZZ_llvv" --cut "SR" --genDoDownload 1 --runLocally 1 --extFilesDir "ZZ_llvv_comparisons" --saveInfoAfterRivet 1 &
python run_chain.py --genJobName "${i_job}" --evtMax 15000 --runRivet 1 --routine "ZZ_llvv" --cut "NO" --genDoDownload 1 --runLocally 1 --extFilesDir "ZZ_llvv_comparisons" --saveInfoAfterRivet 1 &
done