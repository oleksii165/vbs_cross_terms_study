run_1_download_0=1
order_str="--runQUAD 1" # "--runINT 1 --runQUAD 1 --runCROSS 1"
############## WZ
python run_chain_all_rivets.py --tGenProd "WmZ" --tGenDec "lllv" $order_str --runLocally 0 --routine "WZ_lllv" --cut "SR" --runRivet $run_1_download_0
python run_chain_all_rivets.py --tGenProd "WpZ" --tGenDec "lllv" $order_str --runLocally 0 --routine "WZ_lllv" --cut "SR" --runRivet $run_1_download_0
# get back
python run_chain_all_rivets.py --tGenProd "WmZ" --tGenDec "lllv" $order_str --runLocally 0 --routine "WZ_lllv" --cut "SR" --runRivet $run_1_download_0
python run_chain_all_rivets.py --tGenProd "WpZ" --tGenDec "lllv" $order_str --runLocally 0 --routine "WZ_lllv" --cut "SR" --runRivet $run_1_download_0
