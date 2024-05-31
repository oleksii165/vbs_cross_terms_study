import ROOT
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--prod_dec", type="str", default="zvvy")
parser.add_option("--print_bins", type="str", default="no")
parser.add_option("--print_model", type="str", default="no")
opts, _ = parser.parse_args()

bin_pref = {"zvvy": "pred_measurement_", "wz": "pred_measurement_Run2_mTWZvsBDT_"}
n_bins = {"zvvy": 8, "wz": 20}
model = {"zvvy": "clip100_pdf", "wz": "EFT_aQGCnew_allsys_Run2_MtWZ_mJJ_M1_expected_quad_statOnly_unscaled_pdf"}
f = ROOT.TFile(f"../workspaces/{opts.gen_prod_dec}_unclip.root")

# Retrieve workspace from file
w = f.Get("w")

model = w[model[opts.gen_prod_dec]]
if opts.print_model=="yes": model.Print("t")

if opts.print_bins=="yes":
    for binno in range(n_bins[opts.gen_prod_dec]):
        myvar = w[f"{bin_pref[opts.gen_prod_dec]}{binno}"]
        myvar.Print("v")
