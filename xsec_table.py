import lib_utils as lu
import pandas as pd
import glob

ops = ["T0","T1","T2","T5","T6","T7","T8","T9",
       "S02","S1",
       "M0","M1","M2","M3","M4","M5","M7"]

anas = ["ZZ", "Zy", "WZ", "Wy"]
base_dirs = [["ZZ_llll"], ["Zy_vvy"], ["WpZ_lllv","WmZ_lllv"], ["Wpy_lvy","Wmy_lvy"]]


def get_prod_xsec(prod_dec, op, order):
    my_dir = f"/sps/atlas/kurdysh/vbs_cross_terms_study/eft_files/{prod_dec}/"
    xsec_fb = -1
    log_files = glob.glob(f"{my_dir}/*{op}*{order}*/*{op}*{order}*.log/tarball*/log.generate")
    if len(log_files)==1:
        xsec_fb = lu.get_xsec(log_files[0], prod_xsec=True)
    return xsec_fb

def fill_missing(mydf):
    mydf = mydf.replace("-1.000","-")
    mydf = mydf.replace("-2.000","-")
    return mydf

df_quad = pd.DataFrame(index=ops, columns=anas)
df_int = pd.DataFrame(index=ops, columns=anas)
for ana, base_dir_arr in zip(anas, base_dirs):
    for i_op in ops:
        op_xsecs_arr_quad = []
        op_xsecs_arr_int = []
        for i_prod_dec in base_dir_arr:
            op_xsecs_arr_quad.append(get_prod_xsec(i_prod_dec, i_op, "QUAD"))
            op_xsecs_arr_int.append(get_prod_xsec(i_prod_dec, i_op, "INT"))
        df_quad.at[i_op, ana] = '{:.3f}'.format(round(sum(op_xsecs_arr_quad), 3))
        df_int.at[i_op, ana] = '{:.3f}'.format(round(sum(op_xsecs_arr_int), 3))

df_quad = fill_missing(df_quad)
df_int = fill_missing(df_int)
odir = "/exp/atlas/kurdysh/vbs_cross_terms_study/plots/sum_tables/"
lu.save_df(df_quad, f"{odir}/sum_production_xsec_QUAD.pdf", save_csv=False, save_pdf=False, save_latex=True)
lu.save_df(df_int, f"{odir}/sum_production_xsec_INT.pdf", save_csv=False, save_pdf=False, save_latex=True)
print("hi")