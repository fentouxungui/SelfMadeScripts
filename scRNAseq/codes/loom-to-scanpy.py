import os
os.getcwd() # 查看当前工作目录
os.chdir("D:\\Xi_Lab\Paper-Data\Fly-Cell-Atlas\Downloaded\Gut")

import scanpy as sc

adata = sc.read_loom("./s_fca_biohub_gut_10x.loom",
                     validate=False,
                     sparse=True,
                     X_name='',
                     obs_names='CellID',
                     var_names='Gene')

for key in lc.col_attrs.keys():
    print(lc.col_attrs[key].dtype)



for key in lc.row_attrs.keys():
    if len(list(lc.row_attrs[key].dtype) == 1:
        print("Data: " + key)

atest = [('Atf3_(+)-track', '<i8'), ('BEAF-32_(+)-track', '<i8'), ('Blimp-1_(-)-track', '<i8'), ('CG11504_(+)-track', '<i8'), ('CG15601_(+)-track', '<i8'), ('CG2116_(+)-track', '<i8'), ('CG8281_(+)-track', '<i8'), ('CTCF_(+)-track', '<i8'), ('Clamp_(+)-track', '<i8'), ('CrebA_(+)-track', '<i8'), ('D19B_(+)-track', '<i8'), ('Dp_(+)-track', '<i8'), ('E_(bx)_(+)-track', '<i8'), ('ERR_(+)-track', '<i8'), ('Eip74EF_(-)-track', '<i8'), ('Eip93F_(+)-track', '<i8'), ('Hr51_(+)-track', '<i8'), ('Hr78_(+)-track', '<i8'), ('Hr78_(-)-track', '<i8'), ('Jra_(+)-track', '<i8'), ('Kr_(+)-track', '<i8'), ('Myb_(+)-track', '<i8'), ('Myc_(+)-track', '<i8'), ('Pdp1_(+)-track', '<i8'), ('Sox14_(+)-track', '<i8'), ('Sry-delta_(+)-track', '<i8'), ('Stat92E_(-)-track', '<i8'), ('Vsx2_(+)-track', '<i8'), ('Xbp1_(+)-track', '<i8'), ('Zif_(+)-track', '<i8'), ('abd-A_(+)-track', '<i8'), ('br_(+)-track', '<i8'), ('cg_(+)-track', '<i8'), ('cg_(-)-track', '<i8'), ('cic_(+)-track', '<i8'), ('cwo_(+)-track', '<i8'), ('exex_(+)-track', '<i8'), ('exex_(-)-track', '<i8'), ('fkh_(+)-track', '<i8'), ('gsb-n_(+)-track', '<i8'), ('jing_(-)-track', '<i8'), ('kay_(+)-track', '<i8'), ('kni_(+)-track', '<i8'), ('kni_(-)-track', '<i8'), ('l_(3)neo38_(-)-track', '<i8'), ('mirr_(-)-track', '<i8'), ('nerfin-1_(+)-track', '<i8'), ('odd_(+)-track', '<i8'), ('phol_(+)-track', '<i8'), ('pros_(+)-track', '<i8'), ('scrt_(+)-track', '<i8'), ('sd_(+)-track', '<i8'), ('su_(Hw)_(+)-track', '<i8'), ('tj_(+)-track', '<i8')]


for key in lc.row_attrs.keys():
    print("Data: " + key)
    print(lc.row_attrs[key][1])
    print(len(lc.row_attrs[key][1]))
    #print(lc.row_attrs[key].dtype)

lc.row_attrs[lc.row_attrs.keys()[1]].ndim

