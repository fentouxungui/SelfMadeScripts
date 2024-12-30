import os
os.getcwd() # 查看当前工作目录
os.chdir("D:\\Xi_Lab\Paper-Data\Fly-Cell-Atlas\Downloaded\Gut")

import loompy
ds = loompy.connect("./s_fca_biohub_gut_10x.loom", validate=False)
ds.attrs.keys()
