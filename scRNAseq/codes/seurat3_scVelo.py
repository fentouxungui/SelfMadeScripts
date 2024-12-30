import scvelo as scv
import loompy
import os
print(os.getcwd())#获得当前工作目录
os.chdir("D:\\Xi_Lab\\石静远\\escort单细胞\\results_for_sjy\\seurat3\\percent.mt.02_nCountRNA10000_dims30")
print(os.getcwd())

loom_file = loompy.connect("escort.seurat3_percent.mt.02_nCountRNA10000_dims30.loom")
adata = scv.read("escort.seurat3_percent.mt.02_nCountRNA10000_dims30.loom", cache=True)
