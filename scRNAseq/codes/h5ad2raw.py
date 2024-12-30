import os
os.getcwd() # 查看当前工作目录
os.chdir("D:\\Xi_Lab\Paper-Data\Fly-Cell-Atlas\Analysis\Gut\Genomics")

import numpy as np
import pandas as pd
import scanpy as sc

adata = sc.read_h5ad("../../../Downloaded/Gut/gut.h5ad")
adata

adata = adata.raw.to_adata()
adata

adata.write("./gut.raw.h5ad")