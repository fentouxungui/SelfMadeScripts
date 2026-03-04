# Trans H5ad file to Seurat



> Error:
>
> ```
> dat <- Read10X_h5("GSE190889_annoted_fullsize.h5")
> Error in `[[.H5File`(infile, paste0(genome, "/data")) : 
>   An object with name matrix/data does not exist in this group
> ```



## Try 1 worked!

> [Using BPCells with Seurat Objects](https://satijalab.org/seurat/articles/seurat5_bpcells_interaction_vignette)

```
library(BPCells)
library(Seurat)
pbmc.df <- open_matrix_anndata_hdf5('path/to/file.h5ad')
write_matrix_dir(
   mat = pbmc.df,
   dir = 'BPCells')
pbmc.mat <- open_matrix_dir(dir = 'BPCells/')
pbmc.mat <- Azimuth:::ConvertEnsembleToSymbol(mat = pbmc.mat, species = "mouse")
# 新加入的，防止后续整合遇到警告：Warning: Converting BPCells matrix to dgCMatrix for integration as on-disk CCA Integration is not currently supported
pbmc.mat <- as(pbmc.mat, 'dgCMatrix') 
pbmc <- CreateSeuratObject(counts = pbmc.mat, project = project.name, min.cells = 3, min.features = 200)
pbmc
```

原代码（看引用的网页）仅仅适用于单样本分析，如果后续还要做整合，会遇到以下警告信息：

```
obj <- IntegrateLayers(
+   object = obj, method = CCAIntegration,
+   orig.reduction = "pca", new.reduction = 'integrated.cca',
+   verbose = FALSE)
Warning: Converting BPCells matrix to dgCMatrix for integration as on-disk CCA Integration is not currently supported
Warning: Converting to a dense matrix may use excessive memory
This message is displayed once every 8 hours.
Warning: Converting BPCells matrix to dgCMatrix for integration as on-disk CCA Integration is not currently supported
Warning: Converting BPCells matrix to dgCMatrix for integration as on-disk CCA Integration is not currently supported
```

> [CCAIntegration Approach in IntegrateLayers does not use InterableMatrix from BPCells](https://github.com/satijalab/seurat/issues/9408)
>
> 目前的Seurat暂不支持BPCells for IntegrateLayer currently！

经尝试： 加入``pbmc.mat <- as(pbmc.mat, 'dgCMatrix') ``可以解决以上警告！

使用BPCells matrix的好处是，貌似得到的SeuratObject对象很小，几乎不占存储空间！



## Try 2 未尝试

> [在R中读入h5ad文件，并转换为seurat对象](https://blog.csdn.net/weixin_44203980/article/details/146424233)

```
library(reticulate)
# library(anndata)
ad <- reticulate::import("anndata")
adata <- ad$read_h5ad(pbmc10kmono)
 
adata
# AnnData object with n_obs × n_vars = 3782 × 13483
#     obs: 'orig.ident', 'n_genes', 'celltypeL0'
#     var: 'features', 'n_cells'
 
adata$T$X[1:5,1:5]
# 5 x 5 sparse Matrix of class "dgRMatrix"
#            AAACAGCCATCCAGGT-1 AAACCAACACAATGCC-1 AAACCAACAGGAACTG-1
# AL627309.5                  .                  .                  .
# LINC01409                   .                  .                  .
# LINC01128                   .                  .                  .
# LINC00115                   1                  .                  .
# NOC2L                       .                  .                  2
#            AAACCAACATAATCCG-1 AAACCAACATTGTGCA-1
# AL627309.5                  1                  .
# LINC01409                   1                  .
# LINC01128                   .                  .
# LINC00115                   .                  .
# NOC2L                       .                  .
 
# 细胞和基因的元数据
cell_meta <- as.data.frame(adata$obs)
gene_meta <- as.data.frame(adata$var)
# 创建Seurat对象 
seurat_obj <- Seurat::CreateSeuratObject(
  counts = adata$T$X,   ## 这里需要转置，因为seurat中以cell作为列，gene作为行
  meta.data = cell_meta,
  assay = "RNA"
)
# 添加基因注释（如基因类型）
seurat_obj[["RNA"]]@meta.features <- gene_meta
```



## Try 3 未尝试

> [Read_10X_h5 error](https://github.com/satijalab/seurat/issues/6998)

```
library(rhdf5)
library(Seurat)

geo <- H5Fopen("~/Downloads/GSE190889_annoted_fullsize.h5")

h5ls(geo)

mat <- as.sparse(geo$matrix$counts)

barcodes <- geo$matrix$barcodes
features <- geo$matrix$genes

rownames(mat) <- features
colnames(mat) <- barcodes


seu <- CreateSeuratObject(counts = mat)

seu[["cell_type"]] <- geo$metadata$cell_type
seu[["health_status"]] <- geo$metadata$health_status
seu[["index"]] <- geo$metadata$index
seu[["internal_id"]] <- geo$metadata$internal_id
seu[["n_counts"]] <- geo$metadata$n_counts
seu[["n_genes"]] <- geo$metadata$n_genes
```

 适用于以下结构数据：

```
h5ls("GSE190889_annoted_fullsize.h5")
       group          name       otype  dclass           dim
0          /        matrix   H5I_GROUP                      
1    /matrix      barcodes H5I_DATASET  STRING         46199
2    /matrix        counts H5I_DATASET   FLOAT 17111 x 46199
3    /matrix         genes H5I_DATASET  STRING         17111
4          /      metadata   H5I_GROUP                      
5  /metadata     cell_type H5I_DATASET  STRING         46199
6  /metadata health_status H5I_DATASET  STRING         46199
7  /metadata         index H5I_DATASET  STRING         46199
8  /metadata   internal_id H5I_DATASET  STRING         46199
9  /metadata      n_counts H5I_DATASET   FLOAT         46199
10 /metadata       n_genes H5I_DATASET INTEGER         46199
```



若为以下结构，暂不知如何修改代码：

```
h5ls(geo)
                group          name       otype  dclass     dim
0                   /             X   H5I_GROUP                
1                  /X          data H5I_DATASET   FLOAT 5801587
2                  /X       indices H5I_DATASET INTEGER 5801587
3                  /X        indptr H5I_DATASET INTEGER    2178
4                   /        layers   H5I_GROUP                
5                   /           obs   H5I_GROUP                
6                /obs        _index H5I_DATASET  STRING    2177
7                /obs        sample   H5I_GROUP                
8         /obs/sample    categories H5I_DATASET  STRING       1
9         /obs/sample         codes H5I_DATASET INTEGER    2177
10                  /          obsm   H5I_GROUP                
11                  /          obsp   H5I_GROUP                
12                  /           uns   H5I_GROUP                
13                  /           var   H5I_GROUP                
14               /var        _index H5I_DATASET  STRING   33560
15               /var feature_types   H5I_GROUP                
16 /var/feature_types    categories H5I_DATASET  STRING       1
17 /var/feature_types         codes H5I_DATASET INTEGER   33560
18               /var  gene_symbols   H5I_GROUP                
19  /var/gene_symbols    categories H5I_DATASET  STRING   33534
20  /var/gene_symbols         codes H5I_DATASET INTEGER   33560
21               /var gene_versions H5I_DATASET  STRING   33560
22               /var        genome   H5I_GROUP                
23        /var/genome    categories H5I_DATASET  STRING       1
24        /var/genome         codes H5I_DATASET INTEGER   33560
25                  /          varm   H5I_GROUP                
26                  /          varp   H5I_GROUP                
```



# anndata R package

> [Conversion from h5ad to Seurat object](https://github.com/satijalab/seurat/issues/9072)

```
library(Seurat)
library(reticulate)
library(anndata)

data <- read_h5ad("s_fca_biohub_head_10x.h5ad")
cds <- CreateSeuratObject(counts = t(as.matrix(data$X)), meta.data = data$obs)
head(cds@meta.data)


library(anndataR)
library(Seurat)
# 1. Read the H5AD file into an AnnData object
adata <- read_h5ad("s_fca_biohub_head_10x.h5ad")

# 2. Convert AnnData to Seurat object, mapping reductions
# The default mapping handles 'X_pca' to 'PCA', 'X_umap' to 'UMAP', etc.
seurat_obj <- as_Seurat(adata, as = "Seurat")
```

# Converting from AnnData to Seurat via h5Seurat

> [Converting from AnnData to Seurat via h5Seurat](https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html)

```
library(SeuratDisk)
Convert("pbmc3k_final.h5ad", dest = "h5seurat", overwrite = TRUE)
pbmc3k <- LoadH5Seurat("pbmc3k_final.h5seurat")
pbmc3k
```







# SCHARD

>  [SCHARD](https://github.com/cellgeni/schard)

```
# load h5ad as Single Cell Experiment
ba16.sce = schard::h5ad2sce('ba16.h5ad')
# load h5ad as Seurat
snhx = schard::h5ad2seurat('sn.heart.h5ad')
# load all visium samples as single Seurat object
visx = schard::h5ad2seurat_spatial('vis.heart.h5ad')
# or load as list of Seurat objects (per slide
visl = schard::h5ad2seurat_spatial('vis.heart.h5ad',simplify = FALSE)
# or load raw counts
snhr = schard::h5ad2seurat('sn.heart.h5ad',use.raw = TRUE)
# raw counts for visium
visr = schard::h5ad2seurat_spatial('vis.heart.h5ad',use.raw = TRUE)

# check that it works
Seurat::SpatialPlot(visx,features = 'total_counts')
Seurat::SpatialPlot(visx,features = 'total_counts',images = 'HCAHeartST11702009')
Seurat::SpatialPlot(visl$HCAHeartST11702010,features = 'total_counts')
plot(colSums(visx),colSums(visr),pch=16) # raw counts are different from normolized ones
Seurat::DimPlot(snhx,group.by = 'cell_state') # the name of reduction is 'Xumap_' (autotranslated from scanpy to Seurat), somehow DimPlot manages to find it, but probably safier to specify it manually with reduction = 'Xumap_'

```





# 案例：Fly cell atlas

```python
# conda activate scvi-tools
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
adata = sc.read_h5ad('/home/xilab/Data_Backup/Data_From_Paper/2022.03.04_Science_PMID35239393_Fly-Cell-Atlas-A-single-nucleus-transcriptomic-atlas-of-the-adult-fruit-fly/head/s_fca_biohub_head_10x.h5ad')
adata
adata.layers # 以前的anndata没有layers
# 确认有无raw counts
adata.raw.X.shape
adata.X.shape
# 试一下能否将normalized data转为raw counts
raw_counts = np.expm1(adata.raw.X)
raw_counts[:5,:5].toarray()

# df_raw_expression = pd.DataFrame(
#     adata.raw.X.toarray(),
#     index=adata.obs.index,
#     columns=adata.raw.var.index
# )
# df_raw_expression.to_csv('normalized_expression.csv')

adata_full = adata.raw.to_adata()
adata_full.X.shape
adata_full.write('head.h5ad')
```

```R
cds <- schard::h5ad2seurat('head.h5ad')
cds
```

