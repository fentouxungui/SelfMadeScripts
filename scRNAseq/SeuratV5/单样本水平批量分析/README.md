# 单细胞数据分析 单样本水平

## 一些考量

- 分析流程： CellRanger + SeuratV4，所有样本的CellRanger分析结果均位于同一目录下，或者位于不同目录下。
- 不同的样本，分析结果放到不同的目录中，考虑后期不同样本可能有不同的下游分析方案。
- 不同样本的分析参数【代码中要输出所设置的参数】，可以汇总到一个表格中（或网页上）。
- 自动化创建每个样本的分析目录和拷贝相关的分析代码

## 自动化创建分析目录，并拷贝代码

同时保存CellRanger Output Matirx的路径，方便后续Seurat的读取。

``` bash
Rscript prepare.R
```

```bash
# CellRanger outputs 目录
$ tree CellRanger-Outputs -L 1
CellRanger-Outputs
├── E_LI_D115_B1
├── E_LI_D115_B2
├── P_SI_D1_B0
└── P_SI_DAdult_B0

22 directories, 0 files

# prepare.R 运行后的目录结构
$ tree -L 2
.
├── E_LI_D115_B1
│   └── CellRanger_Matrix_path.rds
├── E_LI_D115_B2
│   └── CellRanger_Matrix_path.rds
├── P_SI_DAdult_B0
│   └── CellRanger_Matrix_path.rds
└── test.R

22 directories, 23 files

```

## 每个样本分别运行Seurat

逐个样本运行``Seurat.Rmd``，每个样本分别设定QC cutoff，然后再Knit。

分群，降维，以及计算各群的marker基因。

## 汇总质控结果

```R
Rscript qc-aggregation.R
```





