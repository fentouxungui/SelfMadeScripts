# 预处理从GEO下载好的单细胞数据

整理成每个样本放到一个目录下，并修改文件名。

起始数据结果：

```
.
├── GSM5623768_barcodes.tsv.gz
├── GSM5623768_features.tsv.gz
├── GSM5623768_matrix.mtx.gz
├── GSM5623769_barcodes.tsv.gz
├── GSM5623769_features.tsv.gz
├── GSM5623769_matrix.mtx.gz
├── GSM5623770_barcodes.tsv.gz
├── GSM5623770_features.tsv.gz
├── GSM5623770_matrix.mtx.gz
├── GSM5623776_barcodes.tsv.gz
├── GSM5623776_features.tsv.gz
└── GSM5623776_matrix.mtx.gz
```

整理后：

```
.
├── GSM5623768
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── GSM5623769
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── GSM5623770
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── GSM5623776
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
└── prepare.R
```

所用代码：

```
samples <- gsub("_barcodes.tsv.gz","", list.files(pattern = "barcodes.tsv.gz$"))
samples
for (i in samples) {
  dirname <- gsub("GSM\\d*_", "",i)
  dir.create(dirname )
  s.samples <- grep(i, list.files(pattern = "gz$"),value = TRUE)
  for (j in s.samples) {
    file.copy(j, paste0(dirname,"/",gsub("GSM\\d*_","",j)))
  }
  # 对于旧版数据，需要修改文件名后缀
  # file.copy("GSM\\d*_genes.tsv.gz", paste0(dirname,"/","features.tsv.gz"))
}
```

