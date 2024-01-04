# Fly

## build-fly-reference-from-ensembl.sh

***Recommended!***

关于为什么不使用从Flybase下载的最新reference，下面对两个数据库的gtf文件做了一个对比。需要注意的时Cellranger要求使用protein coding genes! (猜测是能表达出mRNA的基因【flybase】，或者biotype为”protein_coding“的【Ensembl】)

对于Ensembl-release-110:

```bash
$ grep -v "^#" BDGP6.46.110.gtf | cut -f 3 | sort | uniq -c
 163189 CDS
 196625 exon
  46766 five_prime_utr
  24278 gene
      4 Selenocysteine
  30802 start_codon
  30791 stop_codon
  33737 three_prime_utr
  41610 transcript
  
$ grep -v "^#" BDGP6.46.110.gtf | grep "transcript_id" | cut -f 9 | cut -d " " -f 2 | cut -d '"' -f 2 | sort | uniq | wc -l
24278

$ grep -v "^#" BDGP6.46.110.gtf | awk '{if ( $3 == "gene") print $0;}' | grep "protein_coding" | wc -l
13986

$ grep -v "^#" BDGP6.46.110.filtered.gtf | cut -f 3 | sort | uniq -c
 163189 CDS
 183572 exon
  46766 five_prime_utr
  13986 gene
      4 Selenocysteine
  30802 start_codon
  30791 stop_codon
  33737 three_prime_utr
  30789 transcript
```

对于Flybase r6.55:

```
$ cut -f 3  dmel-all-r6.55.gtf | sort | uniq -c
  33741 3UTR
  46806 5UTR
 163253 CDS
 190039 exon
  17873 gene
    485 miRNA
  30802 mRNA
   3060 ncRNA
    262 pre_miRNA
    365 pseudogene
    115 rRNA
    270 snoRNA
     32 snRNA
  30888 start_codon
  30828 stop_codon
    312 tRNA

$ cat dmel-all-r6.55.gtf | grep "transcript_id" | cut -f 9 | cut -d " " -f 2 | cut -d '"' -f 2 | sort | uniq | wc -l
17872

$ cat dmel-all-r6.55.gtf | awk '{if ($3 == "mRNA") print $0;}' | cut -f 9 | cut -d " " -f 2 | cut -d '"' -f 2 | sort | uniq -c | wc -l
13986

```

**结论，flybase中这13986个有mRNA的基因很可能对应Ensembl-release-110中过滤后的基因！所以用Ensembl的gtf就可以了！**

## build-fly-reference-from-flybase.sh

脚本会使用所有基因，包括非蛋白编码基因，这是**cellranger所不推荐的**！



# Human

## build-human-reference-from-ensembl.sh

**version**: Ensembl-Release110-Genecode-v44

```bash
$ cat gencode.v44.primary_assembly.annotation.gtf.filtered | grep -v "^#" | cut -f 3 | sort | uniq -c 
 876939 CDS
1587066 exon
  38592 gene
    130 Selenocysteine
  96996 start_codon
  90699 stop_codon
 226049 transcript
 376802 UTR
```



# Mouse

## build-mouse-reference-from-ensembl.sh

***Recommended!***

```bash
$ cat gencode.vM33.primary_assembly.annotation.gtf.filtered | grep -v "^#" | cut -f 3 | sort | uniq -c 
 524683 CDS
 833846 exon
  33611 gene
     64 Selenocysteine
  59524 start_codon
  55248 stop_codon
 125424 transcript
 184727 UTR
```



## build-mouse-reference-from-ensembl-with-pseudogenes.sh

```bash
$ cat gencode.vM33.primary_assembly.annotation.gtf.filtered | grep -v "^#" | cut -f 3 | sort | uniq -c 
 524683 CDS
 855852 exon
  47150 gene
     64 Selenocysteine
  59524 start_codon
  55248 stop_codon
 139519 transcript
 184727 UTR
```



