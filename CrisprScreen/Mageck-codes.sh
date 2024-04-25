# 1. some fastqc commands

ls *.gz | xargs fastqc -o ../fastqc -t 6

# 判定 R2 reads有没有sgRNA信息

zcat NC_L3_Q3344W0052.R2.fastq.gz | awk '{if (NR%4 == 2) print $1}' | head -n 100

# 2. mageck count

mageck count \
-l ./zj.library.csv \
-n all \
--sample-label D6-L,D6-N,L-18h,N-18h,NC \
--fastq \
../fastq/D6-L_L3_Q1878W0194.R1.fastq.gz \
../fastq/D6-N_L2_Q1878W0194.R1.fastq.gz \
../fastq/L-18h_L4_Q1878W0194.R1.fastq.gz \
../fastq/N-18h_L3_Q1878W0194.R1.fastq.gz \
../fastq/NC_L1_Q1878W0194.R1.fastq.gz \
--fastq-2 \
../fastq/D6-L_L3_Q1878W0194.R2.fastq.gz \
../fastq/D6-N_L2_Q1878W0194.R2.fastq.gz \
../fastq/L-18h_L4_Q1878W0194.R2.fastq.gz \
../fastq/N-18h_L3_Q1878W0194.R2.fastq.gz \
../fastq/NC_L1_Q1878W0194.R2.fastq.gz \
--pdf-report \
--trim-5 17,18

# use code bellow #
mageck count \
-l ../zj.library.csv \
-n all \
--sample-label D6-L,D6-N,L-18h,N-18h,NC \
--fastq \
../fastq/D6-L_L3_Q1878W0194.R1.fastq.gz \
../fastq/D6-N_L2_Q1878W0194.R1.fastq.gz \
../fastq/L-18h_L4_Q1878W0194.R1.fastq.gz \
../fastq/N-18h_L3_Q1878W0194.R1.fastq.gz \
../fastq/NC_L1_Q1878W0194.R1.fastq.gz \
--pdf-report \
--trim-5 17,18,19,20


# 3. mageck test
mageck test -k ../count/all.count.txt -t D6-L -c NC  -n D6-L-vs-NC
mageck test -k ../count/all.count.txt -t D6-N -c NC  -n D6-N-vs-NC
mageck test -k ../count/all.count.txt -t L-18h -c NC  -n L-18h-vs-NC
mageck test -k ../count/all.count.txt -t N-18h -c NC  -n N-18h-vs-NC
# or
mageck test -k ../count/all.count.txt --day0-label NC -n All
Rscript All.R

# mageck test -k ../count/all.count.txt -t D6-L,D6-N -c L-18h,N-18h  -n paired --paired # not used!

####################### Not run
# 4. mageck mle

## day28-day14-high: day28-day14-day0-Togethor
### 4.1 designmat.txt
Samples	baseline	day14High	day28high
day0	1	0	0
day14	1	1	0
day28conR3	1	0	1
day28samR4	1	0	1
### 4.2 run
mageck mle -k ../../count/all.count.txt -d designmat.txt -n day14-day28-high

## day28conR3-vs-day0
### 4.1 designmat.txt
Samples	baseline	day28high
day0	10	0
day28conR3	1	1
### 4.2 run
mageck mle -k ../../count/all.count.txt -d designmat.txt -n day28conR3High


## day28samR4-vs-day0
### 4.1 designmat.txt
Samples	baseline	day28high
day0	10	0
day28samR4	1	1
### 4.2 run
mageck mle -k ../../count/all.count.txt -d designmat.txt -n day28samR4High


################ Not run ################
# 5. Vispr visualization

vispr server day14-vs-day0-mle.vispr.yaml

or

vispr server day14-vs-day0-test.vispr.yaml

# 6. MAGeCKFlute (Optional)

check: https://www.bioconductor.org/packages/release/bioc/vignettes/MAGeCKFlute/inst/doc/MAGeCKFlute.html



