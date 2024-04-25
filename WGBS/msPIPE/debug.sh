#!/bin/bash
conda activate /home/xilab/software/miniconda-envs/msPIPE
ls | while read id
do
    bismark2bedGraph -o CpG.cov --dir ${id}/methylcontext ${id}/methylcontext/CpG_context_*  --buffer_size 30%
    bismark2bedGraph -o CHG.cov --dir ${id}/methylcontext --CX ${id}/methylcontext/CHG_context_*  --buffer_size 30%
    bismark2bedGraph -o CHH.cov --dir ${id}/methylcontext --CX ${id}/methylcontext/CHH_context_*  --buffer_size 30%
    ## split files
    gzip -dc ${id}/methylcontext/CpG.cov.gz.bismark.cov.gz 1> ${id}/methylcontext/${id}_CpG.cov.txt
    gzip -dc ${id}/methylcontext/CHG.cov.gz.bismark.cov.gz 1> ${id}/methylcontext/${id}_CHG.cov.txt
    gzip -dc ${id}/methylcontext/CHH.cov.gz.bismark.cov.gz 1> ${id}/methylcontext/${id}_CHH.cov.txt
    ~/software/msPIPE/msPIPE/bin/script/splitF_bychr.py ${id}/methylcontext/${id}_CpG.cov.txt ${id}/methylcontext/CpG_chr/${id}
    ~/software/msPIPE/msPIPE/bin/script/splitF_bychr.py ${id}/methylcontext/${id}_CHG.cov.txt ${id}/methylcontext/CHG_chr/${id}
    ~/software/msPIPE/msPIPE/bin/script/splitF_bychr.py ${id}/methylcontext/${id}_CHH.cov.txt ${id}/methylcontext/CHH_chr/${id}
done

# cd /home/xilab/zhangyc/WGBS/chicken/PRJNA389197/msPIPE
# ~/software/msPIPE/msPIPE/msPIPE.py -p params_galGal6.conf -c 12 --skip_calling --calling_data test_result/methylCALL -o test_result/