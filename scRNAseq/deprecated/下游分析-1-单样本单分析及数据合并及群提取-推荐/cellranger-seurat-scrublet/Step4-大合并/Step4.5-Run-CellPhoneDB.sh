#!/bin/bash

source ~/.bashrc
conda activate py37

cd ./4.5-Ligand-Receptor-Analysis/CellPhoneDB
cellphonedb method statistical-analysis *cell.meta.txt *normalized.UMI.txt \
--threads=32 --counts-data gene_name

cellphonedb plot dot-plot
cellphonedb plot heatmap-plot *cell.meta.txt

# cellphonedb plot dot-plot --columns ../APC-columns.txt --output-name APC.Hub.pdf
# cellphonedb plot dot-plot --columns ../T-columns.txt --output-name T.Hub.pdf
conda deactivate
