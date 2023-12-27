#!/bin/bash
source ~/.bashrc
conda activate py37
jupyter nbconvert --to notebook --inplace --execute Step3.6.1-scVelo-TSNE.ipynb --ExecutePreprocessor.timeout=300
jupyter nbconvert --to html Step3.6.1-scVelo-TSNE.ipynb

jupyter nbconvert --to notebook --inplace --execute Step3.6.2-scVelo-UMAP.ipynb --ExecutePreprocessor.timeout=300
jupyter nbconvert --to html Step3.6.2-scVelo-UMAP.ipynb
conda deactivate
