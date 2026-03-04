#!/bin/bash
source ~/.bashrc
conda activate scanpy-latest-py38

jupyter nbconvert --to notebook --inplace --execute Step4.4-scVelo-TSNE.ipynb --ExecutePreprocessor.timeout=3000
jupyter nbconvert --to html Step4.4-scVelo-TSNE.ipynb

jupyter nbconvert --to notebook --inplace --execute Step4.4-scVelo-UMAP.ipynb --ExecutePreprocessor.timeout=3000
jupyter nbconvert --to html Step4.4-scVelo-UMAP.ipynb

conda deactivate
