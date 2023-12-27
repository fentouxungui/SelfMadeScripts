#!/bin/bash
source ~/.bashrc
conda activate py37
# Merge loom file
jupyter nbconvert --to notebook --inplace --execute Step4.4-Merge-Loom-files-of-Replicates.ipynb --ExecutePreprocessor.timeout=600
jupyter nbconvert --to html Step4.4-Merge-Loom-files-of-Replicates.ipynb


# Run scVelo
jupyter nbconvert --to notebook --inplace --execute Step4.4.1-scVelo-TSNE.ipynb --ExecutePreprocessor.timeout=600
jupyter nbconvert --to html Step4.4.1-scVelo-TSNE.ipynb

jupyter nbconvert --to notebook --inplace --execute Step4.4.2-scVelo-UMAP.ipynb --ExecutePreprocessor.timeout=600
jupyter nbconvert --to html Step4.4.2-scVelo-UMAP.ipynb
conda deactivate
