#!/bin/bash

source ~/.bashrc
conda activate scrublet-toberemoved
jupyter nbconvert --to notebook --inplace --execute scrublet_basics.ipynb --ExecutePreprocessor.timeout=300
jupyter nbconvert --to html scrublet_basics.ipynb
