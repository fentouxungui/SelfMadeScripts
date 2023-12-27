# %matplotlib inline
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import re

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

def Scrublet_Rmd(mat_dir="./3.1-QC-Passed-Cells-Rds/"):
    for Acounts_matrix in os.listdir(mat_dir):
        if Acounts_matrix.endswith('UMI.mtx'):
            print("***** " + Acounts_matrix + " *****")
            counts_matrix = scipy.io.mmread(mat_dir + Acounts_matrix).T.tocsc()
            scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
            doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,
                                                                      min_cells=3,
                                                                      min_gene_variability_pctl=85,
                                                                      n_prin_comps=30)
            scrub.plot_histogram()
            plt.savefig("scrublet.score.png")
    return doublet_scores, predicted_doublets
