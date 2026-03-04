# 1. Install and Load NOISeq
if (!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
BiocManager::install("NOISeq")
library(NOISeq)

# 2. Prepare Your Data
# Assume 'counts' is your matrix (genes x samples) and 'factors' defines groups
# Example: counts <- matrix(rpois(1000), 100, 10); colnames(counts) <- paste0("S", 1:10)
# Example: factors <- data.frame(condition = rep(c("Control", "Treat"), 5))
counts <- read.table("../6_Counts_featureCounts/final_counts.txt",stringsAsFactors = FALSE,row.names = 1, header = TRUE)
counts <- counts[,6:7]
colnames(counts) <- c('W1118', 'Brat-RNAi')
sample_group <- data.frame(group = c('control', 'case'), row.names = colnames(counts))


# 3. Read Data into NOISeq format (eSet)
# You might need feature length/biotype for full QC, but for DE with no reps, counts/factors are key
mydata <- readData(counts, factors = sample_group)

# 4. Run NOISeq (No Replicates)
# Specify 'replicates = "no"' and set simulation parameters (e.g., nss=5, pnr=0.2)
noiseq_results <- noiseq(mydata, replicates = "no", k = 0.5, pnr = 0.2, nss = 5, factor = 'group')
write.csv(noiseq_results@results[[1]], file = 'NOISeq.csv')


# 5. Get Differentially Expressed Genes (DEGs)
# Filter for significant genes (e.g., with a certain reliability/probability)
DEGs <- degenes(noiseq_results, q = 0.95) # Adjust q for desired confidence

# 6. Visualize Results (e.g., Volcano Plot, MA Plot)
plot(noiseq_results) # Basic plot
plotMA(noiseq_results, n=20) # MA plot
