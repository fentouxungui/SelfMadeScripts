library(edgeR)
counts <- read.table("../6_Counts_featureCounts/final_counts.txt",stringsAsFactors = FALSE,row.names = 1, header = TRUE)
counts <- counts[,6:7]
colnames(counts) <- c('W1118', 'Brat-RNAi')
# Assuming 'counts' is your matrix and 'y' is DGEList object
y <- DGEList(counts=counts, group=c("Control", "Treated")) # Even with 1 sample/group
y <- calcNormFactors(y)
# Estimate dispersion (adjust 'method' and 'robust' as needed)
# y <- estimateCommonDisp(y, verbose=TRUE) # Use this if you have *some* replicates
# If *no* replicates, you might use Poisson or a fixed value:
# y <- estimateGLMCommonDisp(y, design=matrix(1, ncol=2, nrow=2), robust=TRUE) # Example for 2 samples
# et <- exactTest(y, dispersion=y$common.dispersion^2) # Or set manually: dispersion=0.04^2
et <- exactTest(y, dispersion=0.2^2)
write.csv(et,file = 'edgeR-bcv0.2.csv')
