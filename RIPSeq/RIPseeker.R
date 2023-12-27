library("RIPSeeker")
library("biomaRt")
library("curl")
library("org.Dm.eg.db")

outDir <- file.path(getwd(), "TEST")
# Parameters setting
binSize <- NULL     
minBinSize <- NULL
maxBinSize <- NULL
multicore <- TRUE # TRUE will not run properly
strandType <- NULL # "-"   

biomart <- "ensembl"
biomaRt_dataset <- "dmelanogaster_gene_ensembl"
goAnno <- "org.Dm.eg.db"                    

bamfiles <- list.files(".",pattern = "filter.bam$",full.names = TRUE)
bamfiles
mainSeek <- edit(mainSeek)
# line :suppressMessages(nbhGRList <- mclapply(as.list(split(alignGR, 
#                                                     seqnames(alignGR))), mainSeekSingleChrom, runViterbi = runViterbi,
# add: mc.cores = 24,
seekOut.PRC2 <- ripSeek(bamPath = bamfiles, 
                        cNAME = "SRR408617",
                        # reverseComplement = TRUE, 
                        paired = FALSE,
                        # strandType = strandType,
                        uniqueHit = TRUE, 
                        assignMultihits = TRUE,
                        rerunWithDisambiguatedMultihits = TRUE,
                        binSize=binSize, 
                        minBinSize = minBinSize,
                        maxBinSize = maxBinSize,
                        biomart = biomart,
                        biomaRt_dataset = biomaRt_dataset,
                        multicore=multicore, 
                        outDir=outDir)
# if bad results, try set inBinSize and maxBinSize as NULL
