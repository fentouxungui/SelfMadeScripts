library(circlize)
cytoband.df = read.table("../../../../CircosDemo/cytoBandIdeo.txt", colClasses = c("character", "numeric",
                                                       "numeric", "character", "character"), sep = "\t")
head(cytoband.df)
cytoband.df <- cytoband.df[!grepl("NW",cytoband.df$V1),]
circos.clear()
circos.par(gap.after = c(rep(1, length(unique(cytoband.df$V1)) -1 ), 5))
# circos.par("gap.degree" = rep(3, length(unique(cytoband.df$V1))))
circos.initializeWithIdeogram(cytoband.df, plotType = c("axis", "labels"))
# circos.track(ylim = c(0, 1), 
#              bg.col = rep_len(c("#FF000040", "#00FF0040", "#0000FF40"),length((unique(cytoband.df$V1)))), 
#              bg.border = NA, track.height = 0.05)

# set.seed(123)
# circos.initializeWithIdeogram(cytoband.df, plotType = NULL)
# circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
#   chr = CELL_META$sector.index
#   xlim = CELL_META$xlim
#   ylim = CELL_META$ylim
#   circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))
#   circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "white",
#               facing = "inside", niceFacing = TRUE)
# }, track.height = 0.15, bg.border = NA)

# read.cytoband.2(species = "galGal6")
test <- read.delim("./D1.D7_DMR_0.5.anno.withMeanMethylationLevel.txt",stringsAsFactors = FALSE)
head(test)
bed <- data.frame(chr = test$seqnames,
                  start = test$start,
                  end = test$end,
                  value1 = -log10(test$avg_qvalue) * ifelse(test$case.control == "positive",1,-1))
# bed.up <- bed[bed$value1 > 0,]
# bed.down <- bed[bed$value1 < 0,]
# circos.genomicDensity(bed.up, col = c("red"), track.height = 0.1)
circos.genomicTrack(bed,bg.col = "#e2f2fa", 
                    panel.fun = function(region, value, ...) {
                        circos.lines(CELL_META$cell.xlim, c(-20,-20), lty = 3, col = "#AAAAAA")
                        circos.lines(CELL_META$cell.xlim, c(-10,-10), lty = 3, col = "#AAAAAA")
                        circos.lines(CELL_META$cell.xlim, c(10,10), lty = 3, col = "#AAAAAA")
                        circos.lines(CELL_META$cell.xlim, c(20,20), lty = 3, col = "#AAAAAA")
                        circos.genomicPoints(region, value, pch = 16, cex = 0.3, 
                                             col = ifelse(value[[1]] > 0, "#E41A1C", "#377EB8"))
                    })
circos.yaxis(side = "left", at =  c(-20,-10,0,10,20),sector.index = get.all.sector.index()[1], labels.cex = 0.8)

test <- read.delim("./D1.D63_DMR_0.5.anno.withMeanMethylationLevel.txt",stringsAsFactors = FALSE)
bed <- data.frame(chr = test$seqnames,
                  start = test$start,
                  end = test$end,
                  value1 = -log10(test$avg_qvalue) * ifelse(test$case.control == "positive",1,-1))
# bed.up <- bed[bed$value1 > 0,]
# bed.down <- bed[bed$value1 < 0,]
# circos.genomicDensity(bed.up, col = c("red"), track.height = 0.1)
circos.genomicTrack(bed,bg.col = "#fee1ee", 
                    panel.fun = function(region, value, ...) {
                      circos.lines(CELL_META$cell.xlim, c(-20,-20), lty = 3, col = "#AAAAAA")
                      circos.lines(CELL_META$cell.xlim, c(-10,-10), lty = 3, col = "#AAAAAA")
                      circos.lines(CELL_META$cell.xlim, c(10,10), lty = 3, col = "#AAAAAA")
                      circos.lines(CELL_META$cell.xlim, c(20,20), lty = 3, col = "#AAAAAA")
                      circos.genomicPoints(region, value, pch = 16, cex = 0.3, 
                                           col = ifelse(value[[1]] > 0, "#E41A1C", "#377EB8"))
                    })
# circos.yaxis(side = "left", at =  c(-20,-10,0,10,20),sector.index = get.all.sector.index()[1], labels.cex = 0.4)

test <- read.delim("./D7.D63_DMR_0.5.anno.withMeanMethylationLevel.txt",stringsAsFactors = FALSE)
bed <- data.frame(chr = test$seqnames,
                  start = test$start,
                  end = test$end,
                  value1 = -log10(test$avg_qvalue) * ifelse(test$case.control == "positive",1,-1))
# bed.up <- bed[bed$value1 > 0,]
# bed.down <- bed[bed$value1 < 0,]
# circos.genomicDensity(bed.up, col = c("red"), track.height = 0.1)
circos.genomicTrack(bed,bg.col = "#d7e3f9",
                    panel.fun = function(region, value, ...) {
                      circos.lines(CELL_META$cell.xlim, c(-20,-20), lty = 3, col = "#AAAAAA")
                      circos.lines(CELL_META$cell.xlim, c(-10,-10), lty = 3, col = "#AAAAAA")
                      circos.lines(CELL_META$cell.xlim, c(10,10), lty = 3, col = "#AAAAAA")
                      circos.lines(CELL_META$cell.xlim, c(20,20), lty = 3, col = "#AAAAAA")
                      circos.genomicPoints(region, value, pch = 16, cex = 0.3, 
                                           col = ifelse(value[[1]] > 0, "#E41A1C", "#377EB8"))
                    })
# circos.yaxis(side = "left", at =  c(-20,-10,0,10,20),sector.index = get.all.sector.index()[1], labels.cex = 0.4)
# circos.genomicDensity(bed, col = c("green"), track.height = 0.1)

# control <- data.frame(chr = test$seqnames,
#                   start = test$start,
#                   end = test$end,
#                   value1 = test$avg.g1.methyl)
# case <- data.frame(chr = test$seqnames,
#                       start = test$start,
#                       end = test$end,
#                       value1 = test$avg.g2.methyl)
# circos.genomicTrack(control,
#                     panel.fun = function(region, value, ...) {
#                       circos.genomicLines(region, value, ...)
#                    })
text(0, 0, "DMR Circular plot", cex = 1)
circos.clear()


