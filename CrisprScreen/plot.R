library(MAGeCKFlute)
library(ggplot2)
####################### plot the sgRNAs ##################
file2 = "../test/day14-vs-day0.sgrna_summary.txt"
sdata = ReadsgRRA(file2)
sdata

# remove "RPL35" and "DTWD1"
sdata <- sdata[!sdata$Gene %in% c("RPL35","DTWD1"),]


# edit(sgRankView)

df = sdata
gene = NULL
top = 10
bottom = 0
neg_ctrl = NULL
binwidth = 0.3
interval = 0.1
bg.col = NA
filename = NULL 
width = 5
height = 3.5
neg_ctrl = unique(sdata$Gene)



if (!all(c("sgrna", "Gene", "LFC") %in% colnames(df))) 
  stop("Make sure your data contains columns of 'sgrna', 'Gene', and 'LFC' ...")
df = as.data.frame(df, stringsAsFactors = FALSE)
df = df[order(df$LFC), ]
df$Gene = as.character(df$Gene)
tmp = stats::aggregate(df$LFC, by = list(df$Gene), median)
colnames(tmp) = c("Gene", "mid")
tmp = tmp[order(tmp$mid), ]
if (top > 0) {
  idx = max((nrow(tmp) - top + 1), 1)
  gene = c(gene, tmp$Gene[idx:nrow(tmp)])
  
}
if (bottom > 0) {
  gene = c(gene, tmp$Gene[1:min(bottom, nrow(tmp))])
}
gene = unique(gene)


subdf = df[df$Gene %in% gene, ]
if (nrow(subdf) < 2) 
  return(ggplot())
subdf$Gene = factor(subdf$Gene, levels = gene)
subdf = subdf[order(subdf$Gene), ]
subdf$index = rep(1:length(gene), as.numeric(table(subdf$Gene)[gene]))
subdf$yend <- (binwidth + interval) * subdf$index - interval
subdf$y <- (binwidth + interval) * (subdf$index - 1)
color <- c(rep("pos", dim(subdf)[1]))
color[which(subdf[, 3] < 0)] <- "neg"
subdf$color <- color
subdf = subdf[, c("sgrna", "Gene", "LFC", "y", "yend", "color", 
                  "index")]
a <- -Inf
b <- Inf
if (is.na(bg.col)) {
  bg.col <- "white"
}
bindex <- as.vector(sapply(seq(1, max(subdf$index), 1), 
                           function(x) {
                             rep(x, 4)
                           }))
bgcol <- data.frame(as.vector(bindex))
bgcol$color <- c(rep("bg", length(bindex)))
colnames(bgcol) <- c("id", "value")
bgcol$x <- rep(c(a, b, b, a), max(subdf$index))
bgcol$y <- as.vector(sapply(seq(1, max(subdf$index), 1), 
                            function(x) {
                              c((interval + binwidth) * (x - 1), (interval + binwidth) * 
                                  (x - 1), (interval + binwidth) * x - interval, 
                                (interval + binwidth) * x - interval)
                            }))
neg_ctrl <- unique(df$Gene)
if (!is.null(neg_ctrl)) {
  neggene = rep(df[df$Gene %in% neg_ctrl, "Gene"], max(subdf$index))
  negsgrna = rep(df[df$Gene %in% neg_ctrl, "sgrna"], max(subdf$index))
  background <- data.frame(sgrna = as.vector(negsgrna), 
                           Gene = as.vector(neggene))
  background$LFC <- rep(df[df$Gene %in% neg_ctrl, 3], 
                        max(subdf$index))
  seq <- as.vector(sapply(seq(1, max(subdf$index), 1), 
                          function(x) {
                            rep(x, length(df[df$Gene %in% neg_ctrl, 2]))
                          }))
  background$y <- (binwidth + interval) * (seq - 1)
  background$yend <- (binwidth + interval) * seq - interval
  background$color <- rep("tbg", length(neggene))
  background$index = 0
}
cols <- c(pos = "#e41a1c", neg = "#377eb8", tbg = 608, black = "black")
p = ggplot()
p = p + geom_polygon(aes_string("x", "y", fill = "value",group = "id"), color = "white", data = bgcol)
if (!is.null(neg_ctrl)) 
  p = p + geom_segment(aes_string("LFC", "y", xend = "LFC",yend = "yend", color = "color"),alpha = 0.1, data = background) # new added: alpha = 0.1
p = p + geom_segment(aes_string("LFC", "y", xend = "LFC", yend = "yend", color = "color"),size = 0.8, data = subdf)
p = p + scale_color_manual(values = cols)
p = p + scale_fill_manual(values = c(bg = bg.col))
p = p + scale_y_continuous(breaks = bgcol$y[seq(1, nrow(bgcol),4)] + binwidth/2, labels = gene, expand = c(0, 0))
p = p + labs(x = "Log2(Fold change)", y = NULL)
p = p + theme_bw(base_size = 14)
# p = p + geom_vline(xintercept = 0,color = "green",linetype = "dashed",size = 1) # new added
p = p + theme(plot.title = element_text(hjust = 0.5))
p = p + theme(legend.position = "none")
p = p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_blank(), panel.background = element_blank())
p1 <- p

####################### plot the ranked genes ##########
file1 = "../test/day14-vs-day0.gene_summary.txt"
gdata = ReadRRA(file1)
head(gdata)
gdata <- gdata[!gdata$id %in% c("RPL35","DTWD1"),]


gdata$Rank = rank(gdata$Score)
gdata <- gdata[gdata$Score > 0,]

data <- gdata
x = "Rank"
y = "Score"
label = "id"
top = 10
auto_cut_y = TRUE
ylab = "Log2FC"
groups = c("top", "bottom")
model = c("none","ninesquare", "volcano", "rank")[1]
x_cut = NULL
y_cut = NULL
slope = 1
intercept = NULL
auto_cut = FALSE
auto_cut_x = auto_cut
auto_cut_diag = auto_cut
group_col = NULL
groupnames = NULL
label.top = TRUE
toplabels = NULL
display_cut = FALSE
color = NULL
shape = 16
size = 1
alpha = 0.6
main = NULL
xlab = x
legend.position = "none"

requireNamespace("ggplot2", quietly = TRUE) || stop("need ggplot package")
requireNamespace("ggrepel", quietly = TRUE) || stop("need ggrepel package")
data = as.data.frame(data, stringsAsFactors = FALSE)
data = data[!(is.na(data[, x]) | is.na(data[, y])), ]
if (label == 0) {
  data$Label = rownames(data)
}else{ data$Label = as.character(data[, label]) }

model = tolower(model)
if (model == "ninesquare") {
  if (length(x_cut) == 0) 
    x_cut = c(-CutoffCalling(data[, x], 2), CutoffCalling(data[, 
                                                               x], 2))
  if (length(y_cut) == 0) 
    y_cut = c(-CutoffCalling(data[, y], 2), CutoffCalling(data[, 
                                                               y], 2))
  if (length(intercept) == 0) 
    intercept = c(-CutoffCalling(data[, y] - slope * 
                                   data[, x], 2), CutoffCalling(data[, y] - slope * 
                                                                  data[, x], 2))
}
if (model == "volcano") {
  if (length(x_cut) == 0) 
    x_cut = c(-CutoffCalling(data[, x], 2), CutoffCalling(data[, 
                                                               x], 2))
  if (length(y_cut) == 0) 
    y_cut = -log10(0.05)
}
if (model == "rank") {
  if (length(x_cut) == 0) 
    x_cut = c(-CutoffCalling(data[, x], 2), CutoffCalling(data[, 
                                                               x], 2))
}
if (auto_cut_x) 
  x_cut = c(-CutoffCalling(data[, x], auto_cut_x), CutoffCalling(data[, 
                                                                      x], auto_cut_x))
if (auto_cut_y) 
  y_cut = c(-CutoffCalling(data[, y], auto_cut_y), CutoffCalling(data[, 
                                                                      y], auto_cut_y))
if (auto_cut_diag) 
  intercept = c(-CutoffCalling(data[, y] - slope * data[, 
                                                        x], auto_cut_diag), CutoffCalling(data[, y] - slope * 
                                                                                            data[, x], auto_cut_diag))
avail_groups = c("topleft", "topright", "bottomleft", "bottomright", 
                 "midleft", "topcenter", "midright", "bottomcenter", 
                 "midcenter", "top", "mid", "bottom", "left", "center", 
                 "right", "none")
mycolour = c("#1f78b4", "#fb8072", "#33a02c", "#ff7f00", 
             "#bc80bd", "#66c2a5", "#6a3d9a", "#fdb462", "#ffed6f", 
             "#e78ac3", "#fdb462", "#8da0cb", "#66c2a5", "#fccde5", 
             "#fc8d62", "#d9d9d9")
names(mycolour) = avail_groups
if (model == "ninesquare") 
  groups = c("midleft", "topcenter", "midright", "bottomcenter")
if (model == "volcano") 
  groups = c("topleft", "topright")
if (model == "rank") 
  groups = c("left", "right")
groups = intersect(groups, avail_groups)
if (length(x_cut) > 0) {
  idx1 = data[, x] < min(x_cut)
  idx2 = data[, x] > max(x_cut)
}else {
  idx1 = NA
  idx2 = NA
}
if (length(y_cut) > 0) {
  idx3 = data[, y] < min(y_cut)
  idx4 = data[, y] > max(y_cut)
}else {
  idx3 = NA
  idx4 = NA
}
if (length(intercept) > 0) {
  idx5 = data[, y] < slope * data[, x] + min(intercept)
  idx6 = data[, y] > slope * data[, x] + max(intercept)
}else {
  idx5 = NA
  idx6 = NA
}
data$group = "none"
for (gr in groups) {
  if (gr == "topleft") 
    idx = cbind(idx1, idx4, idx6)
  if (gr == "topcenter") 
    idx = cbind(!idx1, !idx2, idx4, idx6)
  if (gr == "topright") 
    idx = cbind(idx2, idx4, idx6)
  if (gr == "midleft") 
    idx = cbind(idx1, idx6, !idx3, !idx4)
  if (gr == "midcenter") 
    idx = cbind(!idx1, !idx2, !idx3, !idx4, !idx5, !idx6)
  if (gr == "midright") 
    idx = cbind(idx2, !idx3, !idx4, idx5)
  if (gr == "bottomleft") 
    idx = cbind(idx1, idx3, idx5)
  if (gr == "bottomcenter") 
    idx = cbind(!idx1, !idx2, idx3, idx5)
  if (gr == "bottomright") 
    idx = cbind(idx2, idx3, idx5)
  if (gr == "top") {
    if (length(y_cut) > 0 & length(intercept) > 0) 
      idx = idx4 & idx6
    else if (length(y_cut) > 0) 
      idx = idx4
    else idx = idx6
  }
  if (gr == "mid") 
    idx = (!idx3) & (!idx4)
  if (gr == "bottom") {
    if (length(y_cut) > 0 & length(intercept) > 0) 
      idx = idx3 & idx5
    else if (length(y_cut) > 0) 
      idx = idx3
    else idx = idx5
  }
  if (gr == "left") {
    if (length(x_cut) > 0 & length(intercept) > 0) 
      if (slope > 0) 
        idx = idx1 & idx6
    else idx = idx1 & idx5
    else if (length(x_cut) > 0) 
      idx = idx1
    else if (slope > 0) 
      idx = idx6
    else idx = idx5
  }
  if (gr == "center") 
    idx = (!idx1) & (!idx2)
  if (gr == "right") {
    if (length(x_cut) > 0 & length(intercept) > 0) 
      if (slope > 0) 
        idx = idx2 & idx5
    else idx = idx2 & idx6
    else if (length(x_cut) > 0) 
      idx = idx2
    else if (slope > 0) 
      idx = idx5
    else idx = idx6
  }
  if (is.null(ncol(idx))) {
    if (sum(!is.na(idx)) > 0) 
      data$group[idx] = gr
    else warning("No cutpoint for group:", gr)
  }
  else {
    idx = idx[, !is.na(idx[1, ])]
    if (is.null(ncol(idx))) 
      warning("No cutpoint for group:", gr)
    else if (ncol(idx) < 4 & gr == "midcenter") 
      warning("No cutpoint for group:", gr)
    else data$group[rowSums(idx) == ncol(idx)] = gr
  }
}
data$group = factor(data$group, levels = unique(c(groups, 
                                                  "none")))
if (length(groupnames) != length(groups)) 
  groupnames = groups
if (length(groups) > 0) 
  names(groupnames) = groups
if (length(group_col) == length(groups)) 
  mycolour[groups] = group_col
if (length(groups) == 0) 
  mycolour["none"] = "#FF6F61"
data$rank = top + 1
for (g in groups) {
  idx1 = data$group == g
  x_symb = 0
  y_symb = 0
  if (g == "topleft") {
    x_symb = 1
    y_symb = -1
  }
  if (g == "topcenter") {
    x_symb = 0
    y_symb = -1
  }
  if (g == "topright") {
    x_symb = -1
    y_symb = -1
  }
  if (g == "midleft") {
    x_symb = 1
    y_symb = 0
  }
  if (g == "midright") {
    x_symb = -1
    y_symb = 0
  }
  if (g == "bottomleft") {
    x_symb = 1
    y_symb = 1
  }
  if (g == "bottomcenter") {
    x_symb = 0
    y_symb = 1
  }
  if (g == "bottomright") {
    x_symb = -1
    y_symb = 1
  }
  if (g == "top") {
    x_symb = 0
    y_symb = -1
  }
  if (g == "bottom") {
    x_symb = 0
    y_symb = 1
  }
  if (g == "left") {
    x_symb = 1
    y_symb = 0
  }
  if (g == "right") {
    x_symb = -1
    y_symb = 0
  }
  tmp = data[, c(x, y)]
  tmp[, x] = (tmp[, x] - min(tmp[, x]))/(max(tmp[, x]) - 
                                           min(tmp[, x]))
  tmp[, y] = (tmp[, y] - min(tmp[, y]))/(max(tmp[, y]) - 
                                           min(tmp[, y]))
  data$rank[idx1] = rank((x_symb * tmp[, x] + y_symb * 
                            tmp[, y])[idx1])
}
data$rank[data$rank == 0] = Inf
if (mode(toplabels) == "list") {
  data$Label[data$rank > top & !(data$Label %in% unlist(toplabels))] = ""
  data$group = data$Label
  if (length(toplabels) > 0) {
    tmp = stack(toplabels)
    tmp = tmp[!duplicated(tmp[, 1]), ]
    rownames(tmp) = tmp[, 1]
    data$group[data$group %in% tmp[, 1]] = as.character(tmp[data$group[data$group %in% 
                                                                         tmp[, 1]], 2])
    data$group[!(data$group %in% tmp[, 2]) & data$group != 
                 ""] = "Top hits"
  }
}else {
  data$Label[data$rank > top & !(data$Label %in% toplabels)] = ""
}
if (is.null(color)) {
  color = "group"
}else if (length(color) == 1) {
  if (!color %in% colnames(data)) {
    data$color = color
    color = "color"
  }
}else {
  data$color = color[1]
  color = "color"
  warning("Only the first color is took.")
}
gg = data
gg = gg[order(gg[, color]), ]
gg$Rank <- max(gg$Rank) - gg$Rank + 1 # reverse the point direction

p = ggplot(gg, aes_string(x, y, label = "Label"))
if (all(c(shape, size) %in% colnames(gg))) {
  p = p + geom_point(aes_string(shape = shape, size = size),alpha = alpha)
}else if (shape %in% colnames(gg)) {
  p = p + geom_point(aes_string(shape = shape), size = size, alpha = alpha)
}else if (size %in% colnames(gg)) {
  p = p + geom_point(aes_string(size = size), shape = shape,  alpha = alpha)
}else {
  p = p + geom_point(size = size, shape = shape, alpha = alpha,color = "gray80")
}
p = p + geom_point(data = gg[gg$id %in% gene,,],size = size, shape = shape, alpha = alpha,color = "red")
# if (color == "group") {
#   if (mode(toplabels) != "list") 
#     p = p + scale_color_manual(values = mycolour[names(groupnames)], labels = groupnames)
#   else p = p + scale_color_manual(values = c("#d9d9d9", "#fb8072", "#80b1d3", "#fdb462", "#bc80bd", "#b3de69", 
#                                              "#bebada", "#8dd3c7", "#ffffb3", "#fccde5", "#ccebc5", "#ffed6f"))
# }else {
#   if (mode(gg[, color]) == "numeric") 
#     p = p + scale_color_gradient2(low = "#377eb8", high = "#e41a1c", midpoint = 0)
#   else if (!"try-error" %in% class(try(col2rgb(gg[1, color]), silent = TRUE))) {
#     mycolour = unique(gg[, color])
#     names(mycolour) = mycolour
#     p = p + scale_color_manual(values = mycolour)
#   }
#   else {
#     p = p + scale_color_brewer(type = "div")
#   }
# }
if (label.top) 
  p = p + ggrepel::geom_text_repel(max.overlaps = 20)
if (display_cut) {
  if (length(x_cut) > 0) 
    p = p + geom_vline(xintercept = x_cut, linetype = "dotted")
  if (length(y_cut) > 0) 
    p = p + geom_hline(yintercept = y_cut, linetype = "dotted")
  if (length(intercept) > 0) 
    p = p + geom_abline(slope = slope, intercept = intercept, 
                        linetype = "dotted")
}
p = p + labs(x = xlab, y = ylab, title = main, color = NULL)
p = p + theme_bw(base_size = 14)
p = p + theme(plot.title = element_text(hjust = 0.5))
p = p + theme(legend.position = legend.position)
p +
  theme(panel.grid =element_blank()) +
  theme(panel.border = element_blank()) +   ## 删去外层边框
  theme(axis.line = element_line(size=1, colour = "black")) -> p2
p2
######## merge two plots ###############
library(patchwork)
wrap_plots(p1,p2)

