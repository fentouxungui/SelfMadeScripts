library(ggplot2)
library(tidyverse)
library(ggforce)
library(ggprism)

library(qqman)

peaks.files <- list.files("./",pattern = "txt$")
peaks.files
i <- peaks.files[3]
i
temp <- read.delim(i,stringsAsFactors = FALSE)
temp <- temp[,c("seqnames","start","end","name","score","strand","annotation","q_val", "meth_diff","SYMBOL", "transcriptId", "GENENAME")]
# temp$start <- temp$start - 1
temp$annotation <- gsub("\\(.*\\)$","",temp$annotation)
temp <- temp[temp$annotation != "Distal Intergenic",]
temp <- temp[!(grepl("^chrUn",temp$seqnames)),]
temp$q_val <- -log10(temp$q_val) + log10(0.5) # if use trans = scales::pseudo_log_trans() in ggplot2, this 
# will remove the white area along y=0
temp$score <- temp$q_val
temp$q_val <- sign(temp$meth_diff) * temp$q_val
temp$seqnames <- gsub("chr","",temp$seqnames)
temp$seqnames <- factor(temp$seqnames,levels =  c(1:28,30:33,"Z"))
###数据处理
head(temp)
temp %>% 
  group_by(seqnames) %>% 
  summarise(df_chr_len=max(end)) %>% #计算染色体长度
  mutate(total = cumsum(df_chr_len) - df_chr_len + sum(df_chr_len) * 0.01 * 0:(length(df_chr_len)-1)) %>%
  select(-df_chr_len) %>% #计算染色体初始位置
  left_join(temp, ., by="seqnames")  %>%
  arrange(seqnames, start) %>%
  mutate( BPcum = end + total) -> df_SNP_position#计算累计SNP的位置


X_axis <-  df_SNP_position %>% group_by(seqnames) %>% summarize(center=( max(BPcum) +min(BPcum) ) / 2 )
X_axis

#添加高亮和注释信息：snpsOfInterest中的rs编号和P值大于10的点
data <- df_SNP_position %>%
  mutate( is_hyper=ifelse(q_val > 2, "yes", "no")) %>%
  mutate( is_hypo=ifelse(q_val < -2, "yes", "no")) %>%
  mutate( is_annotate=ifelse(abs(q_val) > 30, "yes", "no"))
#绘图
ggplot(data, aes(x=BPcum, y=q_val)) +
  geom_point(aes(color=seqnames),alpha=0.8, size=1)+
  # scale_color_manual(values = rep(c('#30A9DE','#EFDC05','#E53A40','#090707'), length(levels(data$seqnames)) ))+#颜色设置
  scale_color_manual(values = rep(c('gray',"darkgray"), length(levels(data$seqnames)) ))+#颜色设置
  scale_x_continuous(label = X_axis$seqnames, breaks= X_axis$center)+#设定X轴
  scale_y_continuous(trans = scales::pseudo_log_trans(), 
                     breaks = c(-20,-10,-6,-4,-2,0,2,4,6,10,20)) +#去除绘图区和X轴之间的gap
  geom_hline(yintercept = c(-2,2), #添加阈值线
             color = 'darkgray',size = 1, linetype = "dashed") + 
  geom_point(data=subset(data, is_hyper=="yes"), color="red", size=1.5)+
  geom_point(data=subset(data, is_hypo=="yes"), color="green",  size=1.5)+
  theme(legend.position = "none",#去除图例
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line()) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 90)) + 
  labs(y= "-log10(q value)", x = "Chrosome")


############# part2

#添加高亮和注释信息：snpsOfInterest中的rs编号和P值大于10的点
data <- df_SNP_position %>%
  mutate( is_hyper=ifelse(q_val > 0, "yes", "no")) %>%
  mutate( is_annotate=ifelse(abs(q_val) > 30, "yes", "no"))
#绘图
ggplot(data, aes(x=BPcum, y=q_val - 0.49)) +
  geom_point(alpha=0.6, size=1,col = "green")+
  scale_x_continuous(label = X_axis$seqnames, breaks= X_axis$center)+#设定X轴
  scale_y_continuous(trans = scales::pseudo_log_trans(), 
                     breaks = c(-20,-10,-6,-4,-2,0,2,4,6,10,20)) +#去除绘图区和X轴之间的gap
  geom_point(data=subset(data, is_hyper=="yes"), color="red", size=1,alpha=0.6)+
  theme(legend.position = "none",#去除图例
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line()) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 90)) + 
  labs(y= "-log10(q value)", x = "Chrosome")

