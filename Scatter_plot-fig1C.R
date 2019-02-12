setwd("~/Google Drive/NYUMC/Skok lab/BiSulphiteGenome_analysis/June_2018_re-doing-some-plots")
library(readr)
library(ggplot2)

make_average_meth_plot <- function (file1, file2){
  file1 <- "atac.merged.Peak_atac.merged.Cov__MERGED.tsv"
  file2 <- "atac.merged.Peak_wgbs.Cov__MERGED.tsv"
  data1 <- data.frame(read_delim(file1, "\t", escape_double = FALSE, trim_ws = TRUE), stringsAsFactors = F)[,c(1:4,9:11)]
  data2 <- data.frame(read_delim(file2, "\t", escape_double = FALSE, trim_ws = TRUE), stringsAsFactors = F)[,c(1:4,9:11)]
  data1$coverage <- data1$meth.counts + data1$un.meth.counts
  data2$coverage <- data2$meth.counts + data2$un.meth.counts
  data1$ID <- paste0(data1$peak.chr,":", data1$peak.start, ":", data1$peak.end)
  data2$ID <- paste0(data2$peak.chr,":", data2$peak.start, ":", data2$peak.end)
  
  outfile <- paste0(strsplit(file1, "__")[[1]][1], "___",
                    strsplit(file2, "__")[[1]][1], ".pdf")
  pdf(outfile, onefile = T)
  tmp1 <- unique(data1[,c(7,9)])
  tmp2 <- unique(data2[,c(7,9)])
  merged <- merge(tmp1, tmp2, by="ID")
  names(merged)[2:3] <- c(strsplit(file1, "__")[[1]][1], strsplit(file2, "__")[[1]][1])
  cor.value <- cor(x=merged[,2], y=merged[,3], method = "pearson")
  less.less <- sum (merged[,2]<50 & merged[,3]<50)
  less.more <- sum (merged[,2]<50 & merged[,3]>50)
  more.less <- sum (merged[,2]>50 & merged[,3]<50)
  more.more <- sum (merged[,2]>50 & merged[,3]>50)
  numbers <- c(less.less, less.more, more.more, more.less)
  print (head(merged))
  p <- ggplot(merged,aes(x=merged[,2], y=merged[,3])) +
    geom_point(color="grey", alpha=0.5, size=0.8) +
    theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) +
    xlab(paste (names(merged)[2], list(numbers))) +
    ylab(names(merged)[3]) +
    ggtitle(paste0("Average methylation in Peaks (coverage cutoff=0)", "  PCC=",round(cor.value,4))) +
    theme(plot.title = element_text(size = 12),
          axis.text.x=element_text(colour="black", size = 12),
          axis.text.y=element_text(colour="black", size = 12)) +
    #stat_smooth(method=lm, level=0.95) +
    stat_smooth(method=loess, level=0.95, fill = "yellow") +
    ylim(c(0,100)) +
    geom_vline(xintercept = 50, col="brown") +
    geom_hline(yintercept = 50, col="brown")
  print (p)
  
  head(merged)
  
  data1 <- data1[data1$coverage>5,]
  data2 <- data2[data2$coverage>5,]
  tmp1 <- unique(data1[,c(7,9)])
  tmp2 <- unique(data2[,c(7,9)])
  merged <- merge(tmp1, tmp2, by="ID")
  names(merged)[2:3] <- c(strsplit(file1, "__")[[1]][1], strsplit(file2, "__")[[1]][1])
  cor.value <- cor(x=merged[,2], y=merged[,3], method = "pearson")
  less.less <- sum (merged[,2]<50 & merged[,3]<50)
  less.more <- sum (merged[,2]<50 & merged[,3]>50)
  more.less <- sum (merged[,2]>50 & merged[,3]<50)
  more.more <- sum (merged[,2]>50 & merged[,3]>50)
  numbers <- c(less.less, less.more, more.more, more.less)
  
  write.table(merged, 
              paste0(names(merged)[2],"_", names(merged)[3],"_cov.cutoff.5",".txt"), 
              col.names = T, row.names = F, quote = F, sep = "\t")
  print (head(merged))
  
  p <- ggplot(merged,aes(x=merged[,2], y=merged[,3])) +
    geom_point(color="grey", alpha=0.5, size=0.8) +
    theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) +
    xlab(paste (names(merged)[2], list(numbers))) +
    ylab(names(merged)[3]) +
    ggtitle(paste0("Average methylation in Peaks (coverage cutoff=5)", "  PCC=",round(cor.value,4))) +
    theme(plot.title = element_text(size = 12),
          axis.text.x=element_text(colour="black", size = 12),
          axis.text.y=element_text(colour="black", size = 12)) +
    #stat_smooth(method=lm, level=0.95) +
    stat_smooth(method=loess, level=0.95, fill = "yellow") +
    ylim(c(0,100)) +
    geom_vline(xintercept = 50, col="brown") +
    geom_hline(yintercept = 50, col="brown")
  print (p)
  
  
  data1 <- data1[data1$coverage>10,]
  data2 <- data2[data2$coverage>10,]
  tmp1 <- unique(data1[,c(7,9)])
  tmp2 <- unique(data2[,c(7,9)])
  merged <- merge(tmp1, tmp2, by="ID")
  names(merged)[2:3] <- c(strsplit(file1, "__")[[1]][1], strsplit(file2, "__")[[1]][1])
  cor.value <- cor(x=merged[,2], y=merged[,3], method = "pearson")
  less.less <- sum (merged[,2]<50 & merged[,3]<50)
  less.more <- sum (merged[,2]<50 & merged[,3]>50)
  more.less <- sum (merged[,2]>50 & merged[,3]<50)
  more.more <- sum (merged[,2]>50 & merged[,3]>50)
  
  write.table(merged, 
              paste0(names(merged)[2],"_", names(merged)[3],"_cov.cutoff.10",".txt"), 
              col.names = T, row.names = F, quote = F, sep = "\t")
  numbers <- c(less.less, less.more, more.more, more.less)
  p <- ggplot(merged,aes(x=merged[,2], y=merged[,3])) +
    geom_point(color="grey", alpha=0.5, size=0.8) +
    theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) +
    xlab(paste (names(merged)[2], list(numbers))) +
    ylab(names(merged)[3]) +
    ggtitle(paste0("Average methylation in Peaks (coverage cutoff=10)", "  PCC=",round(cor.value,4))) +
    theme(plot.title = element_text(size = 12),
          axis.text.x=element_text(colour="black", size = 12),
          axis.text.y=element_text(colour="black", size = 12)) +
    #stat_smooth(method=lm, level=0.95) +
    stat_smooth(method=loess, level=0.95, fill = "yellow") +
    ylim(c(0,100)) +
    geom_vline(xintercept = 50, col="brown") +
    geom_hline(yintercept = 50, col="brown")
  print (p)
  dev.off()
}


make_average_meth_plot(file1="atac.merged.Peak_atac.merged.Cov__MERGED.tsv",
                       file2="atac.merged.Peak_wgbs.Cov__MERGED.tsv")

make_average_meth_plot(file1="ctcf.merged.Peak_ctcf.merged.Cov__MERGED.tsv",
                       file2="ctcf.merged.Peak_wgbs.Cov__MERGED.tsv")

make_average_meth_plot(file1="mESctcf.merged.Peak_mESctcf.merged.Cov__MERGED.tsv",
                       file2="mESctcf.merged.Peak_wgbs.Cov__MERGED.tsv")

make_average_meth_plot(file1="m-atac-merged.Peak_atac-3.merged.Cov__MERGED.tsv",
                       file2="m-atac-merged.Peak_atac-4.merged.Cov__MERGED.tsv")

make_average_meth_plot(file1="ctcf-merged.Peak_ctcf-1.merged.Cov__MERGED.tsv",
                       file2="ctcf-merged.Peak_ctcf-2.merged.Cov__MERGED.tsv")



# get coordinates of 4 groups for anntation
head(merged)
library(splitstackshape)

df <- (cSplit(merged, "ID", sep=":"))
df <- data.frame(df, stringsAsFactors = F)
df <- df[,c(3:5,1:2)]
names(df) <- c("chr", "start", "end", "atac", "wgs")
df$chr <- as.character(df$chr)

m.less.than.50.w.less.than.50 <- df[(df$atac<50 & df$wgs<50) , c(1:3)]
m.more.than.50.w.more.than.50 <- df[(df$atac>50 & df$wgs>50) , c(1:3)]
m.less.than.50.w.more.than.50 <- df[(df$atac<50 & df$wgs>50) , c(1:3)]
m.more.than.50.w.less.than.50 <- df[(df$atac>50 & df$wgs<50) , c(1:3)]

write.table(m.less.than.50.w.less.than.50, "m.less.than.50.w.less.than.50.bed", col.names = F, row.names = F, sep="\t", quote = F)
write.table(m.more.than.50.w.more.than.50,"m.more.than.50.w.more.than.50.bed", col.names = F, row.names = F, sep="\t", quote = F)
write.table(m.less.than.50.w.more.than.50, "m.less.than.50.w.more.than.50.bed", col.names = F, row.names = F, sep="\t", quote = F)
write.table(m.more.than.50.w.less.than.50, "m.more.than.50.w.less.than.50.bed", col.names = F, row.names = F, sep="\t", quote = F)



write.table(m.less.than.50.w.less.than.50, "c.less.than.50.w.less.than.50.bed", col.names = F, row.names = F, sep="\t", quote = F)
write.table(m.more.than.50.w.more.than.50,"c.more.than.50.w.more.than.50.bed", col.names = F, row.names = F, sep="\t", quote = F)
write.table(m.less.than.50.w.more.than.50, "c.less.than.50.w.more.than.50.bed", col.names = F, row.names = F, sep="\t", quote = F)
write.table(m.more.than.50.w.less.than.50, "c.more.than.50.w.less.than.50.bed", col.names = F, row.names = F, sep="\t", quote = F)





