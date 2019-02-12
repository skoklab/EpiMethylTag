library(readr)
library(grDevices)
library(ggplot2)

work.dir <- "~/Google Drive/NYUMC/Skok lab/BiSulphiteGenome_analysis/Jan_2018/"
setwd(work.dir)

meth.vs.coverage <- function(file1, coverage.cutoff){
  
  #file1 <- "atac.merged.Peak_atac.merged.Cov__MERGED.tsv"
  #coverage.cutoff=5
  sample1 <- strsplit(file1,"__",fixed=T)[[1]][1]
  
  title1 <- paste0(sample1,"_","coverage_cutoff_0-groups3")
  title2 <- paste0(sample1,"_","coverage_cutoff_",coverage.cutoff,"-groups3")
  
  data1 <- data.frame(read_delim(file1, "\t", escape_double = FALSE, trim_ws = TRUE, col_names = TRUE),stringsAsFactors = F)[,5:10]
  colnames(data1) <- c("meth.chr", "meth.start", "meth.end", "meth.percent", "meth.counts", "un.meth.counts")
  data1 <- data1[complete.cases(data1),]
  
  data1$coverage <- (data1$meth.counts + data1$un.meth.counts)
  data1$group <- "group"

  data1$group[(data1$meth.percent>=0 & data1$meth.percent<=20.5) ] <- "0-20"
  data1$group[(data1$meth.percent>=20.5 & data1$meth.percent<=80.5) ] <- "21-80"
  data1$group[(data1$meth.percent>=80.5 & data1$meth.percent<=100)] <- "81-100"
  print (table(data1$group))
  print (data1[data1$group=="group",])
  
  data2 <- data1[data1$coverage > coverage.cutoff,]
  print (table(data2$group))
  
  print ("plot1")
  
  pdf(paste0("meth_vs_cov_plots/", title1, ".pdf"))
  p1 <- ggplot(data1, aes(group, log10(coverage))) + geom_violin(aes(fill=group)) +
    xlab (paste("Meth Groups","Medians are:", 
                median(data1$coverage[data1$group=="0-20"]),
                median(data1$coverage[data1$group=="21-80"]),
                median(data1$coverage[data1$group=="81-100"]) )) +
    stat_summary(fun.y=median, geom="point", size=2, color="black") +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggtitle(title1)
  print (p1)
  dev.off()
  
  print("plot2")
  print (dim(data2))
  print(names(data2))
  pdf(paste0("meth_vs_cov_plots/", title2, ".pdf"))
  p2 <- ggplot(data2, aes(group, log10(coverage))) + geom_violin(aes(fill=group)) +
    xlab (paste("Meth Groups","Medians are:", 
                median(data1$coverage[data1$group=="0-20"]),
                median(data1$coverage[data1$group=="21-80"]),
                median(data1$coverage[data1$group=="81-100"]) )) +
    stat_summary(fun.y=median, geom="point", size=2, color="black") +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggtitle (title2)
  print (p2)
  dev.off()
  
  
  wilcoxtest <- function (vector1, vector2) {
    x <- wilcox.test(vector1, vector2, paired=F, exact = F)
    return (x$p.value)
  }  
  
  
  pvals <- data.frame(rbind(c("group1", "group2", wilcoxtest(log10(data1$coverage[data1$group=="0-20"]) , log10(data1$coverage[data1$group=="21-80"]))),
                            c("group1", "group3", wilcoxtest(log10(data1$coverage[data1$group=="0-20"]) , log10(data1$coverage[data1$group=="81-100"]))),
                            c("group2", "group3", wilcoxtest(log10(data1$coverage[data1$group=="21-80"]) , log10(data1$coverage[data1$group=="81-100"]))) ))
  names(pvals) <- c("Group-A", "Group-B", "wilcox.Pvalue")
  write.table(pvals, paste0("meth_vs_cov_plots/",title1,"_pvals.txt"), row.names = F, col.names = T, quote = F, sep= "\t")
  
  pvals <- data.frame(rbind(c("group1", "group2", wilcoxtest(log10(data2$coverage[data2$group=="0-20"]) , log10(data2$coverage[data2$group=="21-80"]))),
                            c("group1", "group3", wilcoxtest(log10(data2$coverage[data2$group=="0-20"]) , log10(data2$coverage[data2$group=="81-100"]))),
                            c("group2", "group3", wilcoxtest(log10(data2$coverage[data2$group=="21-80"]) , log10(data2$coverage[data2$group=="81-100"]))) ))
  names(pvals) <- c("Group-A", "Group-B", "wilcox.Pvalue")
  write.table(pvals, paste0("meth_vs_cov_plots/",title2,"_pvals.txt"), row.names = F, col.names = T, quote = F, sep= "\t")
  

  
}


meth.vs.coverage("atac.merged.Peak_wgbs.Cov__MERGED.tsv", 5)
meth.vs.coverage("atac.merged.Peak_atac.merged.Cov__MERGED.tsv", 5)
meth.vs.coverage("ctcf.merged.Peak_ctcf.merged.Cov__MERGED.tsv", 5)
meth.vs.coverage("ctcf.merged.Peak_wgbs.Cov__MERGED.tsv", 5)
meth.vs.coverage("mESctcf.merged.Peak_mESctcf.merged.Cov__MERGED.tsv", 5)
meth.vs.coverage("mESctcf.merged.Peak_wgbs.Cov__MERGED.tsv", 5)

# coordinate data for each violin 3 groups
head(data2)
data2$ID <- paste0(data2$meth.chr, ":", data2$meth.start, ":", data2$meth.end)
df <- data2


write.table(data2[data2$group=="0-20" & data2$coverage<50,c(1:3)], 
            "atac.peak.atac.cov_0-20_meth_cov_less_50.bed", sep="\t", row.names = F, col.names = F, quote = F)
write.table(data2[data2$group=="0-20" & data2$coverage>50,c(1:3)], 
            "atac.peak.atac.cov_0-20_meth_cov_more_50.bed", sep="\t", row.names = F, col.names = F, quote = F)
write.table(data2[data2$group=="21-80",c(1:3)], 
            "atac.peak.atac.cov_21-80_meth.bed", sep="\t", row.names = F, col.names = F, quote = F)
write.table(data2[data2$group=="81-100",c(1:3)], 
            "atac.peak.atac.cov_81-100_meth.bed", sep="\t", row.names = F, col.names = F, quote = F)





