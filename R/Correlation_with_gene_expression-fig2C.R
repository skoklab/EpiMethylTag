##############################################################################################
#  Authors : Gunjan Sethia, Skok lab, Dept. Pathology, NYU Langone Health  #
##############################################################################################

library(readr)
library(sqldf)
library(ggplot2)
library(splitstackshape)

work.dir <- "~/Google Drive/NYUMC/Skok lab/BiSulphiteGenome_analysis/Jan_2018/gene_expression"
setwd(work.dir)

counts <- data.frame(read_delim("RNAseq-mESC-merged.TPM.unstr.txt", "\t", escape_double = FALSE, col_names = TRUE,  trim_ws = TRUE), stringsAsFactors = F)
colnames(counts) <- c("Gene","counts")
counts <- data.frame(cSplit(counts, "Gene", ":"), stringsAsFactors = F)[,c(2,1)]
colnames(counts) <- c("Gene","counts")

######################################################################################################################################

g2 <- data.frame(read_delim("Group1.txt", "\t", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE), stringsAsFactors=F)
g1 <- data.frame(read_delim("Group2.txt", "\t", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE), stringsAsFactors=F)
g3 <- data.frame(read_delim("Group3.txt", "\t", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE), stringsAsFactors=F)
g4 <- data.frame(read_delim("Group4.txt", "\t", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE), stringsAsFactors=F)

names(g1) <-"Gene"
names(g2) <-"Gene"
names(g3) <-"Gene"
names(g4) <-"Gene"

g1 <- sqldf('SELECT DISTINCT * FROM g1')
g2 <- sqldf('SELECT DISTINCT * FROM g2')
g3 <- sqldf('SELECT DISTINCT * FROM g3')
g4 <- sqldf('SELECT DISTINCT * FROM g4')

c.g1 <- merge(g1, counts, by="Gene")
c.g2 <- merge(g2, counts, by="Gene")
c.g3 <- merge(g3, counts, by="Gene")
c.g4 <- merge(g4, counts, by="Gene")

c.g1$Groups <- "group1"
c.g2$Groups <- "group2"
c.g3$Groups <- "group3"
c.g4$Groups <- "group4"

final <- data.frame(rbind(c.g1, c.g2, c.g3, c.g4), stringsAsFactors = F)
n <- c(dim(c.g1)[1], dim(c.g2)[1], dim(c.g3)[1], dim(c.g4)[1])

ggplot(data=final ,aes(x=Groups,y=log10(counts), fill=Groups)) + 
  geom_boxplot() + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(paste("data in each group:", dim(c.g1)[1], dim(c.g2)[1], dim(c.g3)[1], dim(c.g4)[1]))

wilcoxtest <- function (vector1, vector2) {
  x <- wilcox.test(vector1, vector2, paired=F)
  return (x$p.value)
}  
  
pvals <- data.frame(rbind(c("group1", "group2", wilcoxtest(c.g1$counts , c.g2$counts)),
                          c("group1", "group3", wilcoxtest(c.g1$counts , c.g3$counts)),
                          c("group1", "group4", wilcoxtest(c.g1$counts , c.g4$counts)),
                          c("group2", "group3", wilcoxtest(c.g2$counts , c.g3$counts)),
                          c("group3", "group4", wilcoxtest(c.g2$counts , c.g4$counts)),
                          c("group3", "group4", wilcoxtest(c.g3$counts , c.g4$counts)) ))
names(pvals) <- c("GroupA", "GroupB", "wilcox.Pvalue")
write.table(pvals, "PL_merged_groups_pvals.txt", row.names = F, col.names = T, quote = F, sep= "\t")

######################################################################################################################################

g1  <- NULL
g2  <- NULL
g3  <- NULL
g4  <- NULL
final <- NULL

file2 <- "~/Google Drive/NYUMC/Skok lab/BiSulphiteGenome_analysis/Jan_2018/atac_groups/chipseeker_mar23/g1_atac-atac-meth_00to20-cov_00to50_TSS.bed"
file1 <- "~/Google Drive/NYUMC/Skok lab/BiSulphiteGenome_analysis/Jan_2018/atac_groups/chipseeker_mar23/g2_atac-atac-meth_00to20-cov_50andabove_TSS.bed"
file3 <- "~/Google Drive/NYUMC/Skok lab/BiSulphiteGenome_analysis/Jan_2018/atac_groups/chipseeker_mar23/g3_atac-atac-meth_20to80_TSS.bed"
file4 <- "~/Google Drive/NYUMC/Skok lab/BiSulphiteGenome_analysis/Jan_2018/atac_groups/chipseeker_mar23/g4_atac-atac-meth_80to100_TSS.bed"

g1 <- data.frame(read_delim(file1, "\t", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE), stringsAsFactors=F)[,c(15,16)]
names(g1) <- c("ensemble", "Gene")
g2 <- data.frame(read_delim(file2, "\t", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE), stringsAsFactors=F)[,c(15,16)]
names(g2) <- c("ensemble", "Gene")
g3 <- data.frame(read_delim(file3, "\t", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE), stringsAsFactors=F)[,c(15,16)]
names(g3) <- c("ensemble", "Gene")
g4 <- data.frame(read_delim(file4, "\t", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE), stringsAsFactors=F)[,c(15,16)]
names(g4) <- c("ensemble", "Gene")

g1 <- sqldf('SELECT DISTINCT * FROM g1')
g2 <- sqldf('SELECT DISTINCT * FROM g2')
g3 <- sqldf('SELECT DISTINCT * FROM g3')
g4 <- sqldf('SELECT DISTINCT * FROM g4')

write.table(g1$Gene[order(g1$Gene)], "Group1_Chipseeker.txt", row.names = F, col.names = F, quote = F, sep= "\t")
write.table(g2$Gene[order(g2$Gene)], "Group2_Chipseeker.txt", row.names = F, col.names = F, quote = F, sep= "\t")
write.table(g3$Gene[order(g3$Gene)], "Group3_Chipseeker.txt", row.names = F, col.names = F, quote = F, sep= "\t")
write.table(g4$Gene[order(g4$Gene)], "Group4_Chipseeker.txt", row.names = F, col.names = F, quote = F, sep= "\t")

c.g1 <- merge(g1, counts, by="Gene")
c.g2 <- merge(g2, counts, by="Gene")
c.g3 <- merge(g3, counts, by="Gene")
c.g4 <- merge(g4, counts, by="Gene")

c.g1$Groups <- "group1"
c.g2$Groups <- "group2"
c.g3$Groups <- "group3"
c.g4$Groups <- "group4"

final <- data.frame(rbind(c.g1, c.g2, c.g3, c.g4), stringsAsFactors = F)
n <- c(dim(c.g1)[1], dim(c.g2)[1], dim(c.g3)[1], dim(c.g4)[1])

ggplot(data=final ,aes(x=Groups,y=log10(counts), fill=Groups)) + 
  geom_boxplot() + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  xlab(paste("data in each group:", dim(c.g1)[1], dim(c.g2)[1], dim(c.g3)[1], dim(c.g4)[1]))

  wilcoxtest <- function (vector1, vector2) {
    x <- wilcox.test(vector1, vector2, paired=F)
    return (x$p.value)
  }  
  
  pvals <- data.frame(rbind(c("group1", "group2", wilcoxtest(c.g1$counts , c.g2$counts)),
                            c("group1", "group3", wilcoxtest(c.g1$counts , c.g3$counts)),
                            c("group1", "group4", wilcoxtest(c.g1$counts , c.g4$counts)),
                            c("group2", "group3", wilcoxtest(c.g2$counts , c.g3$counts)),
                            c("group3", "group4", wilcoxtest(c.g2$counts , c.g4$counts)),
                            c("group3", "group4", wilcoxtest(c.g3$counts , c.g4$counts)) ))
  names(pvals) <- c("GroupA", "GroupB", "wilcox.Pvalue")
  write.table(pvals, "Chipseeker_merged_groups_pvals.txt", row.names = F, col.names = T, quote = F, sep= "\t")
  
