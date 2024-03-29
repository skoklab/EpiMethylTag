##############################################################################################
#  Authors : Priscillia Lhoumaud & Sana Badri, Skok lab, Dept. Pathology, NYU Langone Health #
##############################################################################################

setwd("C:/Users/lhoump01/Google Drive/method paper/Samples_December2017/Motifs/December2017_CTCF_MA0139/Violin/")
library(ggplot2)

############# Violin plot after removing CpGs with less than 5 reads
CTCF=as.matrix(read.table("C2_C12_CTCF_WGBS_filt5_forviolin.txt", sep="\t", na.strings = "NA",header=T))
CTCF <- as.data.frame(CTCF)
CTCF$Meth_Percent <- as.numeric(as.character(CTCF$Meth_Percent))

colnames(CTCF)

tab_plot <- CTCF[,1:2]
colnames(tab_plot)
dim(tab_plot)
levels(tab_plot$GROUP)

# medians
CTCFC2_median <- round(median(tab_plot$Meth_Percent[tab_plot$GROUP == "C2_CTCF"]),digits = 5)
CTCFC12_median <- round(median(tab_plot$Meth_Percent[tab_plot$GROUP == "C12_CTCF"]),digits = 5)
WGBSC2_median <- round(median(tab_plot$Meth_Percent[tab_plot$GROUP == "C2_WGBS"]),digits = 5)
WGBSC12_median <- round(median(tab_plot$Meth_Percent[tab_plot$GROUP == "C12_WGBS"]),digits = 5)
CTCFC2_mean <- round(mean(tab_plot$Meth_Percent[tab_plot$GROUP == "C2_CTCF"]),digits = 3)
CTCFC12_mean <- round(mean(tab_plot$Meth_Percent[tab_plot$GROUP == "C12_CTCF"]),digits = 3)
WGBSC2_mean <- round(mean(tab_plot$Meth_Percent[tab_plot$GROUP == "C2_WGBS"]),digits = 3)
WGBSC12_mean <- round(mean(tab_plot$Meth_Percent[tab_plot$GROUP == "C12_WGBS"]),digits = 3)

tab_plot$GROUP <- factor(tab_plot$GROUP,levels = c('C2_CTCF','C12_CTCF','C2_WGBS','C12_WGBS'),ordered = TRUE)
levels(tab_plot$GROUP)

# violin plots
#ggplot(tab_plot, aes(x=GROUP, y=Meth_Percent, fill=GROUP))+xlab(label = "CTCF_MA0139.1") +geom_violin(trim=F)+geom_boxplot(width=0.1,outlier.size = 0.1, fill="white")+ylab(label = "log2(Methylation Percentage per CpG")+scale_fill_manual(values=c("red", "blue","darkred","darkblue"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggplot(tab_plot, aes(x=GROUP, y=(Meth_Percent), fill=GROUP))+xlab(label = "CTCF_MA0139.1") +geom_violin(trim=F)+ylab(label = "log2(Methylation Percentage per CpG)")+scale_fill_manual(values=c("red", "blue","darkred","darkblue"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
#ggplot(tab_plot, aes(x=GROUP, y=log2(Meth_Percent), fill=GROUP))+xlab(label = "CTCF_MA0139.1") +geom_violin(trim=F)+geom_boxplot(width=0.1,outlier.size = 0.1, fill="white")+ylab(label = "log2(Methylation Percentage per CpG")+scale_fill_manual(values=c("red", "blue","darkred","darkblue"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
#ggplot(tab_plot, aes(x=GROUP, y=sqrt(Meth_Percent), fill=GROUP))+xlab(label = "CTCF_MA0139.1") +geom_violin(trim=F)+ylab(label = "sqrt(Methylation Percentage per CpG")+scale_fill_manual(values=c("red", "blue","darkred","darkblue"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))



#ggplot(tab_plot2, aes(x=group, y=values, group=group, fill=group))+xlab(label = "NSD2") + geom_violin(trim = F)+ geom_boxplot(width=0.1,outlier.size = 0.1, fill="white")+ylab(label = "log2FC (RNA)")+scale_fill_manual(values=c("red", "blue", "grey"))+geom_text(data = p_meds, aes(x = group, y = med, label = med), size = 2, vjust = -0.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+geom_signif(comparisons = list(c("A to B", "B to A"),c("A to B", "STABLE"), c("B to A", "STABLE")), map_signif_level=TRUE, test = 'wilcox.test')

# t-test C2 vs C12 CTCF M-ChIP vs WGBS (after filtering for more than 5 reads)
CTCF_ttest=as.matrix(read.table("C2_C12_CTCF_WGBS_filt5_stats.txt", sep="\t", na.strings = "NA",header=T))
CTCF_test2=log2(CTCF_ttest)
CTCF_test3=asin(CTCF_ttest/100)
hist(CTCF_test3)
colnames(CTCF_ttest)
colnames(CTCF_test3)
colnames(CTCF_test2)
View(CTCF_ttest)
boxplot(CTCF_ttest[,1:4])
boxplot(CTCF_test2[,1:4])
boxplot(CTCF_test3[,1:4])
nrow(CTCF_ttest[,1])

t.test(CTCF_ttest[,1],CTCF_ttest[,2])$p.value
t.test(CTCF_ttest[,1],CTCF_ttest[,3])$p.value
t.test(CTCF_ttest[,2],CTCF_ttest[,4])$p.value
t.test(CTCF_ttest[,3],CTCF_ttest[,4])$p.value
wilcox.test(CTCF_ttest[,1],CTCF_ttest[,2])$p.value
wilcox.test(CTCF_ttest[,1],CTCF_ttest[,3])$p.value
wilcox.test(CTCF_ttest[,2],CTCF_ttest[,4])$p.value
wilcox.test(CTCF_ttest[,3],CTCF_ttest[,4])$p.value
t.test(CTCF_ttest[,1],CTCF_ttest[,2], paired=TRUE)$p.value
t.test(CTCF_ttest[,3],CTCF_ttest[,4], paired=TRUE)$p.value
wilcox.test(CTCF_ttest[,1],CTCF_ttest[,2],paired=TRUE)$p.value
wilcox.test(CTCF_ttest[,3],CTCF_ttest[,4],paired=TRUE)$p.value

wilcox.test(CTCF_test2[,1],CTCF_test2[,3])$p.value
wilcox.test(CTCF_test2[,1],CTCF_test2[,2])$p.value
wilcox.test(CTCF_test2[,2],CTCF_test2[,4])$p.value
wilcox.test(CTCF_test2[,3],CTCF_test2[,4])$p.value



t.test(CTCF_test3[,1],CTCF_test3[,3])$p.value
t.test(CTCF_test3[,1],CTCF_test3[,2],paired=TRUE)$p.value
t.test(CTCF_test3[,2],CTCF_test3[,4])$p.value
t.test(CTCF_test3[,3],CTCF_test3[,4],paired=TRUE)$p.value

wilcox.test(CTCF_test3[,1],CTCF_test3[,3])$p.value
wilcox.test(CTCF_test3[,1],CTCF_test3[,2],paired=TRUE)$p.value
wilcox.test(CTCF_test3[,2],CTCF_test3[,4])$p.value
wilcox.test(CTCF_test3[,3],CTCF_test3[,4],paired=TRUE)$p.value

library(ggplot2)
ggplot(as.data.frame(CTCF_test3),aes(C2_CTCF,C2_WGBS))+
  geom_point()
ggplot(as.data.frame(CTCF_test3),aes(C12_CTCF,C12_WGBS))+
  geom_point()

hist(CTCF_test2[,1])
hist(CTCF_test2[,2])
hist(CTCF_test2[,3])
hist(CTCF_test2[,4])
hist(CTCF_ttest[,1])

# Model of data regression to test difference between CTCF M-ChIP vs WGBS (after filtering for more than 5 reads)
CTCF_ttest=as.matrix(read.table("C2_C12_CTCF_WGBS_filt5_stats.txt", sep="\t", na.strings = "NA",header=T))
df = as.data.frame(CTCF_ttest)
df = as.data.frame(CTCF_ttest/100)
df=(tidyr::gather_(df, "key", "value", colnames(df)))
df = tidyr::separate_(df, "key", c("position", "assay"), sep = "_")
df$value = df$value + 0.001 * (df$value == 0)
df$value = df$value - 0.001 * (df$value == 1)
df$logit = log(df$value) - log(1 - df$value)
summary(lm(logit ~ position + assay, data = df))
inv.logit <- function(x) 1 / (1 + exp(-x))
fit=betareg(value ~ position * assay, data = df)
summary(fit)

