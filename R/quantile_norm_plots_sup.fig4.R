##############################################################################################
#  Authors : Gunjan Sethia, Skok lab, Dept. Pathology, NYU Langone Health  #
##############################################################################################

library(readr)
library(ggplot2)
library(dplyr)

setwd("~/Google Drive/NYUMC/Skok lab/BiSulphiteGenome_analysis/Jan_2018/quantile-TSV-MERGED/")

file="atac.merged.Peak_atac.merged.Cov__MERGED.tsv"
file="ctcf.merged.Peak_ctcf.merged.Cov__MERGED.tsv"
file="mESctcf.merged.Peak_mESctcf.merged.Cov__MERGED.tsv"

data <- data.frame(read_delim(paste0("~/Google Drive/NYUMC/Skok lab/BiSulphiteGenome_analysis/Jan_2018/quantile-TSV-MERGED/",file), 
                                                      "\t", escape_double = FALSE, trim_ws = TRUE), stringsAsFactors = F)
head(data)
sample <- strsplit(file,"__")[[1]][1]
data$round.Q <- round(data$quantile)
data.MAIN <- data
data <- data.MAIN[data.MAIN$coverage>5,c(8,14)]
z <- data.frame(table(data$round.Q))
names(z) <- c("quantile", "counts")
#data$round.Q <- factor(data$round.Q)
df <- data.frame()
df <- aggregate(data[,1], list(data$round.Q), mean)
names(df) <- c("round.Q","meanMETH")
library(Rmisc)
se.summary <- summarySE(data = data, measurevar="meth.percent", groupvars = "round.Q")
tmp <- merge(df, se.summary[,c(1,5)], by="round.Q")
df.final <- unique(tmp)

pdf(paste0(sample,"_coverage_cutoff_5.pdf"))
ggplot(df.final, aes(x=round.Q,y=meanMETH))  + 
  geom_point(size=2,alpha=0.6, col="blue")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(sample) +
  xlab("peak quantile") +
  ylab("average (% methylation)") +
  geom_errorbar(aes(ymin=meanMETH-se, ymax=meanMETH+se), width=.4, col="brown", 
                position = "dodge" )
dev.off()
  
