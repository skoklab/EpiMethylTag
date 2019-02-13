##############################################################################################
#  Authors : Priscillia Lhoumaud & Sana Badri, Skok lab, Dept. Pathology, NYU Langone Health #
##############################################################################################


setwd("C:/Users/lhoump01/Google Drive/method paper/Samples_December2017/Gunjan/scatter-plots/scatter-plots_ATAC_and_CTCF-PDFs/methylation/")

############## ATAC versus WGBS #####################
x<-read.delim("atac.merged.Peak_atac.merged.Cov_atac.merged.Peak_wgbs.Cov_cov.cutoff.5.txt")
# scatter plot of x and y variables
# color by groups
library(ggplot2)

blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank())


source("http://bioconductor.org/biocLite.R")
biocLite("hexbin")
library(hexbin)
library(grid)


cor.test(x$atac.merged.Peak_atac.merged.Cov, x$atac.merged.Peak_wgbs.Cov, method = "pearson", conf.level = 0.95)
cor.value = cor.test(x$atac.merged.Peak_atac.merged.Cov, x$atac.merged.Peak_wgbs.Cov, method = "pearson", conf.level = 0.95)
my_breaks = c(2, 10, 50, 250, 1250)
# Number of bins in each direction?
ggplot(x,aes(atac.merged.Peak_atac.merged.Cov, atac.merged.Peak_wgbs.Cov) ) +
  stat_binhex(bins=100)+
  theme_bw()+
  xlim(0,100)+
  ylim(0,100)+
  scale_fill_gradient(name = "count", trans = "log",low="black", high="red",
                      breaks = my_breaks, labels = my_breaks)+
  geom_smooth(method="lm")+
  #geom_smooth()+
  #stat_smooth(method=loess, level=0.95, fill = "yellow") +
  ggtitle(paste0("Average methylation in Peaks (coverage cutoff=5)" , "PCC=",round(cor.value$estimate,4)))



############## CTCF versus WGBS #####################
x<-read.delim("ctcf.merged.Peak_ctcf.merged.Cov_ctcf.merged.Peak_wgbs.Cov_cov.cutoff.5.txt")
# scatter plot of x and y variables
# color by groups
#library(ggplot2)

blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

#library("gridExtra")

#source("http://bioconductor.org/biocLite.R")
#biocLite("hexbin")
#library(hexbin)
#library(grid)


cor.test(x$ctcf.merged.Peak_ctcf.merged.Cov, x$ctcf.merged.Peak_wgbs.Cov, method = "pearson", conf.level = 0.95)
cor.value = cor.test(x$ctcf.merged.Peak_ctcf.merged.Cov, x$ctcf.merged.Peak_wgbs.Cov, method = "pearson", conf.level = 0.95)
my_breaks = c(2, 10, 50, 250, 1250)
# Number of bins in each direction?
ggplot(x,aes(ctcf.merged.Peak_ctcf.merged.Cov, ctcf.merged.Peak_wgbs.Cov) ) +
  stat_binhex(bins=100)+
  theme_bw()+
  scale_fill_gradient(name = "count", trans = "log",low="black", high="red",
                      breaks = my_breaks, labels = my_breaks)+
  geom_smooth(method="lm")+
  #stat_smooth(method=loess, level=0.95, fill = "yellow") +
  ggtitle(paste0("Average methylation in Peaks (coverage cutoff=5)" , "PCC=",round(cor.value$estimate,4)))




############## CTCF Schubeler versus WGBS #####################
x<-read.delim("mESctcf.merged.Peak_mESctcf.merged.Cov_mESctcf.merged.Peak_wgbs.Cov_cov.cutoff.5.txt")
# scatter plot of x and y variables
# color by groups
library(ggplot2)
colnames(x)
        

blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

library("gridExtra")

#source("http://bioconductor.org/biocLite.R")
#biocLite("hexbin")
library(hexbin)
library(grid)

cor.test(x$mESctcf.merged.Peak_mESctcf.merged.Cov, x$mESctcf.merged.Peak_wgbs.Cov, method = "pearson", conf.level = 0.95)
cor.value = cor.test(x$mESctcf.merged.Peak_mESctcf.merged.Cov, x$mESctcf.merged.Peak_wgbs.Cov, method = "pearson", conf.level = 0.95)
my_breaks = c(2, 10, 50, 250, 1250)
# Number of bins in each direction?
ggplot(x,aes(mESctcf.merged.Peak_mESctcf.merged.Cov, mESctcf.merged.Peak_wgbs.Cov) ) +
  stat_binhex(bins=100)+
  theme_bw()+
  scale_fill_gradient(name = "count", trans = "log",low="black", high="red",
                      breaks = my_breaks, labels = my_breaks)+
                      geom_smooth(method="lm")+
#stat_smooth(method=loess, level=0.95, fill = "yellow") +
ggtitle(paste0("Average methylation in Peaks (coverage cutoff=5)" , "PCC=",round(cor.value$estimate,4)))


############## M-ATAC replicate 1 versus 2 #####################
x<-read.delim("m-atac-merged.Peak_atac-3.merged.Cov_m-atac-merged.Peak_atac-4.merged.Cov_cov.cutoff.5.txt")
# scatter plot of x and y variables
# color by groups
#library(ggplot2)
colnames(x)


blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

#library("gridExtra")

#source("http://bioconductor.org/biocLite.R")
#biocLite("hexbin")
#library(hexbin)
#library(grid)

cor.test(x$m.atac.merged.Peak_atac.3.merged.Cov, x$m.atac.merged.Peak_atac.4.merged.Cov, method = "pearson", conf.level = 0.95)
cor.value = cor.test(x$m.atac.merged.Peak_atac.3.merged.Cov, x$m.atac.merged.Peak_atac.4.merged.Cov, method = "pearson", conf.level = 0.95)
my_breaks = c(2, 10, 50, 250, 1250)
# Number of bins in each direction?
ggplot(x,aes(m.atac.merged.Peak_atac.3.merged.Cov, m.atac.merged.Peak_atac.4.merged.Cov) ) +
  stat_binhex(bins=100)+
  theme_bw()+
  scale_fill_gradient(name = "count", trans = "log",low="black", high="red",
                      breaks = my_breaks, labels = my_breaks)+
  geom_smooth(method="lm")+
  #stat_smooth(method=loess, level=0.95, fill = "yellow") +
  ggtitle(paste0("Average methylation in Peaks (coverage cutoff=5)" , "PCC=",round(cor.value$estimate,4)))


############## CTCF M-ChIP replicate 1 versus 2 #####################
x<-read.delim("ctcf-merged.Peak_ctcf-1.merged.Cov_ctcf-merged.Peak_ctcf-2.merged.Cov_cov.cutoff.5.txt")
# scatter plot of x and y variables
# color by groups
#library(ggplot2)
colnames(x)


blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

#library("gridExtra")

#source("http://bioconductor.org/biocLite.R")
#biocLite("hexbin")
#library(hexbin)
#library(grid)

cor.test(x$ctcf.merged.Peak_ctcf.1.merged.Cov, x$ctcf.merged.Peak_ctcf.2.merged.Cov, method = "pearson", conf.level = 0.95)
cor.value = cor.test(x$ctcf.merged.Peak_ctcf.1.merged.Cov, x$ctcf.merged.Peak_ctcf.2.merged.Cov, method = "pearson", conf.level = 0.95)
my_breaks = c(2, 10, 50, 250, 1250)
# Number of bins in each direction?
ggplot(x,aes(ctcf.merged.Peak_ctcf.1.merged.Cov, ctcf.merged.Peak_ctcf.2.merged.Cov) ) +
  stat_binhex(bins=100)+
  theme_bw()+
  scale_fill_gradient(name = "count", trans = "log",low="black", high="red",
                      breaks = my_breaks, labels = my_breaks)+
  geom_smooth(method="lm")+
  #stat_smooth(method=loess, level=0.95, fill = "yellow") +
  ggtitle(paste0("Average methylation in Peaks (coverage cutoff=5)" , "PCC=",round(cor.value$estimate,4)))




