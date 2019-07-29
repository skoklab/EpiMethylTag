library(ggplot2)
library(readr)
library(scales)
library(RColorBrewer)

setwd("~/Google Drive/NYUMC/Skok lab/BiSulphiteGenome_analysis/only_scatter")

scatter_plot <- function(file, col1, col2){
  message (file)
  type <- strsplit(strsplit(file,"_")[[1]][3],".",fixed=T)[[1]][1]
  data <- data.frame(read_delim(file,"\t", escape_double = FALSE, trim_ws = TRUE),stringsAsFactors = F)
  print (dim(data))
  data <- data[ , c(col1,col2)]
  print (head(data))
  
  cor.value <- round(cor(data[,1],data[,2],method = "pearson") ,3)
  print (cor.value)
  
  out.file <- paste0(type,"--",names(data)[1], "--vs--", names(data)[2],".pdf")
  print (out.file)
  
  pdf(out.file)
  q <- ggplot(data, aes(x=data[,1], y=data[,2])) +
    geom_point(alpha=0.35,size=1.0, col="lightgrey") + 
    stat_smooth(method=lm, level=0.99) + 
    theme_bw() + 
    #annotate("text", x=10, y=80, label=paste0("PCC=",round(cor.value,3)),colour="red") +
    xlab(paste0(names(data)[1], "  PCC=", cor.value)) +
    ylab(names(data)[2]) +
    ggtitle(paste(names(data)[1],"vs", names(data)[2])) 
  #xlim(c(0,350)) +
  #ylim(c(0,400))
  
  #print (p)
  
  p <- ggplot(data, aes(x=log2(data[,1]), y=log2(data[,2]))) +
    #geom_point(size=0.2, alpha=0.2) + 
    stat_density2d(aes(fill = ..level.., alpha = ..level..), geom = 'polygon', bins=100) +
    scale_fill_gradient(low="blue", high = "red") +
    #scale_color_distiller(palette = "RdPu") +
    scale_alpha() +
    theme(legend.position = "none", axis.title = element_blank(), text = element_text(size = 12)) +
    theme_bw()+
    xlab(paste0(names(data)[1], "  PCC=", cor.value)) +
    ylab(names(data)[2]) +
    ggtitle(paste(names(data)[1],"vs", names(data)[2])) +
    stat_smooth(method=lm, level=0.99, geom="line", linetype="dashed",col="darkgrey") +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
  print (p)
  
  dev.off()
  
}

scatter_plot("Readscount_6samples_ATAC.txt", 2, 3)
scatter_plot("Readscount_6samples_ATAC.txt", 4, 5)
#scatter_plot("Readscount_6samples_ATAC.txt", 6, 7)
scatter_plot("Readscount_6samples_ATAC.txt", 2, 6)
scatter_plot("Readscount_6samples_ATAC.txt", 3, 7)
scatter_plot("Readscount_6samples_ATAC.txt", 2, 4)
scatter_plot("Readscount_6samples_ATAC.txt", 3, 5)
#scatter_plot("Readscount_6samples_ATAC.txt", 4, 6)
scatter_plot("Readscount_6samples_ATAC.txt", 5, 7)
scatter_plot("Readscounts_7samples_CTCF.txt", 2, 3)
scatter_plot("Readscounts_7samples_CTCF.txt", 4, 5)
scatter_plot("Readscounts_7samples_CTCF.txt", 6, 7)
scatter_plot("Readscounts_7samples_CTCF.txt", 7, 8)
scatter_plot("Readscounts_7samples_CTCF.txt", 2, 6)
scatter_plot("Readscounts_7samples_CTCF.txt", 3, 7)
scatter_plot("Readscounts_7samples_CTCF.txt", 2, 4)
scatter_plot("Readscounts_7samples_CTCF.txt", 3, 5)
scatter_plot("Readscounts_7samples_CTCF.txt", 4, 6)
scatter_plot("Readscounts_7samples_CTCF.txt", 5, 7)