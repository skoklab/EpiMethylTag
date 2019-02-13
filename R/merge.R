##############################################################################################
#  Authors : Gunjan Sethia, Skok lab, Dept. Pathology, NYU Langone Health  #
##############################################################################################

library(optparse)
args = commandArgs(trailingOnly=TRUE)
sample=args[1]
peak.file=args[2]
meth.file=args[3]

print (sample)
print (peak.file)
print (meth.file)

# peaks <- read.table("mESC_ATAC-BS_peakdeck.SORTED.MERGED.bed")
# meth <- read.table("mESC_ATAC-BS.bismark.cov")

peaks <- read.table(peak.file)
meth <- read.table(meth.file)

colnames(meth) <- c("m.chr","m.start","m.end","m.percent","m.counts","un.m.counts")
colnames(peaks) <- c("p.chr","p.start","p.end")
str(meth)
str(peaks)
peaks$p.chr <- as.character(peaks$p.chr)

str(meth)
meth$m.chr <- as.character(meth$m.chr)
meth <- meth[meth$m.chr!="chrM",]

unique(meth$m.chr)
unique(peaks$p.chr)

FINAL <- data.frame()

for (each in unique(meth$m.chr)){
  tmp.meth <- meth[meth$m.chr==each,]
  tmp.peaks <- peaks[peaks$p.chr==each,]
  for (i in 1:nrow(tmp.peaks)){
    chr <- tmp.peaks$p.chr[i]
    print (paste(chr,":",which(unique(meth$m.chr)==chr),"of", length(unique(meth$m.chr))))
    print (dim(tmp.peaks))
    print (i)
    tmp <- tmp.meth[tmp.meth$m.start>=tmp.peaks$p.start[i] & tmp.meth$m.end <= tmp.peaks$p.end[i] ,]
    if (nrow(tmp)!=0){
      ID <- NULL
      all <- NULL
      all.df <- NULL
      ID <- paste0("peak.",i)
      all <- cbind(ID,
                   tmp.peaks$p.chr[i],
                   tmp.peaks$p.start[i],
                   tmp.peaks$p.end[i],
                   tmp)
      all.df <- as.data.frame(all,stringAsFactors=F)
      colnames(all.df) <- c("peakID","peak.chr","peak.start","peak.end",
                            "meth.chr","meth.start","meth.end","meth.percent",
                            "meth.counts","un.meth.counts")
      all.df$avg.meth <- round(mean(all.df$meth.percent),3)
      FINAL <- data.frame(rbind(FINAL,all.df))

    }
  }
}

filename=paste0(sample,"__MERGED.tsv")

write.table(FINAL,filename,col.names = T,row.names = F,quote = F,sep="\t")

