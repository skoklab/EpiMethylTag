library(methylKit)
library(readr)
library(ggplot2)
library(reshape2)

setwd("/Volumes/DS/Franco/PriscilliaProject/RawData/")
wd = getwd()

# 1- Load ATAC data

atac <- read_delim("/Volumes/DS/Franco/PriscilliaProject/RawData/atac.merged.Peak_atac.merged.Cov__MERGED.tsv", 
                                                      "\t", escape_double = FALSE, trim_ws = TRUE)
atac = atac[atac$meth.counts+atac$un.meth.counts >= 5,]
quantile(atac$meth.counts+atac$un.meth.counts)

atac <- GRanges(seqnames = atac$peak.chr,
                ranges = IRanges(start = atac$meth.start, end = atac$meth.end),
                strand = "*",
                meth = atac$meth.percent)

wgbs <- read_delim("/Volumes/DS/Franco/PriscilliaProject/RawData/atac.merged.Peak_wgbs.Cov__MERGED.tsv", 
                   "\t", escape_double = FALSE, trim_ws = TRUE)

wgbs = wgbs[wgbs$meth.counts+wgbs$un.meth.counts >= 5,]
quantile(wgbs$meth.counts+wgbs$un.meth.counts)

wgbs <- GRanges(seqnames = wgbs$peak.chr,
                ranges = IRanges(start = wgbs$meth.start, end = wgbs$meth.end),
                strand = "*",
                meth = wgbs$meth.percent)

data = list(atac, wgbs)


# Load motifs #
TFs = dir("/Volumes/DS/Franco/PriscilliaProject/Motifs/MotifSites/", pattern = ".bed")
MotifOcurrence = list()
for(i in TFs){
  MotifOcurrence[[i]] = as.data.frame(read_delim(paste0("/Volumes/DS/Franco/PriscilliaProject/Motifs/MotifSites/",i),
                                                 "\t", escape_double = FALSE, col_names = FALSE,
                                                 trim_ws = TRUE))
  colnames(MotifOcurrence[[i]]) = c("Chr","Start","End","TF","Score","Strand")
  MotifOcurrence[[i]] = MotifOcurrence[[i]][complete.cases(MotifOcurrence[[i]]),]
  
  MotifOcurrence[[i]] = GRanges(seqnames = MotifOcurrence[[i]]$Chr,
                                ranges = IRanges(start = MotifOcurrence[[i]]$Start, end = MotifOcurrence[[i]]$End),
                                strand = Rle(values = MotifOcurrence[[i]]$Strand),
                                tf = MotifOcurrence[[i]]$TF,
                                score = MotifOcurrence[[i]]$Score)
}

names(MotifOcurrence) = unlist(lapply(strsplit(names(MotifOcurrence), split = "_"), function(x) paste0(x[1],"_",x[2])))
barplot(log10(unlist(lapply(MotifOcurrence, function(x) length(x@seqnames)))), las=2)

# Intersect meth with motif sites
MotifMeth = lapply(data, function(x) lapply(MotifOcurrence, function(y) {
  subsetByOverlaps(x,y)
}))
names(MotifMeth) = c("ATAC","WGBS")

MotifMeth$ATAC$Klf4_denovo$meth

m = melt(lapply(MotifMeth$ATAC, function(x) x$meth))

ggplot(m, aes(x = L1, y =value, fill = L1)) + geom_violin() +
 theme_classic() + 
  scale_fill_manual(values = c("red","red","blue","blue","darkgreen","darkgreen"))+
  labs(x = "Motif", y = "Methylation (%)", title = "ATAC-Bseq")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

mean(m$value[m$L1 == "Klf4_denovo"]); mean(m$value[m$L1 == "Klf4_known"])
wilcox.test(m$value[m$L1 == "Klf4_denovo"],m$value[m$L1 == "Klf4_known"])
mean(m$value[m$L1 == "Nanog_denovo"]); mean(m$value[m$L1 == "Nanog_Homer"])
wilcox.test(m$value[m$L1 == "Nanog_denovo"],m$value[m$L1 == "Nanog_Homer"])
mean(m$value[m$L1 == "Oct4_denovo"]); mean(m$value[m$L1 == "Oct4_known9"])
wilcox.test(m$value[m$L1 == "Oct4_denovo"],m$value[m$L1 == "Oct4_known9"])

w = melt(lapply(MotifMeth$WGBS, function(x) x$meth))

ggplot(w, aes(x = L1, y =value, fill = L1)) + geom_violin() +
  theme_classic() + 
  scale_fill_manual(values = c("red","red","blue","blue","darkgreen","darkgreen"))+
  labs(x = "Motif", y = "Methylation (%)", title = "WGBS")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

mean(w$value[w$L1 == "Klf4_denovo"]); mean(w$value[w$L1 == "Klf4_known"])
wilcox.test(w$value[w$L1 == "Klf4_denovo"],w$value[w$L1 == "Klf4_known"])
mean(w$value[w$L1 == "Nanog_denovo"]); mean(w$value[w$L1 == "Nanog_Homer"])
wilcox.test(w$value[w$L1 == "Nanog_denovo"],w$value[w$L1 == "Nanog_Homer"])
mean(w$value[w$L1 == "Oct4_denovo"]); mean(w$value[w$L1 == "Oct4_known9"])
wilcox.test(w$value[w$L1 == "Oct4_denovo"],w$value[w$L1 == "Oct4_known9"])

# Compare between WGBS and ATAC

m$Technique = "ATAC"
w$Technique = "WGBS"

all = rbind(m,w)

ggplot(all, aes(x = Technique, y = value, fill = Technique)) + geom_violin() + facet_wrap(~L1, nrow=3) + 
  scale_fill_manual(values = c("red","blue"))+ theme_classic()

mean(all$value[all$L1 == "Klf4_denovo" & all$Technique == "ATAC"]); mean(all$value[all$L1 == "Klf4_denovo" & all$Technique == "WGBS"])
mean(all$value[all$L1 == "Klf4_known" & all$Technique == "ATAC"]); mean(all$value[all$L1 == "Klf4_known" & all$Technique == "WGBS"])

mean(all$value[all$L1 == "Nanog_denovo" & all$Technique == "ATAC"]); mean(all$value[all$L1 == "Nanog_denovo" & all$Technique == "WGBS"])
mean(all$value[all$L1 == "Nanog_Homer" & all$Technique == "ATAC"]); mean(all$value[all$L1 == "Nanog_Homer" & all$Technique == "WGBS"])

mean(all$value[all$L1 == "Oct4_denovo" & all$Technique == "ATAC"]); mean(all$value[all$L1 == "Oct4_denovo" & all$Technique == "WGBS"])
mean(all$value[all$L1 == "Oct4_known9" & all$Technique == "ATAC"]); mean(all$value[all$L1 == "Oct4_known9" & all$Technique == "WGBS"])

wilcox.test(all$value[all$L1 == "Klf4_denovo" & all$Technique == "ATAC"],all$value[all$L1 == "Klf4_denovo" & all$Technique == "WGBS"])
wilcox.test(all$value[all$L1 == "Klf4_known" & all$Technique == "ATAC"],all$value[all$L1 == "Klf4_known" & all$Technique == "WGBS"])

wilcox.test(all$value[all$L1 == "Nanog_denovo" & all$Technique == "ATAC"],all$value[all$L1 == "Nanog_denovo" & all$Technique == "WGBS"])
wilcox.test(all$value[all$L1 == "Nanog_Homer" & all$Technique == "ATAC"],all$value[all$L1 == "Nanog_Homer" & all$Technique == "WGBS"])

wilcox.test(all$value[all$L1 == "Oct4_denovo" & all$Technique == "ATAC"],all$value[all$L1 == "Oct4_denovo" & all$Technique == "WGBS"])
wilcox.test(all$value[all$L1 == "Oct4_known9" & all$Technique == "ATAC"],all$value[all$L1 == "Oct4_known9" & all$Technique == "WGBS"])



sites = data.frame(
All_in_mm10 = unlist(lapply(MotifOcurrence, function(x) length(x@seqnames))),
All_in_ATAC_peaks = unlist(lapply(MotifMeth$ATAC, function(x) length(x@seqnames))),
All_in_WGBS_peaks = unlist(lapply(MotifMeth$WGBS, function(x) length(x@seqnames))))

write.table(sites, file = "/Volumes/DS/Franco/PriscilliaProject/Motifs/NumberOfSites.csv", quote = F, row.names = T, col.names = T)













