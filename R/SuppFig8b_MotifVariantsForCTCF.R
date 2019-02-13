
##############################################################################################
#  Authors : Franco Izz, Landau Lab, Weill Cornell Medicine #
##############################################################################################

library(methylKit)
library(genomation)
library(GenomicRanges)
library(readr)
library(ggplot2)
library(reshape2)

files = dir(path = "/Volumes/DS/Franco/PriscilliaProject/Motifs/", pattern = "fimo_")

motifs = list()
for(i in files){
  motifs[[i]] <- as.data.frame(read_delim(paste0("/Volumes/DS/Franco/PriscilliaProject/Motifs/",i), 
                           "\t", escape_double = FALSE, trim_ws = TRUE))
}
  chr = lapply(motifs, function(x){
  chr = strsplit(x = x$sequence_name, split = ":")
  unlist(lapply(chr, function(x) x[1]))
})
    
chrStart = lapply(motifs, function(x){
  chrStart = strsplit(x = x$sequence_name, split = ":")
  chrStart = unlist(lapply(chrStart, function(x) unlist(x[[2]])))
  chrStart = strsplit(chrStart, split = "-")
  chrStart = unlist(lapply(chrStart, function(x) x[1]))
})

chrEnd = lapply(motifs, function(x){
  chrEnd = strsplit(x = x$sequence_name, split = ":")
  chrEnd = unlist(lapply(chrEnd, function(x) unlist(x[[2]])))
  chrEnd = strsplit(chrEnd, split = "-")
  chrEnd = unlist(lapply(chrEnd, function(x) x[2]))
})


for(i in names(motifs)){
  motifs[[i]]$chr = chr[[i]]
  motifs[[i]]$chrStart = as.numeric(chrStart[[i]])
  motifs[[i]]$chrEnd = as.numeric(chrEnd[[i]])
}

for(i in names(motifs)){
  motifs[[i]]$motifStart = motifs[[i]]$chrStart + motifs[[i]]$start
  motifs[[i]]$motifEnd = motifs[[i]]$chrStart + motifs[[i]]$stop
}

motifRangesList = lapply(motifs, function(x){
  x = GRanges(seqnames = x$chr,
          ranges = IRanges(start = x$motifStart,end = x$motifEnd),
          strand = x$strand,
          score = x$score,
          p.value = x$`p-value`,
          q.value = x$`q-value`,
          sequence = x$matched_sequence)
  return(x)
})

motifs.starts = lapply(motifs, function(x){
  x = GRanges(seqnames = x$chr,
              ranges = IRanges(start = x$motifStart,end = x$motifStart),
              strand = x$strand,
              score = x$score,
              p.value = x$`p-value`,
              q.value = x$`q-value`,
              sequence = x$matched_sequence)
  return(x)
})


ctcf <- as.data.frame(read_delim("ctcf.merged.Peak_ctcf.merged.Cov__MERGED.tsv", 
                                                      "\t", escape_double = FALSE, trim_ws = TRUE))

ctcf = GRanges(seqnames = ctcf$peak.chr,
               ranges = IRanges(start = ctcf$meth.start, end = ctcf$meth.end),
               strand = "*",
               meth.percent = ctcf$meth.percent,
               meth.counts = ctcf$meth.counts,
               un.meth.counts = ctcf$un.meth.counts
               )

plot(density(ctcf$meth.counts+ctcf$un.meth.counts))
abline(v = 5, lty =2)

ctcf = ctcf[ctcf$meth.counts+ctcf$un.meth.counts >= 5,]

methylation = lapply(motifRangesList, function(x) subsetByOverlaps(ctcf,x))

distanceToStart = list()
for(i in names(motifs.starts)){
    distanceToStart[[i]] = distanceToNearest(x = methylation[[i]], subject = motifs.starts[[i]])
}

cpgCall = lapply(distanceToStart, function(x){
   x@elementMetadata$distance[x@elementMetadata$distance %in% c(2,3)] = "C2"
   x@elementMetadata$distance[x@elementMetadata$distance %in% c(13,14)] = "C12"
   return(x)
})

for(i in names(motifs.starts)){
motifs.starts[[i]]$cpgCall = cpgCall[[i]]@elementMetadata$distance[cpgCall[[i]]@to]
}

cpgsCaptured = as.data.frame(unlist(lapply(methylation, function(x) length(x))))
cpgsCaptured$Motif = rownames(cpgsCaptured)
colnames(cpgsCaptured) = c("CpGs","Motif")
cpgsCaptured = cpgsCaptured[order(cpgsCaptured$CpGs, decreasing = F),]
cpgsCaptured$Motif = gsub(cpgsCaptured$Motif, pattern = ".txt", replacement = "")
cpgsCaptured$Motif = gsub(cpgsCaptured$Motif, pattern = "fimo_", replacement = "")

cpgsCaptured$Motif = factor(cpgsCaptured$Motif, levels = cpgsCaptured$Motif)
write.table(x = cpgsCaptured, file = "/Volumes/DS/Franco/PriscilliaProject/NumberOfCpGsDetected.txt", quote = F, sep = "\t", col.names = T, row.names = T)

ggplot(cpgsCaptured, aes(x = Motif, y = CpGs)) + geom_bar(stat = "identity", width = 0.75, col = "black", fill = "black") +
  theme_classic(base_size = 14)+
  labs(y = "Number of CpGs detected", x = "CTCF motif variation")+
  coord_flip()+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))

methylationSub = lapply(methylation, function(x) x$meth.percent)
melted = melt(methylationSub)
colnames(melted) = c("meth.percent","Motif")
melted$Motif = gsub(melted$Motif, pattern = ".txt", replacement = "")
melted$Motif = gsub(melted$Motif, pattern = "fimo_", replacement = "")


melted$Motif = factor(x = melted$Motif, 
                      levels = c("CG_CA","CG_CT","CG_CC","CG_CG","CA_CG","CT_CG","CC_CG"))

ggplot(melted, aes(x = meth.percent, fill = Motif)) +
  geom_histogram(aes(y=..density..)) +
  facet_wrap(~Motif,nrow=1)+
  labs(y = "Frequency of methylation (density)", x = "Methylation (%)")+
  theme_classic()

ggplot(melted[melted$Motif %in% c("CG_CA","CG_CT","CG_CC","CG_CG"),], aes(x = meth.percent, color = Motif)) +
  geom_density() +
  labs(y = "Density", x = "Methylation (%)")+
  theme_classic()

ggplot(melted[melted$Motif %in% c("CG_CG","CA_CG","CT_CG","CC_CG"),], aes(x = meth.percent, color = Motif)) +
  geom_density() +
  labs(y = "Density", x = "Methylation (%)")+
  theme_classic()

pvals = list()
for(i in levels(melted$Motif)){
  for(j in levels(melted$Motif)){
      pvals[[i]][[j]] = wilcox.test(melted$meth.percent[melted$Motif == i],melted$meth.percent[melted$Motif == j])$p.value
}}
pvals = unlist(pvals)
write.table(as.data.frame(pvals), file = "/Volumes/DS/Franco/PriscilliaProject/pvaluesWilcoxonMotifVariants.txt", quote = F, sep = "\t", col.names = T, row.names = T)
