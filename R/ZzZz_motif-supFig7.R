##############################################################################################
#  Authors : Gunjan Sethia, Skok lab, Dept. Pathology, NYU Langone Health  #
##############################################################################################

#setwd("~/Google Drive/NYUMC/Skok lab/BiSulphiteGenome_analysis/June_2018_re-doing-some-plots/CTCF_motif_plots")

library(readr)

tsv <- data.frame(read_delim("~/Google Drive/NYUMC/Skok lab/BiSulphiteGenome_analysis/Jan_2018/motif-files/PL-ctcf-merged-table.tsv", 
                                   "\t", escape_double = FALSE, trim_ws = TRUE), stringsAsFactors = F)

cpg <- data.frame(read_delim("~/Google Drive/NYUMC/Skok lab/BiSulphiteGenome_analysis/Jan_2018/motif-files/CpG_context_CTCF-merged_q30_rmdup_sorted_by_read_name.txt", 
                             "\t", escape_double = FALSE, trim_ws = TRUE, col_names = F), stringsAsFactors = F)
names(cpg) <- c("read","strand", "chr", "start", "meth.status")
cpg$ID <- paste0(cpg$chr, ":", cpg$start)

tsv <- tsv[(tsv$position==5 | tsv$position==15), ]
tsv$ID <- paste0(tsv$meth.chr,":",tsv$meth.start)
tsv <- tsv[tsv$coverage > 5, ]

merged <- merge(tsv, cpg, by="ID")

dim(cpg)
dim(merged)
dim(tsv)

unique.reads <- unique(merged$read)

final <- data.frame()
for (i in 1:length(unique.reads)){
  #message (unique.reads[i])
  message (paste (i, "of", length(unique.reads)))
  tmp <- NULL
  tmp <- merged[merged$read==unique.reads[i], ]
  if (length(unique(tmp$position))==2){
    final <- data.frame(rbind(final, tmp), stringsAsFactors = F)
  }
}

#write.table(final, "final.tmp.tsv", sep="\t", row.names = F, col.names = T, quote = F)

final <- data.frame(read_delim("final.tmp.tsv", 
                               "\t", escape_double = FALSE, trim_ws = TRUE),stringsAsFactors = F)

data <- unique(final[,c(23,22,25:27)])

#data <- data[data$read!="D00595:307:CAK5VANXX:3:2210:19185:3497_1:N:0:AGTTGGAA",]

unique.reads <- unique(data$read)

final.2 <- data.frame()
final.3 <- data.frame()
final.4 <- data.frame()

for (i in 1:length(unique.reads)){
#for (i in 5159:length(unique.reads)){
  tmp <- NULL
  tmp <- data[data$read==unique.reads[i], ]
  if (nrow(tmp)==2){
    message (paste (i, "of", length(unique.reads)))
    message ("Its a 2...")
    print (tmp)
    holder <- data.frame(cbind(tmp[tmp$position==5,], tmp[tmp$position==15,]), stringsAsFactors = F)
    final.2 <- data.frame(rbind(final.2, holder), stringsAsFactors = F)
  }
  
  if (nrow(tmp)==3){
    message (paste (i, "of", length(unique.reads)))
    message ("Its a 3...")
    #print (tmp)
    
    diff <- abs(c(tmp$start[1]-tmp$start[2], tmp$start[2]-tmp$start[3], tmp$start[1]-tmp$start[3]))
    
    if (10 %in% diff) {
      print ("5:15 on same read")
      x <- tmp
      freq <- data.frame(table(tmp$position), stringsAsFactors = F)
      uniq <- freq$Var1[freq$Freq==1]
      holder <- tmp[tmp$position==uniq,]
      holder <- data.frame(rbind(holder, tmp[abs(holder$start-tmp$start)==10,]), stringsAsFactors = F)
      holder <- holder[order(holder$position),]
      holder <- data.frame(cbind(holder[holder$position==5,], holder[holder$position==15,]), stringsAsFactors = F)
      #print ("-----------------")
      #print (holder)
      final.3 <- data.frame(rbind(final.3, holder), stringsAsFactors = F)
      
    }
    
    if ((10 %in% diff)==FALSE){
      print ("PASS..nothing on same read")
    }
    
    

  }

  if (nrow(tmp)==4){
    message (paste (i, "of", length(unique.reads)))
    message ("Its a 4...")
    #print (tmp)
    x <- tmp
    tmp <- tmp[order(tmp$start),]
    
    if (length(unique(data.frame(table(tmp$position))$Freq))==1){
      message ("4...but unique")
      h12 <- tmp[1:2,]
      h12 <- data.frame(cbind(h12[h12$position==5,], h12[h12$position==15,]), stringsAsFactors = F)
      h34 <- tmp[3:4,]
      h34 <- data.frame(cbind(h34[h34$position==5,], h34[h34$position==15,]), stringsAsFactors = F)
      holder <- data.frame(rbind(h12,h34), stringsAsFactors = F)
      final.4 <- data.frame(rbind(final.4, holder), stringsAsFactors = F)
    }
    
    if (length(unique(data.frame(table(tmp$position))$Freq))!=1){
      message ("4...but NOT unique")
      freq <- data.frame(table(tmp$position))
      n <- freq$Var1[freq$Freq==1]
      x1 <- tmp[tmp$start==(tmp$start[abs(tmp$start[tmp$position==n]-tmp$start)==10]),]
      x2 <- tmp[tmp$position==n,]
      tmp1 <- data.frame(rbind(x1,x2), stringsAsFactors = F)
      holder <- data.frame(cbind(tmp1[tmp1$position==5,] , tmp1[tmp1$position==15,] ), stringsAsFactors = F)
      final.4 <- data.frame(rbind(final.4, holder), stringsAsFactors = F)
      
    }
  }
 
}


df <- data.frame(rbind(final.2, final.3, final.4), stringsAsFactors = F)
df$combination <- paste0(df$meth.status, df$meth.status.1)

#write.table(df, "ZzZz_motif_FINAL.tsv", sep="\t", row.names = F, col.names = T, quote = F)


for (i in 1:length(unique.reads)){
  tmp <- NULL
  tmp <- data[data$read==unique.reads[i], ]
  if (nrow(tmp)>4){
    print (tmp)
  }
  
}


observed <- data.frame((table(df$combination)), stringsAsFactors = F)
names(observed) <- c("meth5.meth15", "counts")
observed$percent <- (observed$counts/sum(observed$counts))*100
observed

data.frame(table(df$meth.status))
data.frame(table(df$meth.status.1))

