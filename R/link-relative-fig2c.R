##############################################################################################
#  Authors : Gunjan Sethia, Skok lab, Dept. Pathology, NYU Langone Health  #
##############################################################################################

# get gene coordinates
counts <- data.frame(read_delim("~/Google Drive/NYUMC/Skok lab/BiSulphiteGenome_analysis/Jan_2018/gene_expression/RNAseq-mESC-merged.TPM.unstr.txt", "\t", escape_double = FALSE, col_names = TRUE,  trim_ws = TRUE), stringsAsFactors = F)
colnames(counts) <- c("Gene","counts")
counts <- data.frame(cSplit(counts, "Gene", ":"), stringsAsFactors = F)[,c(2,1)]
colnames(counts) <- c("Gene","counts")
head(counts)
dim(counts)

genes <- data.frame(read_delim("~/Skok_Lab/EpiMethylTag_pipeline/ref_files/mm10_genes.bed", 
                               "\t", escape_double = FALSE, col_names = FALSE, 
                               trim_ws = TRUE), stringsAsFactors = F)
head(genes)
names(genes)[4] <- "Gene"

merged <- merge(counts, genes, by = "Gene")
merged <- merged[complete.cases(merged),]
dim(merged)
merged <- merged[!(duplicated(merged)),]
head(merged)
dim(merged)
merged$ID <- paste0(merged$X1,":", merged$X2, ":", merged$X3)

gene.coord <- merged[,c(3:5)]
write.table(gene.coord, "matac-gene-tpm-coord.bed", col.names = F, row.names = F, quote = F, sep="\t")
