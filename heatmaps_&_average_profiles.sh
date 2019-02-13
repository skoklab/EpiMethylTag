
##############################################################################################
#  Authors : Priscillia Lhoumaud & Gunjan Sethia, Skok lab, Dept. Pathology, NYU Langone Health #
##############################################################################################

# Fig. 2d and Supplementary Fig. 5b
# extract 2 groups:
#   1.promoter
#   2.non promoter: intragenic (exon, intron, 5UTR, 3UTR) + intergenic
# 4 groups of ATAC CpGs based on methylation and coverage (cov <5 reads removed)
#   1-meth_00to20-cov_50andabove
#   2-meth_00to20-cov_5to50
#   3-meth_20to80
#   4-meth_80to100
# Path on cluster- /ifs/home/lhoump01/Method/deeptools_heatmaps/genomic-annotations

################## loading module ################
module unload python
module unload deeptools
module load deeptools/2.3.3

################## bigWig files ################
# BW_1="/ifs/home/lhoump01/Method/deeptools_heatmaps/ATAC-merged.deduplicated.sort.Q30.bw"
# BW_2="/ifs/home/lhoump01/mESC_ChIP-BS/Public/K4me1_GSM1000121_q30_rmdup_sorted.bw"
# BW_2="/ifs/home/lhoump01/mESC_ChIP-BS/Public/K27ac_GSM1000126_q30_rmdup_sorted.bw"
# BW_2="/ifs/home/lhoump01/mESC_ChIP-BS/Public/K4me3_GSM1000124_q30_rmdup_sorted.bw"

################## for promoter ################
cat 1-promoter.bed file1.txt 2-promoter.bed file2.txt 3-promoter.bed file3.txt 4-promoter.bed file4.txt > 4groups_Promoter.bed

#Try using sort by ATAC signal then run in one block:
BED="/ifs/home/lhoump01/Method/deeptools_heatmaps/genomic-annotations/4groups_Promoter.bed"

group_plot(){
computeMatrix reference-point --referencePoint center --scoreFileName ${BW_1} --regionsFileName ${BED} --outFileName ${OUTPUT_FILE_PREFIX}_matrix.gz --sortRegions descend --sortUsingSamples 1 --missingDataAsZero --beforeRegionStartLength ${Midpoint} --afterRegionStartLength ${Midpoint}
 
plotHeatmap --matrixFile ${OUTPUT_FILE_PREFIX}_matrix.gz --outFileName ${OUTPUT_FILE_PREFIX}_heatmap.pdf --colorMap ${color} --sortRegions descend --sortUsingSamples 1 --heatmapHeight 50 --legendLocation none --whatToShow 'heatmap and colorbar' --zMin ${zMin} --zMax ${zMax}
 
plotProfile --matrixFile ${OUTPUT_FILE_PREFIX}_matrix.gz --outFileName ${OUTPUT_FILE_PREFIX}_profilePlot.pdf --colors red blue green purple --legendLocation 'upper-center' --yMin ${y_scale_min} --yMax ${y_scale_max}
 
plotProfile --matrixFile ${OUTPUT_FILE_PREFIX}_matrix.gz --outFileName ${OUTPUT_FILE_PREFIX}_pergroupprofilePlot.pdf --colors red --legendLocation 'upper-center' --yMin ${y_scale_min} --yMax ${y_scale_max} --perGroup
}
 
#then run in one bloc for ATAC:
RUNDIR=/ifs/home/lhoump01/Method/deeptools_heatmaps
Midpoint=1000
y_scale_min=0
y_scale_max=150
color=Oranges
zMin=0
zMax=75
BW_1="/ifs/home/lhoump01/Method/deeptools_heatmaps/ATAC-merged.deduplicated.sort.Q30.bw"
cd ${RUNDIR}
OUTPUT_FILE_PREFIX="ATAC_meth-cov_4groups_ATAC_promoter"
group_plot
 
#### for histone marks
group_plot(){
computeMatrix reference-point --referencePoint center --scoreFileName ${BW_1} ${BW_2} --regionsFileName ${BED} --outFileName ${OUTPUT_FILE_PREFIX}_matrix.gz --sortRegions descend --sortUsingSamples 1 --missingDataAsZero --beforeRegionStartLength ${Midpoint} --afterRegionStartLength ${Midpoint}
 
plotHeatmap --matrixFile ${OUTPUT_FILE_PREFIX}_matrix.gz --outFileName ${OUTPUT_FILE_PREFIX}_heatmap.pdf --colorMap ${color} --sortRegions descend --sortUsingSamples 1 --heatmapHeight 50 --legendLocation none --whatToShow 'heatmap and colorbar' --zMin ${zMin} --zMax ${zMax}
 
plotProfile --matrixFile ${OUTPUT_FILE_PREFIX}_matrix.gz --outFileName ${OUTPUT_FILE_PREFIX}_profilePlot.pdf --colors red blue green purple --legendLocation 'upper-center' --yMin ${y_scale_min} --yMax ${y_scale_max}
 
plotProfile --matrixFile ${OUTPUT_FILE_PREFIX}_matrix.gz --outFileName ${OUTPUT_FILE_PREFIX}_pergroupprofilePlot.pdf --colors red blue --legendLocation 'upper-center' --yMin ${y_scale_min} --yMax ${y_scale_max} --perGroup
}
  
#then run in one bloc for K4me1:
BW_1="/ifs/home/lhoump01/Method/deeptools_heatmaps/ATAC-merged.deduplicated.sort.Q30.bw"
BW_2="/ifs/home/lhoump01/mESC_ChIP-BS/Public/K4me1_GSM1000121_q30_rmdup_sorted.bw"
RUNDIR=/ifs/home/lhoump01/Method/deeptools_heatmaps
Midpoint=1000
y_scale_min=0
y_scale_max=0.20
color=Greens
zMin=0
zMax=0.25
cd ${RUNDIR}
OUTPUT_FILE_PREFIX="ATAC_meth-cov_4groups_K4me1_promoter"
group_plot
 
#then run in one bloc for K27ac:
RUNDIR=/ifs/home/lhoump01/Method/deeptools_heatmaps
BW_1="/ifs/home/lhoump01/Method/deeptools_heatmaps/ATAC-merged.deduplicated.sort.Q30.bw"
BW_2="/ifs/home/lhoump01/mESC_ChIP-BS/Public/K27ac_GSM1000126_q30_rmdup_sorted.bw"
Midpoint=1000
y_scale_min=0
y_scale_max=0.4
color=Greens
zMin=0
zMax=0.5
cd ${RUNDIR}
OUTPUT_FILE_PREFIX="ATAC_meth-cov_4groups_K27ac_promoter"
group_plot
 
#then run in one bloc for K4me3:
BW_1="/ifs/home/lhoump01/Method/deeptools_heatmaps/ATAC-merged.deduplicated.sort.Q30.bw"
BW_2="/ifs/home/lhoump01/mESC_ChIP-BS/Public/K4me3_GSM1000124_q30_rmdup_sorted.bw"
RUNDIR=/ifs/home/lhoump01/Method/deeptools_heatmaps
Midpoint=1000
y_scale_min=0
y_scale_max=2
color=Greens
zMin=0
zMax=1
cd ${RUNDIR}
OUTPUT_FILE_PREFIX="ATAC_meth-cov_4groups_K4me3_promoter"
group_plot

 
################## for inter_intra ################
### combine intergenic and intragenic together ###
 
# Path on cluster/ifs/home/lhoump01/Method/deeptools_heatmaps/genomic-annotations
 cat 1-Intragenic.bed 1-Intergenic.bed file1.txt 2-Intragenic.bed 2-Intergenic.bed file2.txt 3-Intragenic.bed 3-Intergenic.bed file3.txt 4-Intragenic.bed 4-Intergenic.bed file4.txt > 4groups_InterIntra.bed

BED="/ifs/home/lhoump01/Method/deeptools_heatmaps/genomic-annotations/4groups_InterIntra.bed"
 
################### sort by ATAC signal
#then run in one bloc:

group_plot(){
computeMatrix reference-point --referencePoint center --scoreFileName ${BW_1} --regionsFileName ${BED} --outFileName ${OUTPUT_FILE_PREFIX}_matrix.gz --sortRegions descend --sortUsingSamples 1 --missingDataAsZero --beforeRegionStartLength ${Midpoint} --afterRegionStartLength ${Midpoint}
 
plotHeatmap --matrixFile ${OUTPUT_FILE_PREFIX}_matrix.gz --outFileName ${OUTPUT_FILE_PREFIX}_heatmap.pdf --colorMap ${color} --sortRegions descend --sortUsingSamples 1 --heatmapHeight 50 --legendLocation none --whatToShow 'heatmap and colorbar' --zMin ${zMin} --zMax ${zMax}
 
plotProfile --matrixFile ${OUTPUT_FILE_PREFIX}_matrix.gz --outFileName ${OUTPUT_FILE_PREFIX}_profilePlot.pdf --colors red blue green purple --legendLocation 'upper-center' --yMin ${y_scale_min} --yMax ${y_scale_max}
 
#plotProfile --matrixFile ${OUTPUT_FILE_PREFIX}_matrix.gz --outFileName ${OUTPUT_FILE_PREFIX}_pergroupprofilePlot.pdf --colors red --legendLocation 'upper-center' --yMin ${y_scale_min} --yMax ${y_scale_max} --perGroup
}
 
#then run in one bloc for ATAC:
RUNDIR=/ifs/home/lhoump01/Method/deeptools_heatmaps
Midpoint=1000
y_scale_min=0
y_scale_max=150
color=Oranges
zMin=0
zMax=75
BW_1="/ifs/home/lhoump01/Method/deeptools_heatmaps/ATAC-merged.deduplicated.sort.Q30.bw"
cd ${RUNDIR}
OUTPUT_FILE_PREFIX="ATAC_meth-cov_4groups_ATAC_InterIntra"
group_plot
 
#### for histone marks
group_plot(){
computeMatrix reference-point --referencePoint center --scoreFileName ${BW_1} ${BW_2} --regionsFileName ${BED} --outFileName ${OUTPUT_FILE_PREFIX}_matrix.gz --sortRegions descend --sortUsingSamples 1 --missingDataAsZero --beforeRegionStartLength ${Midpoint} --afterRegionStartLength ${Midpoint}
 
plotHeatmap --matrixFile ${OUTPUT_FILE_PREFIX}_matrix.gz --outFileName ${OUTPUT_FILE_PREFIX}_heatmap.pdf --colorMap ${color} --sortRegions descend --sortUsingSamples 1 --heatmapHeight 50 --legendLocation none --whatToShow 'heatmap and colorbar' --zMin ${zMin} --zMax ${zMax}
 
plotProfile --matrixFile ${OUTPUT_FILE_PREFIX}_matrix.gz --outFileName ${OUTPUT_FILE_PREFIX}_profilePlot.pdf --colors red blue green purple --legendLocation 'upper-center' --yMin ${y_scale_min} --yMax ${y_scale_max}
 
#plotProfile --matrixFile ${OUTPUT_FILE_PREFIX}_matrix.gz --outFileName ${OUTPUT_FILE_PREFIX}_pergroupprofilePlot.pdf --colors red blue --legendLocation 'upper-center' --yMin ${y_scale_min} --yMax ${y_scale_max} --perGroup
}
 
 
#then run in one bloc for K4me1:
RUNDIR=/ifs/home/lhoump01/Method/deeptools_heatmaps
Midpoint=1000
y_scale_min=0
y_scale_max=0.15
color=Greens
zMin=0
zMax=0.25
BW_1="/ifs/home/lhoump01/Method/deeptools_heatmaps/ATAC-merged.deduplicated.sort.Q30.bw"
BW_2="/ifs/home/lhoump01/mESC_ChIP-BS/Public/K4me1_GSM1000121_q30_rmdup_sorted.bw"
cd ${RUNDIR}
OUTPUT_FILE_PREFIX="ATAC_meth-cov_4groups_K4me1_InterIntra"
group_plot
 
#then run in one bloc for K27ac:
RUNDIR=/ifs/home/lhoump01/Method/deeptools_heatmaps
Midpoint=1000
y_scale_min=0
y_scale_max=0.4
color=Greens
zMin=0
zMax=0.5
BW_1="/ifs/home/lhoump01/Method/deeptools_heatmaps/ATAC-merged.deduplicated.sort.Q30.bw"
BW_2="/ifs/home/lhoump01/mESC_ChIP-BS/Public/K27ac_GSM1000126_q30_rmdup_sorted.bw"
cd ${RUNDIR}
OUTPUT_FILE_PREFIX="ATAC_meth-cov_4groups_K27ac_InterIntra"
group_plot
 
#then run in one bloc for K4me3:
RUNDIR=/ifs/home/lhoump01/Method/deeptools_heatmaps
Midpoint=1000
y_scale_min=0
y_scale_max=0.5
color=Greens
zMin=0
zMax=0.5
BW_1="/ifs/home/lhoump01/Method/deeptools_heatmaps/ATAC-merged.deduplicated.sort.Q30.bw"
BW_2="/ifs/home/lhoump01/mESC_ChIP-BS/Public/K4me3_GSM1000124_q30_rmdup_sorted.bw"
cd ${RUNDIR}
OUTPUT_FILE_PREFIX="ATAC_meth-cov_4groups_K4me3_InterIntra"
group_plot



#Fig. 3c
#Of note: location of CTCF motif (MA0139.1) into each cytosine group (as per the fig. 2a) were obtained using FIMO on MEME suite.

#### merge CTCF motifs from the 4 groups
cat 1-CpGs0to20_cov50andabove_intoCTCFmotif_uniq.bed file1.txt 2-CpGs0to20_cov5to50_intoCTCFmotif_uniq.bed file2.txt 3-CpGs20to80intoCTCFmotif_uniq.bed file3.txt 4-CpGs80to100intoCTCFmotif_uniq.bed file4.txt > fimoCTCF_4groups.bed

BED="/ifs/home/lhoump01/Method/deeptools_heatmaps/fimoCTCF_4groups.bed"

group_plot(){
  computeMatrix reference-point --referencePoint center --scoreFileName ${BW_1} ${BW_2} ${BW_3} --regionsFileName ${BED} --outFileName ${OUTPUT_FILE_PREFIX}_matrix.gz --sortRegions descend --sortUsingSamples 1 --missingDataAsZero --beforeRegionStartLength ${Midpoint} --afterRegionStartLength ${Midpoint}
  plotHeatmap --matrixFile ${OUTPUT_FILE_PREFIX}_matrix.gz --outFileName ${OUTPUT_FILE_PREFIX}_heatmap.pdf --colorMap ${color} --sortRegions descend --sortUsingSamples 1 --heatmapHeight 50 --legendLocation none --whatToShow 'heatmap and colorbar' --zMin ${zMin} --zMax ${zMax}
  plotProfile --matrixFile ${OUTPUT_FILE_PREFIX}_matrix.gz --outFileName ${OUTPUT_FILE_PREFIX}_profilePlot.pdf --colors red blue green purple --legendLocation 'upper-center' --yMin ${y_scale_min} --yMax ${y_scale_max}

  plotProfile --matrixFile ${OUTPUT_FILE_PREFIX}_matrix.gz --outFileName ${OUTPUT_FILE_PREFIX}_pergroupprofilePlot.pdf --colors red blue green --legendLocation 'upper-center' --yMin ${y_scale_min} --yMax ${y_scale_max}  --perGroup
}

############# Plot CTCF
#then run in one bloc:
RUNDIR=/ifs/home/lhoump01/Method/deeptools_heatmaps
Midpoint=500
y_scale_min=0
y_scale_max=15
color=Blues
zMin=0
zMax=15
BW_1="/ifs/home/lhoump01/Method/deeptools_heatmaps/ATAC-merged.deduplicated.sort.Q30.bw"
BW_2="/ifs/home/lhoump01/Method/deeptools_heatmaps/CTCF-merged.deduplicated.sort.Q30.bw"
BW_3="/ifs/home/lhoump01/Method/deeptools_heatmaps/mES_CTCF.deduplicated.sort.Q30.bw"
cd ${RUNDIR}
OUTPUT_FILE_PREFIX="ATACpeaks_4groups_CTCFmotifs_CTCF"
group_plot

### plot ATAC
RUNDIR=/ifs/home/lhoump01/Method/deeptools_heatmaps
Midpoint=500
y_scale_min=0
y_scale_max=200
color=Oranges
zMin=0
zMax=75
BW_1="/ifs/home/lhoump01/Method/deeptools_heatmaps/ATAC-merged.deduplicated.sort.Q30.bw"
BW_2="/ifs/home/lhoump01/Method/deeptools_heatmaps/CTCF-merged.deduplicated.sort.Q30.bw"
BW_3="/ifs/home/lhoump01/Method/deeptools_heatmaps/mES_CTCF.deduplicated.sort.Q30.bw"
cd ${RUNDIR}
OUTPUT_FILE_PREFIX="ATACpeaks_4groups_CTCFmotifs_ATAC"
group_plot

## END
