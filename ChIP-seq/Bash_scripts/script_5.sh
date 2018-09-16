############################
# ChIP-sequencing analysis #
############################
# Create directory in home directory (ChIP-Seq_$date) a directory of RNAseq that must contain the file of RNAseq data and 
# the file of the Reference genome in hg19 (RefSeq_Genes.hg19.strand.bed)


# HeLa RNAseq expression groups
# For ChIP-seq average signal by groups of gene expression, files of NM must be created
cd RNAseq
awk '{print $2, $4}' HeLa_rep2_10882.final.rpkm | sed 's/ /\t/g' | sort -k2,2n > HeLaS3_RNAseq_sorted.rpkm
# keep gene name and RPKM value. Then, sort by RPKM value

awk 'FILENAME=="RefSeq_Genes.hg19.strand.bed" {A[$6]=$0} FILENAME=="HeLaS3_RNAseq_sorted.rpkm" {if ($1 in A) {print $0, A[$1]}}' RefSeq_Genes.hg19.strand.bed HeLaS3_RNAseq_sorted.rpkm | sed 's/ /\t/g' | awk '{print $3, $4, $5, $6, $1, $2, $7}' | sed 's/ /\t/g' | sort -k6,6n > HeLaS3_RNAseq_sorted.strand.bed
# RefSeq file is a Bed file from all genes in hg19. From this file, I take columns of interest
# and bind them (by gene name coincidence) to my RNA seq file I have from previous line of code.

awk '{if ($6==0) {print $0}}' HeLaS3_RNAseq_sorted.strand.bed | sed 's/ /\t/g' > ExpressionGroups/HeLaS3_RNAseq_sorted_NotExpressed.strand.bed
awk '{if ($6!=0) {print $0}}' HeLaS3_RNAseq_sorted.strand.bed | sed 's/ /\t/g' > ExpressionGroups/HeLaS3_RNAseq_sorted_Expressed.strand.bed
# Divide my file with all information into two bed files; one containing all genes not expressed (RPKM value = 0)
# and another bed file with genes that are expressed.


cd ExpressionGroups
split -l 1279 --additional-suffix=_HeLaS3_RNAseq_sorted_Expressed.strand.bed HeLaS3_RNAseq_sorted_Expressed.strand.bed
# divide file with expressed genes into 10 groups with the same number of genes per group.

mkdir NM_files
for file in *.bed
do
awk '{print $7}' $file > NM_files/"$(basename "$file")_NM"
done
# to keep NM identification for all of the 10 groups created.

# Plot in R. Boxplot of expression groups can be found in script '2.BoxplotExpression.R'
