############################
# ChIP-sequencing analysis #
############################

# Further analyses
mkdir CEAS

## some settings
CeasRefGene=/home/andrea/Softwares/hg19_ceas.refGene

################################################
# ChIP-seq average signal by expression groups #
################################################

# Expression data comes from RNAseq (See directory RNAseq). 

cd CEAS
for file in ../wigFiles/*inputMerged_subtracted.wig
do
ceas --name="$(basename "$file")_CEAS_all_groups" --pf-res=50 --gn-groups=../../RNAseq/ExpressionGroups/NM_files/HeLaS3_RNAseq_sorted_NotExpressed.strand.bed_NM,../../RNAseq/ExpressionGroups/NM_files/xaa_HeLaS3_RNAseq_sorted_Expressed.strand.bed_NM,../../RNAseq/ExpressionGroups/NM_files/xab_HeLaS3_RNAseq_sorted_Expressed.strand.bed_NM,../../RNAseq/ExpressionGroups/NM_files/xac_HeLaS3_RNAseq_sorted_Expressed.strand.bed_NM,../../RNAseq/ExpressionGroups/NM_files/xad_HeLaS3_RNAseq_sorted_Expressed.strand.bed_NM,../../RNAseq/ExpressionGroups/NM_files/xae_HeLaS3_RNAseq_sorted_Expressed.strand.bed_NM,../../RNAseq/ExpressionGroups/NM_files/xaf_HeLaS3_RNAseq_sorted_Expressed.strand.bed_NM,../../RNAseq/ExpressionGroups/NM_files/xag_HeLaS3_RNAseq_sorted_Expressed.strand.bed_NM,../../RNAseq/ExpressionGroups/NM_files/xah_HeLaS3_RNAseq_sorted_Expressed.strand.bed_NM,../../RNAseq/ExpressionGroups/NM_files/xai_HeLaS3_RNAseq_sorted_Expressed.strand.bed_NM,../../RNAseq/ExpressionGroups/NM_files/xaj_HeLaS3_RNAseq_sorted_Expressed.strand.bed_NM --gn-group-names='Not Expr,1,2,3,4,5,6,7,8,9,10' -g $CeasRefGene -w $file
done 
cd ..

############################
# H1 abundance at the LADs #
############################
# Create directory in home directory (ChIP-Seq_$date) a directory called LADs that must contain the files of LADs data 

# 1. LADs coordinated from UCSC
# Sort by chromosome and position LADs and bdg files (this latter one is already done)
cd LADs
sort -k1,1 -k2,2n LADS.hg19.converted.bed > tmp
mv tmp LADS.hg19.converted.bed
cd ..

for file in bdgFiles/*bdg
do
bedtools map  -a LADs/LADS.hg19.converted.bed -b $file -c 4 -o mean > LADs/"$(basename "$file")$.LADsOverlap.bed"
done

# bedtools map allows to map overlapping features in B file onto features in A file and
# to compute the average score of the bdg file that overlap at LADs start and end positions.


for file in LADs/*Overlap.bed
do
awk '{print $5}' $file > "$(basename "$file")$.mean"
done

###################################
# H1 abundance at the centromeres #
###################################
# Centromeres position plus down/upstream regions if they are closed enough to the centromere.
# Create directory in home directory (ChIP-Seq_$date) a directory called Centromeres that must contain the files of centromere position data 

# 1. Create directory
mkdir Centromeres
cd Centromeres
# 2. tab delimiter and +50/-50 bp inside centromere positions
awk '{print $0}' centromeres_position_ext.bed | sed 's/ /\t/g' > centromeres_position.bed
awk '{print $1, $2+50, $3-50}' centromeres_position.bed | sed 's/ /\t/g' > centromeres_position_ext_50.bed

# extend centromeric positions 50 bp down- and upstream

# 3. p arm (spand 25000)
awk '{print $1, $2-25000, $2}' centromeres_position_ext_50.bed | sed 's/ /\t/g' > p_centrom.bed
sort -k1,1 -k2,2n p_centrom.bed > tmp
mv tmp p_centrom.bed

# spand centromere regions 25000 bp down- and upstream to define pericentromeric regions.

# 3.2 remove acrocentric chromosomes
sed -e '5,7d;14,15d;24d' p_centrom.bed | sed 's/ /\t/g' > tmp
mv tmp p_centrom.bed

# 4. q arm
awk '{print $1, $3, $3+25000}' centromeres_position_ext_50.bed | sed 's/ /\t/g' > q_centrom.bed
sort -k1,1 -k2,2n q_centrom.bed > tmp
mv tmp q_centrom.bed
cd ..

# 5. bedtools map with our ChIP-Signal of the different variants
## p arm
for file in bdgFiles/*bdg
do
bedtools map  -a Centromeres/p_centrom.bed -b $file -c 4 -o mean > Centromeres/"$(basename "$file")$.CentrOverlap.bed"
done

## q arm
for file in bdgFiles/*bdg
do
bedtools map  -a Centromeres/q_centrom.bed -b $file -c 4 -o mean > Centromeres/"$(basename "$file")$.CentrOverlap.bed"
done

###################################
# H1 abundance at FAIRE-Seq peaks #
###################################

# Create directory in home directory (ChIP-Seq_$date) a directory called FAIRE-seq that must contain the peak file (bed file)
mkdir FAIRE-seq

for file in bdgFiles/*bdg
do
bedtools map  -a FAIRE-seq/ENCFF001UYM.bed -b $file -c 4 -o mean > FAIRE-seq/"$(basename "$file")$.FAIRE_peaks_overlap"
done

for file in FAIRE-seq/*overlap
do
awk '{print $11}' $file | sed 's/ /\t/g' > "$(basename "$file")$.mean"
done


#############################
# H1 abundance at enhancers #
#############################

# Create directory in home directory (ChIP-Seq_$date) a directory called enhancers that must contain the enhancer positions (bed file)
mkdir enhancers

for file in bdgFiles/*bdg
do
bedtools map  -a enhancers/enhancers.bed -b $file -c 4 -o mean > enhancers/"$(basename "$file")$.enhancers_peaks_overlap"
done

# PLOTS DONE WITH R
# look in R script called '3.AbundSpecificRegions.R'

##############################################################
# Occupancy of H1 variants at specific histone modifications #
##############################################################

# Create directory in home directory (ChIP-Seq_$date) a directory called PTMs that must contain the PTMs peak files (bed file)
mkdir PTMs
cd PTMs

# H3K27ac_GSM733684/
sort -k1,1 -k2,2n GSM733684_hg19_wgEncodeBroadHistoneHelas3H3k27acStdPk.broadPeak > tmp
mv tmp GSM733684_hg19_wgEncodeBroadHistoneHelas3H3k27acStdPk.broadPeak

for file in ../bdgFiles/*bdg
do
bedtools map -a GSM733684_hg19_wgEncodeBroadHistoneHelas3H3k27acStdPk.broadPeak -b $file -c 4 -o mean > "H3K27ac$(basename "$file")"
done

#H3K27me3_GSM733696
sort -k1,1 -k2,2n GSM733696_hg19_wgEncodeBroadHistoneHelas3H3k27me3StdPk.broadPeak > tmp
mv tmp GSM733696_hg19_wgEncodeBroadHistoneHelas3H3k27me3StdPk.broadPeak

for file in ../bdgFiles/*bdg
do
bedtools map -a GSM733696_hg19_wgEncodeBroadHistoneHelas3H3k27me3StdPk.broadPeak -b $file -c 4 -o mean > "H3K27me3$(basename "$file")"
done

#H3K4me3_GSM733682
sort -k1,1 -k2,2n GSM733682_hg19_wgEncodeBroadHistoneHelas3H3k4me3StdPk.broadPeak > tmp
mv tmp GSM733682_hg19_wgEncodeBroadHistoneHelas3H3k4me3StdPk.broadPeak

for file in ../bdgFiles/*bdg
do
bedtools map -a GSM733682_hg19_wgEncodeBroadHistoneHelas3H3k4me3StdPk.broadPeak -b $file -c 4 -o mean > "H3K4me3$(basename "$file")"
done

#/H3K9me3
sort -k1,1 -k2,2n GSM1003480_hg19_wgEncodeBroadHistoneHelas3H3k09me3Pk.broadPeak > tmp
mv tmp GSM1003480_hg19_wgEncodeBroadHistoneHelas3H3k09me3Pk.broadPeak

for file in cd PTMsbdgFiles/*bdg
do
bedtools map -a GSM1003480_hg19_wgEncodeBroadHistoneHelas3H3k09me3Pk.broadPeak -b $file -c 4 -o mean > "H3K9me3$(basename "$file")"
done

# H3K4me1_GSM798322/
sort -k1,1 -k2,2n GSM798322_hg19_wgEncodeBroadHistoneHelas3H3k04me1StdPk.broadPeak > tmp
mv tmp GSM798322_hg19_wgEncodeBroadHistoneHelas3H3k04me1StdPk.broadPeak

for file in ../bdgFiles/*bdg
do
bedtools map -a GSM798322_hg19_wgEncodeBroadHistoneHelas3H3k04me1StdPk.broadPeak -b $file -c 4 -o mean > "H3K4me1$(basename "$file")"
done

#H3K36me3_GSM733711
sort -k1,1 -k2,2n GSM733711_hg19_wgEncodeBroadHistoneHelas3H3k36me3StdPk.broadPeak > tmp
mv tmp GSM733711_hg19_wgEncodeBroadHistoneHelas3H3k36me3StdPk.broadPeak

for file in ../bdgFiles/*bdg
do
bedtools map -a GSM733711_hg19_wgEncodeBroadHistoneHelas3H3k36me3StdPk.broadPeak -b $file -c 4 -o mean > "H3K36me3$(basename "$file")"
done

#H3K9ac_GSM733756
sort -k1,1 -k2,2n GSM733756_hg19_wgEncodeBroadHistoneHelas3H3k9acStdPk.broadPeak > tmp
mv tmp GSM733756_hg19_wgEncodeBroadHistoneHelas3H3k9acStdPk.broadPeak

for file in ../bdgFiles/*bdg
do
bedtools map -a GSM733756_hg19_wgEncodeBroadHistoneHelas3H3k9acStdPk.broadPeak -b $file -c 4 -o mean > "H3K9ac$(basename "$file")"
done
#
cd ..

# PLOTS DONE WITH R
# look in R script called '4.OccupAtPTMs.R'


# Average ChIP signal profiling of H1 variants was performed for other genomic features

############################
# RNA Pol II binding sites #
############################
# Mean ChIP-seq signal of H1 variants around the center of RNA Pol II binding sites

# Create directory in home directory (ChIP-Seq_$date) a directory of RNAPOLL_II that must contain the files of RNAPOL II peaks data 

for file in RNAPOL_II/*Peak
do
awk '{print $1, $2, $3}' $file | sed 's/ /\t/g' > RNAPOL_II/"$(basename "$file")_for_sitepro"
done

cd RNAPOL_II
for file in *sitepro
do
sitepro --name="$(basename "$file")" --pf-res=50 -b $file -w ../wigFiles/H1_2.inputMerged_subtracted.wig -w ../wigFiles/H1X.inputMerged_subtracted.wig -w ../wigFiles/2_HA.inputMerged_subtracted.wig -w ../wigFiles/0_HA.inputMerged_subtracted.wig
done
cd ..

########
# CHD1 #
########
# Mean ChIP-seq signal of H1 variants around the center of CHD1 peaks

# Create directory in home directory (ChIP-Seq_$date) a directory of CHD1 that must contain the files of CHD1 peaks data 
cd CHD1
awk '{print $1, $2, $3}' GSE91728_ENCFF512MCR_peaks_hg19_FDR_0.01.bed | sed 's/ /\t/g' > GSE91728_ENCFF512MCR_peaks_hg19_FDR_0.01_for_site_pro.bed
sitepro --name=H1_all_fromBam_sitepro_CHD1_FDR_0.01_peaks --pf-res=50 -b GSE91728_ENCFF512MCR_peaks_hg19_FDR_0.01_for_site_pro.bed -w ../wigFiles/H1_2.inputMerged_subtracted.wig -w ../wigFiles/H1X.inputMerged_subtracted.wig -w ../wigFiles/2_HA.inputMerged_subtracted.wig -w ../wigFiles/0_HA.inputMerged_subtracted.wig --span 3000
cd ..

########
# CTCF #
########
# Mean ChIP-seq signal of H1 variants around the center of CTCF_GSM733785 peaks

# Create directory in home directory (ChIP-Seq_$date) a directory of CTCF that must contain the files of CTCF peaks data 
cd CTCF
awk '{print $1, $2, $3}' GSM733785_hg19_wgEncodeBroadHistoneHelas3CtcfStdPk.broadPeak | sed 's/ /\t/g' > GSM733785_hg19_wgEncodeBroadHistoneHelas3CtcfStdPk.broadPeak_for_site_pro
sitepro --name=H1_all_fromBam_sitepro_CTCF_GSM733785_peaks --pf-res=50 -b CTCF_GSM733785/GSM733785_hg19_wgEncodeBroadHistoneHelas3CtcfStdPk.broadPeak_for_site_pro -w ../wigFiles/H1_2.inputMerged_subtracted.wig -w ../wigFiles/H1X.inputMerged_subtracted.wig -w ../wigFiles/2_HA.inputMerged_subtracted.wig -w ../wigFiles/0_HA.inputMerged_subtracted.wig --span 3000 --span 3000
cd ..




