########################
# Cytoband aggregation #
########################

# Required directories
# Once one is located in home directory, create directory where all processed data will be stored (for instance, in my case directory called Natalia)
# In home directory there is a directory called Softwares where many files and softwares for further analysis are required
# In H1_Cytobands_$date directory, bedfiles corresponding to the coordinates of each of the cytoband groups should be located.

########################################
# Abundance of 12 factors at cytobands #
########################################

HomeDirectory=/home/andrea/Natalia/H1_Cytobands_$date

# Abundance of H3K4me3, H3K9me3 and H3K27me3 and abundance of H1 variants at cytobands are already done.
# They can be found in directory H1atCytobands > 1.H1Cytobands.sh

cd PTMs

####################
# P300_ENCFF001VIZ #
####################
cd PTMs
#sort bed file of the peaks
sort -k1,1 -k2,2n ENCFF001VIZ.bed > tmp
mv tmp ENCFF001VIZ.bed

for file in ../cytoBand*
do
intersectBed -a $file -b ENCFF001VIZ.bed -c > "P300_ENCFF001VIZ$(basename "$file")"
done

####################
# EZH2_GSM1003520  #
####################

#sort bed file of the peaks
sort -k1,1 -k2,2n GSM1003520_hg19_wgEncodeBroadHistoneHelas3Ezh239875Pk.broadPeak > tmp
mv tmp GSM1003520_hg19_wgEncodeBroadHistoneHelas3Ezh239875Pk.broadPeak

for file in ../cytoBand*
do
intersectBed -a $file -b GSM1003520_hg19_wgEncodeBroadHistoneHelas3Ezh239875Pk.broadPeak -c > "EZH2_GSM1003520$(basename "$file")"
done

#############
# RNAPOL II #
#############

#sort bed file of the peaks
sort -k1,1 -k2,2n RNA_pol_II_HeLa.narrowPeak > tmp
mv tmp RNA_pol_II_HeLa.narrowPeak

for file in ../cytoBand*
do
intersectBed -a $file -b RNA_pol_II_HeLa.narrowPeak -c > "RNAPOLII$(basename "$file")"
done


#############
#  H3K36me3 #
#############

#sort bed file of the peaks
sort -k1,1 -k2,2n GSM733711_hg19_wgEncodeBroadHistoneHelas3H3k36me3StdPk.broadPeak > tmp
mv tmp GSM733711_hg19_wgEncodeBroadHistoneHelas3H3k36me3StdPk.broadPeak

for file in ../cytoBand*
do
intersectBed -a $file -b GSM733711_hg19_wgEncodeBroadHistoneHelas3H3k36me3StdPk.broadPeak -c > "H3K36me3_GSM733711$(basename "$file")"
done

#################
# CHD1 filtered #
#################

#sort bed file of the peaks
sort -k1,1 -k2,2n GSE91728_ENCFF512MCR_peaks_hg19_FDR_0.01.bed > tmp
mv tmp GSE91728_ENCFF512MCR_peaks_hg19_FDR_0.01.bed

for file in ../cytoBand*
do
intersectBed -a $file -b GSE91728_ENCFF512MCR_peaks_hg19_FDR_0.01.bed -c > "CHD1_FDR_0.01$(basename "$file")"
done


#####################
# H3K4me1_GSM798322 #
#####################

#sort bed file of the peaks
sort -k1,1 -k2,2n GSM798322_hg19_wgEncodeBroadHistoneHelas3H3k04me1StdPk.broadPeak > tmp
mv tmp GSM798322_hg19_wgEncodeBroadHistoneHelas3H3k04me1StdPk.broadPeak

for file in ../cytoBand*
do
intersectBed -a $file -b GSM798322_hg19_wgEncodeBroadHistoneHelas3H3k04me1StdPk.broadPeak -c > "H3K4me1_GSM798322$(basename "$file")"
done

#####################
# H3K9ac_GSM733756 #
#####################

#sort bed file of the peaks
sort -k1,1 -k2,2n GSM733756_hg19_wgEncodeBroadHistoneHelas3H3k9acStdPk.broadPeak > tmp
mv tmp GSM733756_hg19_wgEncodeBroadHistoneHelas3H3k9acStdPk.broadPeak

for file in ../cytoBand*
do
intersectBed -a $file -b GSM733756_hg19_wgEncodeBroadHistoneHelas3H3k9acStdPk.broadPeak -c > "H3K9ac_GSM733756$(basename "$file")"
done

########################
# H3K27ac at cytobands #
########################

#sort bed file of the peaks
sort -k1,1 -k2,2n GSM733684_hg19_wgEncodeBroadHistoneHelas3H3k27acStdPk.broadPeak > tmp
mv tmp GSM733684_hg19_wgEncodeBroadHistoneHelas3H3k27acStdPk.broadPeak

for file in ../cytoBand*
do
intersectBed -a $file -b GSM733684_hg19_wgEncodeBroadHistoneHelas3H3k27acStdPk.broadPeak -c > "H3K27ac$(basename "$file")"
done

#####################
# CTCF at cytobands #
#####################

#sort bed file of the peaks
cd ../CTCF_GSM733785
sort -k1,1 -k2,2n GSM733785_hg19_wgEncodeBroadHistoneHelas3CtcfStdPk.broadPeak > tmp
mv tmp GSM733785_hg19_wgEncodeBroadHistoneHelas3CtcfStdPk.broadPeak

for file in ../cytoBand*
do
intersectBed -a $file -b GSM733785_hg19_wgEncodeBroadHistoneHelas3CtcfStdPk.broadPeak -c > "CTCF$(basename "$file")"
done

cd ..

###############################
# Division by cytoband length #
###############################
mkdir normalized_by_cytoband_length

for file in PTMs/*bed
do
awk '{print $0, $6/($3-$2)}' $file | sed 's/ /\t/g' > normalized_by_cytoband_length/"$(basename "$file")_norm_cytob_length"
done

