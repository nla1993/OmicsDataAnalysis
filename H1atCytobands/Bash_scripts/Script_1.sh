###########################
# H1 studies at cytobands #
###########################

# Required directories
# Once one is located in home directory, create directory where all processed data will be stored (for instance, in my case directory called Natalia)
# In home directory there is a directory called Softwares where many files and softwares for further analysis are required

# Some settings
date=180714
# Some settings
mkdir H1_Cytobands_$date

HomeDirectory=/home/andrea/Natalia/H1_Cytobands_$date
hgGcPercentDirectory=/home/andrea/Softwares/hgGcPercent
GenomeVersion=hg19

# In H1_Cytobands_$date directory, bedfiles corresponding to the coordinates of each of the cytoband groups should be located.

cd $HomeDirectory
mkdir GC_content
##############
# GC content #
##############
# GC content in each band of gpos25 set of bands
for file in *bed
do
../../Softwares/hgGcPercent -bedRegionIn=$file -bedRegionOut=GC_content/"$(basename "$file")_GCcontent.bed" -noLoad $GenomeVersion ../../Softwares/hg19.2bit
done

# to do the mean value
cd GC_content
for file in *GC_content.bed
do
awk '{print $0, $4/($3 -$2)}' $file | sed 's/ /\t/g' > tt
mv tt "$(basename "$file")"
done

# now I need to bind this with the names of each band which are in the original file (cytoBand_gpos25.bed)
paste -d' ' ../cytoBand_gpos25.bed cytoBand_gpos25_GCcontent.bed | sed 's/ /\t/g' | awk '{print $4, $11, $4$1}' | sed 's/ /\t/g' > tt
mv tt cytoBand_gpos25_GCcontent.bed # Repeat this proccess for all cytoband groups.
cd ..


########################
# H3K4me3 at cytobands #
########################
mkdir PTMs
cd PTMs

#sort bed file of the peaks
sort -k1,1 -k2,2n ENCFF387SID.bed > tmp
mv tmp ENCFF387SID.bed


for file in ../cytoBand*
do
intersectBed -a $file -b ENCFF387SID.bed -c > "H3K4me3$(basename "$file")"
done

# intersectBed is used to get the overlaping between two genomic features. -c reports, for each entry in A, report the number of hits in B

########################
# H3K27me3 at cytobands #
########################

#sort bed file of the peaks
sort -k1,1 -k2,2n ENCFF252BLX.bed > tmp
mv tmp ENCFF252BLX.bed

for file in ../cytoBand*
do
intersectBed -a $file -b ENCFF252BLX.bed -c > "H3K27me3$(basename "$file")"
done

########################
# H3K9me3 at cytobands #
########################

#sort bed file of the peaks
sort -k1,1 -k2,2n GSM1003480_hg19_wgEncodeBroadHistoneHelas3H3k09me3Pk.broadPeak > tmp
mv tmp GSM1003480_hg19_wgEncodeBroadHistoneHelas3H3k09me3Pk.broadPeak

for file in ../cytoBand*
do
intersectBed -a $file -b GSM1003480_hg19_wgEncodeBroadHistoneHelas3H3k09me3Pk.broadPeak -c > "H3K9me3$(basename "$file")"
done

cd ..

#########################################
# Abundance of H1 variants at cytobands #
#########################################
mkdir H1abundance

for file in cytoBand*
sort -k1,1 -k2,2n $file > tmp
mv tmp "$(basename "$file")"
done

#########################################
# Abundance of H1 variants at cytobands #
#                gpos25                 #     
#########################################

for file in ../ChIP-Seq_$date/bdgFiles/*bdg
do
bedtools map -a cytoBand_gpos25.bed -b $file -c 4 -o mean > H1abundance/"$(basename "$file")_gpos25_overlap.bed"
done

for file in H1abundance/*gpos25*
do
awk '{print $0, $4$1}' $file | sed 's/ /\t/g' > tt
mv tt "$(basename "$file")"
done


#########################################
# Abundance of H1 variants at cytobands #
#                gpos50                 #     
#########################################

for file in ../ChIP-Seq_$date/bdgFiles/*bdg
do
bedtools map -a cytoBand_gpos50.bed -b $file -c 4 -o mean > H1abundance/"$(basename "$file")_gpos50_overlap.bed"
done

for file in H1abundance/*gpos50*
do
awk '{print $0, $4$1}' $file | sed 's/ /\t/g' > tt
mv tt "$(basename "$file")"
done

#########################################
# Abundance of H1 variants at cytobands #
#                gpos75                 #     
#########################################

for file in ../ChIP-Seq_$date/bdgFiles/*bdg
do
bedtools map -a cytoBand_gpos75.bed -b $file -c 4 -o mean > H1abundance/"$(basename "$file")_gpos75_overlap.bed"
done

for file in H1abundance/*gpos75*
do
awk '{print $0, $4$1}' $file | sed 's/ /\t/g' > tt
mv tt "$(basename "$file")"
done

#########################################
# Abundance of H1 variants at cytobands #
#                gpos100                 #     
#########################################

for file in ../ChIP-Seq_$date/bdgFiles/*bdg
do
bedtools map -a cytoBand_gpos100.bed -b $file -c 4 -o mean > H1abundance/"$(basename "$file")_gpos100_overlap.bed"
done

for file in H1abundance/*gpos100*
do
awk '{print $0, $4$1}' $file | sed 's/ /\t/g' > tt
mv tt "$(basename "$file")"
done

# PLOTS DONE WITH R
# look in R script called '1.H1Cytobands.R'
