############################
# ChIP-sequencing analysis #
############################


# 7. 
#############################
# Input Merged Subtraction #
#############################

# Input subtraction (IP-IN) This is done by MACS2 but remembering that in these analyses we are NOT producing the .bdg files from MACS but
# with beedtools


for file in bdgFiles/*bdg
do
macs2 bdgcmp -t $file -c bdgFiles/Input_merged.sorted_filtered.bdg -m subtract -o bdgFiles/"$(basename "$file")$.inputMerged_subtracted.bdg"
done
#
for file in bdgFiles/*_subtracted.bdg
do
mv -- "$file" "${file%.hg19.sorted_filtered.bdg$.inputMerged_subtracted.bdg}.inputMerged_subtracted.bdg"
done
#
rm bdgFiles/I*inputMerged_subtracted.bdg #Input - Input
#


# 8.
###############################
# bdg to wig (wig correlate)  #
# with InputMerged_subtracted #
#      need for Ceas          #
###############################

# Convert .bdg to wig

echo "Convert bdg files to wig files with inputsubtracted files which will be used for further analysis"
for file in bdgFiles/*inputMerged_subtracted.bdg
do
perl ../../Softwares/bedgraph_to_wig.pl --bedgraph $file --wig wigFiles/"$(basename "$file")$.wig" --step 50
done

for file in wigFiles/*wig
do
mv -- "$file" "${file%.bdg$.wig}.wig"
done


# 9. Wig to BigWig
#BigWig files are used to visualize wig files (signal files) in the genome browser

echo "Wig to BigWig Files"
for file in wigFiles/*wig
do
../../Softwares/wigToBigWig $file $ChromSize BigWigFiles/"$(basename "$file").bigwig"
done


for file in BigWigFiles/*bigwig
do
mv -- "$file" "${file%.wig.bigwig}.bigwig"
done




