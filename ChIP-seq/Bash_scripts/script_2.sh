############################
# ChIP-sequencing analysis #
############################
#ALWAYS RUN FROM CHIP-SEQ DIRECTORY

# 6. bdg to wig to do wig correlate previously to decide in merging input files or not

echo "bdg files to wig files. This is to do wig correlation and decide if using input merge file"
for file in bdgFiles/*bdg
do
perl ../../Softwares/bedgraph_to_wig.pl --bedgraph $file --wig wigFiles/"$(basename "$file")$.wig" --step 50
done

for file in wigFiles/*wig
do
mv -- "$file" "${file%.bdg$.wig}.wig"
done

# Wig correletation (wig files separated by space)
# go to directory wigfiles
cd wigFiles
../../../Softwares/wigCorrelate # NAMES OF ALL THE WIG FILES SEPARATED BY A SPACE
cd ..
# WIG CORRELATION RESULTS #
# Wig correlation results are printed directly in the console.
# If the correlation values of the Inputs are the most similar to each other, we decide to merge them and have an unique Input file

# Once you have wig correlation results, remove these wig files not to infere further analyses where these wig files are not required
echo "Removing wig files after obtaining wig correlation results"
cd wigFiles
rm *
cd..

# PLOTS DONE WITH R
# look in R script called '1.WigCorrelationPlots.R'