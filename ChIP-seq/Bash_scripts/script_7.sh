############################
# ChIP-sequencing analysis #
############################

# Further analyses

# ChIP regions

# Some settings
mkdir bedFiles_from_BamFiles
mkdir BedFiles
cd BedFiles
mkdir annotations
cd ..

HomeDirectory=/home/andrea/Natalia/ChIP-Seq_$date
SICERDirectory=/home/andrea/Softwares/SICER_V1.1/SICER

date=180714
GenomeVersion=hg19
RedundancyThreshold=1
WindowSize=200
FragmentSize=150
EffectiveGenomeFraction=0.75
GapSize=200
FDR=0.01

# 1. Bam to Bed files
echo 'Convert Bam Files to Bed Files'
for file in BamFiles/*bam
do
bamToBed -i $file > bedFiles_from_BamFiles/"$(basename "$file").bed"
done

for file in bedFiles_from_BamFiles/*bed
do
mv -- "$file" "${file%.bam.bed}.bed"
done

# 2. Run SICER.sh with specific parameters
# DO INDIVIDUALLY
# EXAMPLE: H1_2.hg19.sorted_filtered.bed
echo "Generation of peak files with SICER"
cd $SICERDirectory

./SICER.sh $HomeDirectory/bedFiles_from_BamFiles/ H1_2.hg19.sorted_filtered.bed Input_merged.sorted_filtered.bed $HomeDirectory/BedFiles/ $GenomeVersion $RedundancyThreshold $WindowSize $FragmentSize $EffectiveGenomeFraction $GapSize $FDR

cd $HomeDirectory/BedFiles
rm !(*island.bed)
cd $HomeDirectory

# Annotations with Homer

# Not in the cluster. Run in own computer where HOMER is installed

##################
# Annotate peaks #
##################

# Add strand information required for Homer        H1X.hg19.sorted_filtered-W200-G200-FDR0.01-island.bed
echo "Add strand information required for Homer"
for file in BedFiles/*bed
do
awk '{print "peak_"i++, $1,$2,$3, "+"}' $file | sed 's/ /\t/g' > BedFiles/annotations/"$(basename "$file").annot"
done

for file in BedFiles/annotations/*annot
do
mv -- "$file" "${file%.hg19.sorted_filtered-W200-G200-FDR0.01-island.bed.annot}._peaks.bed"
done

echo "Run annotations"
for file in BedFiles/annotations/*peaks.bed
do
annotatePeaks.pl $file $GenomeVersion > BedFiles/annotations/"$(basename "$file").txt"
done

# Plots are done in Excel              



