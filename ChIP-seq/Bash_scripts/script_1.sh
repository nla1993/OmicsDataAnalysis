############################
# ChIP-sequencing analysis #
############################

# Required directories
# Once one is located in home directory, create directory where all processed data will be stored (for instance, in my case directory called Natalia)
# In home directory there is a directory called Softwares where many files and softwares for further analysis are required

# Some settings
date=180714
bowtieIndex=/usr/local/bin/bowtie/indexes/hg19/hg19 # for alignment with Bowtie
GenomeVersion=hg19
ChromSize=/home/andrea/Softwares/hg19.ChromSize # for genomeCoveragedBed.

#
calc () {
    bc -l <<< "$@"
}
#

mkdir ChIP-Seq_$date
cd ChIP-Seq_$date

mkdir fastq  # Raw data should be located in this directory
mkdir SamFiles
mkdir BamFiles
mkdir bdgFiles
mkdir wigFiles
mkdir BigWigFiles


# 1. Mapping with bowtie
echo "Alignment using bowtie: bowtie --best --strata -m1 --sam -n 3 --chunkmbs 200 $bowtieIndex -q $file > ../SamFiles/"$(basename "$file")$.$GenomeVersion.sam""
# --best to just show those aligments in the best alignment stratum, - m1 is used to keep only unique alignment
# -- sam specifies the output file extension (SAM file), -- n3 allows a maximum of 3 mismatches, -q specifies that
# the input file is fastq

for file in fastq/*fastq
do
bowtie --best --strata -m1 --sam -n 3 --chunkmbs 200 $bowtieIndex -q $file > SamFiles/"$(basename "$file")$.$GenomeVersion.sam" 
done

for file in SamFiles/*sam
do
mv -- "$file" "${file%.fastq$.$GenomeVersion.sam}.$GenomeVersion.sam"
done

echo "Remove fastq files"
cd fastq
rm *fastq
cd ..

# 2. Create bam file from sam files
echo "Converting sam file to bam file: samtools view -bS -o BamFiles/"$(basename "$file")$.bam" $file"
for file in SamFiles/*sam
do
samtools view -bS -o BamFiles/"$(basename "$file")$.bam" $file
done

for file in BamFiles/*bam
do
mv -- "$file" "${file%.sam$.bam}.bam"
done


echo "Remove sam files"
cd SamFiles
rm *
cd ..


# 3. Sort bam files by coordinate and index
echo "Sort BAM file by coordinate and create index using samtools:samtools sort $file -o BamFiles/"$(basename "$file")$.sorted.bam""
for file in BamFiles/*bam
do
samtools sort $file -o BamFiles/"$(basename "$file")$.sorted.bam"
done

for file in BamFiles/*sorted.bam
do
mv -- "$file" "${file%.bam$.sorted.bam}.sorted.bam"
done


echo "Remove unsorted Bam files"
cd BamFiles
rm *$GenomeVersion.bam
cd ..

#4. Create Input_merged that will be required if the decision of merging the Inputs is taken
################
# Input Merged #
################
# Input substration of the input merged
# 1. Merge of the inputs
echo "Merging inputs with: samtools merge Input_merged.sorted.bam I_fHA.hg19.sorted.bam I_HA.hg19.sorted.bam I.hg19.sorted.bam"
cd BamFiles
samtools merge Input_merged.sorted.bam I_fHA.hg19.sorted.bam I_HA.hg19.sorted.bam I.hg19.sorted.bam
cd ..


# 5. Filtering Bam Files with -F 3972 

#######################################
#      filter bam files -F 3972       #
# .bdg files from BamFiles generation #
#######################################

# 1. Filter by mapped reads
# 2. genomeCoverageBed to create .bdg
# 2.2 Normalized by million reads

echo "Filtering Bam Files with -F 3972: samtools view -b -F 3972 $file > BamFiles/"$(basename "$file")$.filtered.bam""
# Flag -F 3972 is used to remove reads that: To remove reads: unmapped, second in pair, not primary alignments, 
# reads that fail quality checks, duplicated reads and non-unique reads 


for file in BamFiles/*bam
do
samtools view -b -F 3972 $file > BamFiles/"$(basename "$file")$.filtered.bam"
done

for file in BamFiles/*.filtered.bam
do
mv -- "$file" "${file%.bam$.filtered.bam}_filtered.bam"
done

echo "Remove unfiltered Bam files"
cd BamFiles
rm *sorted.bam
cd ..

echo "1. Calculating Mapped Reads"
echo "2. Calculating Scaling Factor"
echo "3. Bam file to bdg file with genomeCoveragedBed"
# genomeCoveragedBed needs the next parameters; -ibam to specify that the input file is in Bam format
# -g to specify the chromosome sizes, -bga reports depth in BedGraph format, -scale to scale the coverage
# by a constant factor

for file in BamFiles/*bam
do
mappedReads=$(samtools view -c $file)
ScaleFactor=$(calc 1000000/$mappedReads)
genomeCoverageBed -ibam $file -g $ChromSize -bga -scale $ScaleFactor > bdgFiles/"$(basename "$file")$.bdg"
done

for file in bdgFiles/*bdg
do
mv -- "$file" "${file%.bam$.bdg}.bdg"
done

echo "Sorting bdg files by chromosome and start position" #will be required for further analysis
cd bdgFiles
for file in *bdg
do
sort -k1,1 -k2,2n $file > tmp
mv tmp $file
done
cd ..

# Until here one gets the bdg files of each of the samples you have (IPs and Inputs)

# OPTIONAL
# TIP:
# Save Bam files filtered (in your own computer) if you don't want to re-do the whole analysis again
# Once saved in own computer, remove them to save space
# Do not remove them if you want to do ChIP-seq annotations

# echo "Removing Bam files filtered"
# cd BamFiles
# rm *
# cd ..











