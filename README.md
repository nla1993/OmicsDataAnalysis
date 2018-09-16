# OmicsDataAnalysis
Scripts used for 1) ChIP-sequencing Analysis and 2) Cytogenetic Band Classification Analysis

Scripts: ChIP-sequencing analysis

BASH Scripts

Required software: Bowtie, SamTools, BedTools, bedgraph_to_wig.pl, wigCorrelate, MACS2, wigToBigWig, CEAS and SICER.sh
Required files: hg19.ChromSize, RefSeq_Genes.hg19.strand.bed and hg19_ceas.refGene

Script 1: 
-	Mapping with Bowtie
-	Convert Sam to Bam files
-	Sort Bam fles by coordinate and index
-	Merging of the input files
-	Filtering Bam files
-	Convert Bam to Bedgraph (bdg) files

Script 2:
-	Convert bdg to wig files (before input merge subtraction)
-	Wig correlation

Script 3: 
-	Input merged subtraction
-	Convert input merged subtracted bdg files to wig files
-	Convert Wig To BigWig files

Script 4:
-	Wig correlation (after input merge subtraction)

Script 5:
-	 Create gene expression groups for HeLa RNA-seq data

Script 6; Downstream analyses:
-	 ChIP-seq average signal profile by expression groups at TSS and gene body
-	H1 abundace at LADs, centromeres, FAIRE-seq peaks and enhancers
-	Occupancy of H1 variants at specific histone modifications
-	 ChIP-seq average signal profile at RNA Pol II binding sites, CHD1 peaks and CTCF peaks.

Script 7; Peak calling analysis. Not shown in the project:
-	Convert Bam to Bed files
-	Run SICER.sh to create peak files
-	Annotate peaks with Homer

R scripts have been used to create the visualization of the analysis results. In each of the Bash scripts is specified which R script correspond to each section.


Scripts: Cytogenetic Band Classification

BASH Scripts

Required software: hgGcPercent and BedTools

Script 1: 
-	GC content calculation
-	Abundance of PTMs at cytogenetic bands
-	Abundance of H1 variants at cytogenetic bands
-	Representation of the results are done in R. It can be found in R script called ‘1.H1Cytobands.R’

Script 2; Cytogenetic band classification:
-	Intersect cytogenetic band BED file with with all of the Peak files corresponding to the factors of the study
-	Scale by cytogenetic band length

R scripts

Script 1: Reading of the generated data scaled by cytogenetic band length and generation of unique table containing all the data

Script 2: Boxplot showing the length of each of the cytogenetic band groups

Script 3: Reading of the files containing H1 abundance at cytogenetic bands and generation of unique table containing all the data and scale of the data.

Script 4: Generation of the heatmap

Script 5: Defining new 8 groups of cytogenetic bands

Script 6: Representation in pie chart plot of the distribution of cytogenetic bands from groups gpos 25-100 in the new defined groups

Script 7: Boxplots showing H1 variants abundance at the new defined groups

Script 8: Scatter plots showing correlation between the different variables of study
