############################
# ChIP-sequencing analysis #
############################


# Uncompress fastq.gz files
cd fastq.qz
# CHANGE THIS TO HELA #gunzip 2018-05-04/*fastq.gz
# CHANGE THIS TO HELA #gunzip 2018-04-25/*fastq.gz


# Merge fastq files (from two sequencing runs)
# CHECK CHECK CHECK CHECK CHECK CHECK CHECK CHECK CHECK CHECK CHECK CHECK CHECK
cat I_fHA_22035_ATTCTC.fastq I_fHA_22035_ATTCTCT.fastq > ../fastq/I_fHA.fastq
rm I_fHA_22035_ATTCTC.fastq
rm I_fHA_22035_ATTCTCT.fastq

cat 0_fHA_22036_TGACTT.fastq 0_fHA_22036_TGACTTG.fastq > ../fastq/0_fHA.fastq
rm 0_fHA_22036_TGACTT.fastq
rm 0_fHA_22036_TGACTTG.fastq

cat 2_fHA_22037_TAAGAT.fastq 2_fHA_22037_TAAGATG.fastq > ../fastq/2_fHA.fastq
rm 2_fHA_22037_TAAGAT.fastq
rm 2_fHA_22037_TAAGATG.fastq 

cat 4_fHA_22038_ACTTAC.fastq 4_fHA_22038_ACTTACG.fastq > ../fastq/4_fHA.fastq
rm 4_fHA_22038_ACTTAC.fastq
rm 4_fHA_22038_ACTTACG.fastq

cat X_fHA_22039_AATCCG.fastq X_fHA_22039_AATCCGT.fastq > ../fastq/X_fHA.fastq
rm X_fHA_22039_AATCCG.fastq
rm X_fHA_22039_AATCCGT.fastq

cat I_HA_22040_GCAACG.fastq I_HA_22040_GCAACGC.fastq > ../fastq/I_HA.fastq
rm I_HA_22040_GCAACG.fastq
rm I_HA_22040_GCAACGC.fastq

cat 0_HA_22041_GCCTGA.fastq 0_HA_22041_GCCTGAA.fastq > ../fastq/0_HA.fastq
rm 0_HA_22041_GCCTGA.fastq
rm 0_HA_22041_GCCTGAA.fastq

cat 2_HA_22042_AGGTCG.fastq 2_HA_22042_AGGTCGG.fastq > ../fastq/2_HA.fastq
rm 2_HA_22042_AGGTCG.fastq
rm 2_HA_22042_AGGTCGG.fastq

cat I_22032_CGTTGG.fastq I_22032_CGTTGGT.fastq > ../fastq/I.fastq
rm I_22032_CGTTGG.fastq
rm I_22032_CGTTGGT.fastq

cat H1_2_22033_TTCTGG.fastq H1_2_22033_TTCTGGT.fastq > ../fastq/H1_2.fastq
rm H1_2_22033_TTCTGG.fastq
rm H1_2_22033_TTCTGGT.fastq

cat H1X_1_22034_AACCGG.fastq H1X_1_22034_AACCGGT.fastq > ../fastq/H1X.fastq 
rm H1X_1_22034_AACCGG.fastq
rm H1X_1_22034_AACCGGT.fastq

cd .. # to be in ChIP-Seq_$date directory