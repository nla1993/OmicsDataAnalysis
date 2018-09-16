############################
# ChIP-sequencing analysis #
############################

# Wig correletation (after input subtraction)
# go to directory wigfiles
cd wigFiles
../../../Softwares/wigCorrelate # NAMES OF ALL THE WIG FILES SEPARATED BY AN SPACE
cd ..
# WIG CORRELATION RESULTS #
# Wig correlation results are printed directly in the console.

# PLOTS DONE WITH R
# look in R script called '1.WigCorrelationPlots.R'