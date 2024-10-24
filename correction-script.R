# Script for multiple hypothesis correction of mammalian PhyloP scores.
# Example call: #Rscript ./phyloP_BH_correction.R $MAMMAL_NAMA_PARAM$
# For every mammal

args <- commandArgs(trailingOnly = TRUE)
mammal <- args[1]

obsPhyloP <- read.table(paste0("/orange/kgraim/data/phyloP/Input_files/cat_runs/",mammal,".phyloP"))
nonParaPhyloP <- read.table(paste0("/orange/kgraim/data/phyloP/Input_files/cat_shufs/", mammal,".phyloP"))
colnames(obsPhyloP) <- c("chr",	"start",	"end",	"name",	"null_scale",	"alt_scale",	"alt_subscale",	"lnlratio",	"pval")
colnames(nonParaPhyloP) <- c("chr",	"start",	"end",	"name",	"null_scale",	"alt_scale",	"alt_subscale",	"lnlratio",	"pval")

empirical.pval <- function(x, dist) { sum(x <= dist)/length(dist) }
nonParaPval <- sapply(obsPhyloP$lnlratio, empirical.pval, nonParaPhyloP$lnlratio)
#nonParaPval <-sapply(obsPhyloP$V8, empirical.pval, nonParaPhyloP$V8)
nonParaFDR <- p.adjust(nonParaPval, method="BH")

final_df <- cbind(obsPhyloP,nonParaFDR)
write.table(final_df, paste0("/orange/kgraim/data/phyloP/Input_files/corrected_runs/",mammal,"_corrected.phyloP"), sep = "\t", col.names = TRUE, row.names = FALSE)

#Also you need to join the nonParaFDR to the obsPhyloP table as an additional column and before you output the final result wiht a BH corrected file.



