if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsubread")
library(Rsubread)


args <- commandArgs( trailingOnly = TRUE )
bamFile <- args[1]
gtfFile <- args[2]
outFilePref <- args[3]

outStatsFilePath  <- paste(outFilePref, '.stat',  sep = '');
outCountsFilePath <- paste(outFilePref, '.count', sep = '');

fCountsList = featureCounts(bamFile, annot.ext=gtfFile, isGTFAnnotationFile=TRUE, GTF.attrType="Name", GTF.featureType="miRNA")

write.table(fCountsList$stat, outStatsFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

featureCounts = cbind(fCountsList$annotation[,1], fCountsList$counts)
write.table(featureCounts, outCountsFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
