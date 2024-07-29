# Read arguments from command line
args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=0) {
    data_dir_name <- args[1]
    code_dir_name <- args[2]
    source(paste(code_dir_name,"universal_variables.R",sep=""))
    source(paste(code_dir_name,"functions.R",sep=""))
}

install.packages("BiocManager", dependencies = TRUE, repos = "http://cran.us.r-project.org")
BiocManager::install(version = "3.16")
using(c("DESeq2", "pheatmap", "EnhancedVolcano", "plotly","tibble","GO.db",
    "gcrma", "SCAN.UPC", "foreach","doParallel","pathview","GOplot","DT","webshot","clusterProfiler",
    "tidyverse", "showtext","plyr","knitr","ggvenn","VennDiagram", "org.Hs.eg.db"), bio=TRUE)


library(BiocManager)
BiocManager::install("DESeq2")
library(DESeq2)
library(pheatmap)
library(org.Mm.eg.db)
library(EnhancedVolcano)
library(plotly)
#library(devtools)
library(tibble)
library(GO.db)
library(gcrma)
library(SCAN.UPC)
library(foreach)
library(doParallel)
library(pathview)
library(GOplot)
library(DT)
library(webshot)
library(clusterProfiler)
library(tidyverse)
library(showtext)
library(plyr)
library(knitr)
library(ggvenn)
library(VennDiagram)
library(htmlwidgets)

# Creating deseq2 object
dds <- DESeqDataSetFromMatrix(countData = inputData,
                            colData = samples,
                            design = design)

dds <- DESeq(dds, betaPrior = betaPrior)
# Regularized log transformation for different analysis (clustering, heatmaps, etc)
rld <- rlogTransformation(dds)
pca <- plotPCA(rld, intgroup = c(colGroups))
