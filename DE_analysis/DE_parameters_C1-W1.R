#r script to set parameters for DE analysis
count_matrix_file = ""
samples_file = ""
output_directory = ""
save.tables = T

#set samples to include in the analysis
keep_samples = c()


#set thresholds for edgeR
lFC_threshold_edgeR = 1
FDR_threshold_edgeR = 0.05

#set thresholds for DESeq2
lFC_threshold_deseq = 1
FDR_threshold_deseq = 0.05

#set tests
#indicate each test as a separate binary vector. The reference condition has to be indicated in the second position.
test <- c("", "")

test_name <- paste0(test[1],"-",test[2])

#annotation (tabular file with gene names and/ordescriptions)
annotation_file <- ""

#set parameters for GO enrichment
keyType = "TAIR"
ont = "BP"     # Options are "BP", "MF", "CC", and "ALL"
pAdjustMethod = "BH"    # Benjamini-Hochberg adjustment method
pvalueCutoff = 0.05
qvalueCutoff = 0.2

#parameters for GO simplification
simCutoff = 0.7
simBy = "p.adjust"
simSelect_fun = min
simMeasure = "Wang"
simSemData = NULL