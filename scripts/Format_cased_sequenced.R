# File: Format_cased_sequenced.R
# Goal: Save cased_sequenced.csv (dimensions: 101 x 4).

library(data.table)

# args <- commandArgs(trailingOnly = TRUE)
# input_dir <- args[1]
# output_dir <- args[2]

# Load expression list from RDS file
# expr_list <- readRDS(file.path(input_dir, 'expr_list.rds'))
expr_list <- readRDS('files/expr_list.rds')
expr <- expr_list[['expr_gene_tpm']]
expr_patient <- sort(colnames(expr))  # Sort patient names from expression data

# Load clinical data
# clin <- read.csv(file.path(input_dir, "CLIN.txt"), stringsAsFactors=FALSE, sep="\t")
clin <- read.csv("files/CLIN.txt", stringsAsFactors=FALSE, sep="\t")
patient <- sort(unique(clin$Patient))  # Sort unique patients from clinical data

# Load SNV data
snv <- read.table("files/SNV.txt.gz", header=TRUE, stringsAsFactors=FALSE, sep=";")
snv_patient <- sort(unique(snv$Patient))  # Sort unique patients from SNV data

# Initialize dataframe to store patient sequencing data
case <- as.data.frame(cbind(patient, rep(0, length(patient)), rep(0, length(patient)), rep(0, length(patient))))
colnames(case) <- c("patient", "snv", "cna", "expr")
rownames(case) <- patient

# Convert snv, cna, and expr columns to numeric
case$snv <- as.numeric(as.character(case$snv))
case$cna <- as.numeric(as.character(case$cna))
case$expr <- as.numeric(as.character(case$expr))

# Update sequencing data for each patient based on SNV and expression data
for(i in 1:nrow(case)) {
  if(rownames(case)[i] %in% snv_patient) {
    case$snv[i] = 1
  }
  if(rownames(case)[i] %in% expr_patient) {
    case$expr[i] = 1
  }
}

# Write the result to a CSV file
#write.table(case, file=file.path(output_dir, "cased_sequenced.csv"), quote=FALSE, sep=";", col.names=TRUE, row.names=FALSE)
write.table(case, "files/cased_sequenced.csv", quote=FALSE, sep=";", col.names=TRUE, row.names=FALSE)
