# Format_downloaded_data.R

# Formats and cleans clinical and expression data.
# Outputs: 
# - "CLIN.txt" (101 x 14)
# - 'expr_list.rds' including:
# - expr_gene_tpm: 61544  x 51
# - expr_gene_counts: 61544 x 51
# - expr_isoform_tpm:  246624 x 51
# - expr_isoform_counts: 246624 x 51

library(data.table)
library(readxl) 
library(stringr)
library(tximport)
library(rhdf5)

# args <- commandArgs(trailingOnly = TRUE)
# work_dir <- args[1]
# annot_dir <- args[2]

# 1. Create CLIN.txt

# clin <- read_excel(file.path(work_dir, '1-s2.0-S0092867417311224-mmc2.xlsx'), sheet='Table S2')
clin <- read_excel("files/1-s2.0-S0092867417311224-mmc2.xlsx") # Load clinical data (73 rows, 12 columns)
mapping <- read.csv("files/filereport_read_run_PRJNA356761_tsv-DESKTOP-LJP6LNF.txt", sep = "\t") # Load mapping data


# Clean column names by replacing non-word characters with '.'
colnames(clin) <- str_replace_all(clin[2, ], '\\W', '.')
clin <- clin[-c(1:2), ] 
clin[clin == 'NA'] <- NA 

# Convert selected columns to numeric types
numcols <- c('Time.to.Death...weeks.', 'Mutation.Load', 'Neo.antigen.Load', 'Neo.peptide.Load', 'Cytolytic.Score')
clin[, numcols] <- sapply(clin[, numcols], as.numeric)

# Convert death status to logical type
clin$Dead.Alive...Dead...True. <- sapply(clin$Dead.Alive...Dead...True., as.logical)

# Extract 'PtXX' from 'sample_title' and store it as 'patient_no'
mapping$patient_no <- str_extract(mapping$sample_title, "Pt\\d+")

# Rename 'Patient' column to 'patient_no' in 'clin'
colnames(clin)[colnames(clin) == "Patient"] <- "patient_no"

# Merge 'clin' and 'mapping' dataframes by 'patient_no'
clin_merged <- merge(clin, mapping, by = "patient_no")

# Rename 'run_accession' to 'patient' in the merged dataframe
colnames(clin_merged)[colnames(clin_merged) == "run_accession"] <- "Patient"

# Save the final clinical data to a file
write.table(clin_merged, "files/CLIN.txt", col.names = TRUE, sep = '\t')

# 2. Create expr_list.rds
source('https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/process_kallisto_output.R')

# Load annotation data
load("files/Gencode.v40.annotation.RData")

# Function to process kallisto output for expression data
process_kallisto_output <- function(work_dir, tx2gene) {
  samples <- list.dirs(file.path(work_dir, 'rnaseq'), full.names = FALSE, recursive = FALSE)
  files <- file.path(work_dir, 'rnaseq', samples, "abundance.h5")
  names(files) <- samples
  
  # Import transcript and gene-level expression data
  expr_tx <- tximport(files, type = "kallisto", txOut = TRUE, ignoreAfterBar = TRUE)
  expr_gene <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
  
  # Create list of processed expression data
  expr_list <- list(
    expr_gene_tpm = log2(expr_gene$abundance + 0.001),
    expr_gene_counts = log2(expr_gene$counts + 1),
    expr_isoform_tpm = log2(expr_tx$abundance + 0.001),
    expr_isoform_counts = log2(expr_tx$counts + 1)
  )
  
  return(expr_list)
}

# Define working directory and process kallisto output
work_dir <- "Gide_kallisto_v0.46.1_GRCh38.40/files"
expr_list <- process_kallisto_output(work_dir, tx2gene)

# Save expression list to an RDS file
saveRDS(expr_list, file = "files/expr_list.rds")

# 3. Create SNV.txt.gz
#snv <- read_excel(file.path(work_dir, '1-s2.0-S0092867417311224-mmc3.xlsx'), sheet='Table S3')
snv <- read_excel("files/1-s2.0-S0092867417311224-mmc3.xlsx")

colnames(snv) <- snv[3, ] 
snv <- snv[-c(1:3), ]     
numcols <- c('Start', 'End', 'Tcov', 'Tac', 'Taf')

# Convert columns to numeric
snv[, numcols] <- sapply(snv[, numcols], as.numeric)

# Extract 'run_accession' from the 'mapping' dataframe
snv$run_accession <- mapping$run_accession[match(snv$Patient, mapping$patient_no)]  # Match 'Patient' in SNV to 'patient_no' in mapping

#: Replace 'PtXX' in snv with 'run_accession'
snv$Patient <- snv$run_accession

# Write SNV data to a gzipped text file
#gz <- gzfile(file.path(work_dir, 'SNV.txt.gz'), "w")
gz <- gzfile("files/SNV.txt.gz", "w")  
write.table(snv, file = gz, quote = FALSE, sep = ";", col.names = TRUE, row.names = FALSE)
close(gz)

