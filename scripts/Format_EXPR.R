# File: Format_cased_sequenced.R
# Save expr files: "EXPR__gene_tpm.csv" and "EXPR__gene_counts.csv" dim 61544 x 51
# "EXPR__isoform_tpm.csv" and "EXPR__isoform_counts.csv" dim 246624 x 51

library(data.table)
library(stringr)

# args <- commandArgs(trailingOnly = TRUE)
# input_dir <- args[1]
# output_dir <- args[2]

# Load expression list
# expr_list <- readRDS(file.path(input_dir, 'expr_list.rds'))
expr_list <- readRDS('files/expr_list.rds')

# Loop through expression list and save each as CSV
for(assay_name in names(expr_list)) {
  write.table( 
    expr_list[[assay_name]], 
    file = file.path("files", paste0('EXPR_', str_replace(assay_name, 'expr_', ''), '.csv')), 
    quote = FALSE, 
    sep = ";", 
    col.names = TRUE, 
    row.names = TRUE 
  )
}


