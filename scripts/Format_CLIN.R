# File: Format_cased_sequenced.R
# Goal: Save cased_sequenced.csv (dimensions: 101  30)

library(stringr)
library(tibble)

# Source utility functions
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/Get_Response.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/format_clin_data.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_tissue.R")

# Read clinical data
clin_original <- read.csv("files/CLIN.txt", stringsAsFactors=FALSE, sep="\t") # dim 101 x 14

# Select relevant columns for clinical data
selected_cols <- c("Patient", "Response", "Dead.Alive...Dead...True.", "Time.to.Death...weeks.", "Subtype", "M.Stage")
clin <- cbind(clin_original[, selected_cols], "Melanoma", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

# Rename columns
colnames(clin) <- c("patient", "recist", "os", "t.os", "histo", "stage", "primary", "drug_type", 
                    "pfs", "t.pfs", "sex", "age", "dna", "dna_info", "rna", "rna_info", "response.other.info", "response")

# Update stage information
clin$stage <- ifelse(str_detect(clin$stage, 'M0'), "III", 
                     ifelse(str_detect(clin$stage, 'M1'), "IV", NA))

# Replace "NE" with NA in recist column
clin$recist[clin$recist == "NE"] <- NA 

# Generate response data
clin$response <- Get_Response(data = clin)

# Convert "TRUE" to 1 and everything else to 0 in the os column
clin$os <- ifelse(clin$os == "TRUE", 1, 0)

# Convert weeks to months for t.os
clin$t.os <- clin$t.os / 4 

# Read case sequenced data
case <- read.csv("files/cased_sequenced.csv", stringsAsFactors=FALSE, sep=";") # dim 101 x 4

# RNA-related columns
clin$rna <- ifelse(clin$patient %in% case[case$expr == 1, ]$patient, "rnaseq", NA)
clin$rna_info <- ifelse(clin$patient %in% case[case$expr == 1, ]$patient, "tpm", NA)

# DNA-related columns
clin$dna <- ifelse(clin$patient %in% case[case$snv == 1, ]$patient, "wes", NA)
clin$dna_info <- ifelse(clin$patient %in% case[case$snv == 1, ]$patient, "snv", NA)

# Reorder columns
clin <- clin[, c(
  "patient", "sex", "age", "primary", "histo", "stage", 
  "response.other.info", "recist", "response", "drug_type", "dna", "dna_info", "rna", "rna_info", "t.pfs", 
  "pfs", "t.os", "os"
)]

# Use the format_clin_data function for further formatting
clin <- format_clin_data(clin_original, "Patient", selected_cols, clin)

# Read curation tissue data
path <- "https://raw.githubusercontent.com/BHKLAB-DataProcessing/ICB_Common/main/data/curation_tissue.csv"
annotation_tissue <- read.csv(path)

# Annotate 'clin' using tissue annotation
clin <- annotate_tissue(clin = clin, study = 'Riaz', annotation_tissue = annotation_tissue, check_histo = TRUE)

# Add treatmentid based on Cohort
clin <- add_column(clin, treatmentid = ifelse(clin_original$Cohort == "NIV3-PROG", "Combo", "PD-1/PD-L1"), .after = 'tissueid')

# Replace empty strings with NA
clin[clin == "-"] <- NA

# Save the formatted clinical data to a CSV file
write.table(clin, "files/CLIN.csv", quote=FALSE, sep=";", col.names=TRUE, row.names=FALSE)
