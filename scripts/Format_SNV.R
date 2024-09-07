# File: Format_cased_sequenced.R
# Save Format_SNV.csv (dimensions: 16706 x 8).

library(data.table)

# args <- commandArgs(trailingOnly = TRUE)
# input_dir <- args[1]
# output_dir <- args[2]

# Load SNV data from gzipped file
# snv <- as.data.frame(fread(file.path(input_dir, "SNV.txt.gz"), stringsAsFactors=FALSE, sep=";"))
snv <- as.data.frame(fread("files/SNV.txt.gz", stringsAsFactors=FALSE, sep=";"))

# Format the snv data by extracting and transforming relevant columns
snv <- cbind(snv[, c("Start", "Patient", "Hugo Symbol", "Variant Classification")],
             sapply(snv[, "Chromosome"], function(x) { paste("chr", x, sep="") }),
             sapply(snv[, "HGVS_c"], function(x) { z = unlist(strsplit(x, ">", fixed=TRUE))[1]; substr(z, nchar(z), nchar(z)) }),
             sapply(snv[, "HGVS_c"], function(x) { unlist(strsplit(x, ">", fixed=TRUE))[2] })
)

# Rename the columns
colnames(snv) <- c("Pos", "Sample", "Gene", "Effect", "Chr", "Ref", "Alt")

# Replace '-' with an empty string in Ref and Alt columns
snv$Ref <- ifelse(snv$Ref %in% "-", "", snv$Ref)
snv$Alt <- ifelse(snv$Alt %in% "-", "", snv$Alt)

# Add a mutation type column based on the length of Ref and Alt
snv <- cbind(snv,
             apply(snv[, c("Ref", "Alt")], 1, function(x) { ifelse(nchar(x[1]) != nchar(x[2]), "INDEL", "SNV") })
)

# Rename columns again to include the mutation type
colnames(snv) <- c("Pos", "Sample", "Gene", "Effect", "Chr", "Ref", "Alt", "MutType")

# Load cased sequenced data
# case <- read.csv(file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE, sep=";")
case <- read.csv("files/cased_sequenced.csv", stringsAsFactors=FALSE, sep=";")

# Filter snv to include only samples that have SNV data
snv <- snv[snv$Sample %in% case[case$snv %in% 1, ]$patient, c("Sample", "Gene", "Chr", "Pos", "Ref", "Alt", "Effect", "MutType")]

# Reformat effect column to standard mutation classifications
snv$Effect <- ifelse(snv$Effect %in% c("missense_variant", "missense_variant&splice_region_variant"), "Missense_Mutation", snv$Effect)
snv$Effect <- ifelse(snv$Effect %in% c("splice_acceptor_variant&intron_variant", "splice_acceptor_variant&splice_donor_variant&intron_variant", 
                                       "splice_donor_variant&intron_variant", "start_lost&splice_region_variant", "stop_gained&splice_region_variant"), "Splice_Site", snv$Effect)
snv$Effect <- ifelse(snv$Effect %in% "stop_gained", "Nonsense_Mutation", snv$Effect)
snv$Effect <- ifelse(snv$Effect %in% "stop_lost", "Stop_Codon_Del", snv$Effect)
snv$Effect <- ifelse(snv$Effect %in% "start_lost", "Start_Codon_Ins", snv$Effect)

# Write the formatted SNV data to a CSV file
# write.table(snv, file=file.path(output_dir, "SNV.csv"), quote=FALSE, sep=";", col.names=TRUE, row.names=FALSE)
write.table(snv, "files/SNV.csv", quote=FALSE, sep=";", col.names=TRUE, row.names=FALSE)

