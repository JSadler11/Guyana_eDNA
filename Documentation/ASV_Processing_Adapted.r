# =============================================================================
# SECTION 1: LIBRARIES
# Load all packages required for the full analysis pipeline.
# Core wrangling: tidyverse (dplyr, ggplot2, tidyr), readr, data.table, scales.
# Community ecology: vegan (diversity indices, NMDS, PERMANOVA, CCA),
#   cluster, pairwiseAdonis (post-hoc pairwise PERMANOVA).
# Modelling: mgcv (Generalised Additive Models for alpha-diversity responses).
# Visualisation: plotly (interactive ggplot widgets), corrmorant
#   (pairwise covariate correlation matrices).
# Additional packages loaded mid-script: microDecon (blank-based
#   decontamination), betapart (beta-diversity partitioning), ggpubr
#   (multi-panel figure assembly), gam/mixOmics (GAMs and sparse PLS).
# =============================================================================

library(tidyverse)
library(RColorBrewer)
library(biomformat)
library(readr)
library(data.table)
library(vegan)
library(plotly)
library(cluster)
library(pairwiseAdonis)
library(corrmorant)
library(scales)
library(mgcv)
library(ggpubr)
library(effectsize)
library(see)
library(lme4)
library(lmerTest)
library(emmeans)
library(lsmeans)
library(readr)
library(gratia)
setwd("/Users/jacksadler/Desktop/project_thesis/taxonomy/12SV5_BLAST/Full_Outputs")

#________________________________________________________
# SECTION 2: ASV SEQUENCE LENGTH FILTERING → BLAST INPUT
# Reads repseqs FASTA, filters to valid 12SV5 amplicon lengths,
# exports length-filtered FASTA for BLAST on the cluster.
# NOTE: This stream is independent of Stream 2 below.
# The two streams converge downstream when BLAST output is
# merged back with the abundance table.
#________________________________________________________

# -----------Read and parse FASTA files
fasta_path <- "qza_conversion_full/repseqs_12SV5_full.fasta"
lines <- readLines(fasta_path)

#------------Convert fasta to txt--------------------------------

header_idx <- grep("^>", lines)
headers <- sub("^>", "", lines[header_idx])  # remove the ">" prefix

sequences <- sapply(seq_along(header_idx), function(i) {
  start <- header_idx[i] + 1
  end <- ifelse(i < length(header_idx), header_idx[i + 1] - 1, length(lines))
  paste(lines[start:end], collapse = "")
})

# ---- Build data frame and write to .tsv ----
df <- data.frame(header = headers, sequence = sequences, stringsAsFactors = FALSE)
write.table(df, "repseqs_12SV5.txt", sep = "\t", row.names = FALSE, quote = FALSE)

repseqs_12SV5 <- read_table("repseqs_12SV5.txt")
colnames(repseqs_12SV5)[colnames(repseqs_12SV5) == "header"] <- "fasta_id"

# Calculate the character length of each unique ASV sequence
repseqs_12SV5_length<-repseqs_12SV5
repseqs_12SV5_length$seq_length <- nchar(repseqs_12SV5$sequence)
hist(repseqs_12SV5_length$seq_length)

# Remove ASVs >20% longer than the expected Riaz 12SV5 amplicon size (110 bp → upper limit 132 bp)
repseqs_12SV5_length<-repseqs_12SV5_length[!(repseqs_12SV5_length$seq_length>132),]

# Remove ASVs >20% shorter than the expected Riaz 12SV5 amplicon size (110 bp → lower limit 88 bp)
repseqs_12SV5_length<-repseqs_12SV5_length[!(repseqs_12SV5_length$seq_length<88),]
hist(repseqs_12SV5_length$seq_length)

# Save lookup table: links SEQ IDs back to their original sequence lengths (needed for joining later)
write.table(repseqs_12SV5_length,"R_outputs/full_12SV5_length_filtered_SEQ_to_sequence_lookup_table.txt",quote = FALSE)

# Format each entry as a two-line FASTA record (header\nsequence) by replacing the space separator with a newline
repseqs_12SV5_length$fasta_entry <- paste(repseqs_12SV5_length$fasta_id, 
                                          repseqs_12SV5_length$sequence, 
                                          sep = "\n")

# Retain a standalone lookup dataframe (FASTA ID ↔ sequence length) for use after BLAST
length_adjust_fasta_id<-repseqs_12SV5_length$fasta_id
length_adjust_fasta_id<-data.frame(length_adjust_fasta_id)
length_adjust_fasta_id$seq_length<-repseqs_12SV5_length$seq_length

repseqs_12SV5_length_min <- repseqs_12SV5_length
repseqs_12SV5_length_min$fasta_id<-NULL
repseqs_12SV5_length_min$seq_length<-NULL

# Export the length-filtered sequences as a FASTA-formatted CSV: this is the BLAST input file
write.csv(repseqs_12SV5_length_min,"R_outputs/full_12SV5_length_filtered_repseqs_12SV5_BLAST_INPUT.fasta", row.names=FALSE,quote = FALSE)

# --> Now BLAST this on the cluster using the "12SV5 gene" database:

# INSERT HERE Sophie's Classifier. May need to install on cluster in new conda
#________________________________________________________
# SECTION 3: SAMPLE ABUNDANCE FILTERING FROM table.qza
# Reads QIIME2 artifact, extracts BIOM abundance table,
# drops low-depth samples (<1,000 reads).
# NOTE: This stream is independent of Stream 1 above.
# The two streams converge downstream when BLAST output is
# merged back with the abundance table. Extract and read table.qza
# is an alternative way to get your FASTA data if need be, 
# without having to rely on QIIME2 View----HOWEVER, here I am using a 
# tsv file that was converted in qiime itself and with biom to a tsv. 
# I'm much comfier with this approach.

#qza_path <- "qza_conversion_s25_test/s25_table.tsv"
#tmp_dir <- tempdir()
#unzip(qza_path, exdir = tmp_dir)
# Find the BIOM file inside
#biom_path <- list.files(tmp_dir, pattern = "\\.biom$", 
                     #   recursive = TRUE, full.names = TRUE)[1]
# Read and convert to matrix
#biom_data <- read_biom(biom_path)
#seqtab_12SV5 <- as.data.frame(as.matrix(biom_data(biom_data)))
# BIOM format is ASVs as rows, samples as columns — transpose to match your original format
#seqtab_12SV5 <- t(seqtab_12SV5)
# Make it a proper data frame with sample names retained
#seqtab_12SV5 <- as.data.frame(seqtab_12SV5)
# Load the DADA2 ASV abundance table (samples as rows, sequences as columns)

seqtab_12SV5 <- table_full

# We want our sequences to become rows and samples to become columns, then tidy rownames. 
# My samples are already in this format, and so I don't need to run things this way:
# seqtab_12SV5_transpose<-t(seqtab_12SV5)
# seqtab_12SV5_transpose<-as.data.frame(seqtab_12SV5_transpose) 
# seqtab_12SV5_transpose<-setDT(seqtab_12SV5_transpose, keep.rownames = TRUE) 
# names(seqtab_12SV5_transpose)[names(seqtab_12SV5_transpose) == "rn"] <- "fasta_entry" 

# If you need to strip all double quote characters from every cell, run the following.
# seqtab_12SV5_transpose[] <- lapply(seqtab_12SV5_transpose, gsub, pattern='"', replacement='')

seqtab_12SV5<-as.data.frame(seqtab_12SV5)
seqtab_12SV5<-seqtab_12SV5 %>% rename(fasta_id = OTU.ID)

# Convert read counts to numeric
seqtab_12SV5<-data.frame(seqtab_12SV5)
i <- c(2:208) 
seqtab_12SV5[ , i] <- apply(seqtab_12SV5[ , i], 2, function(x) as.numeric(as.character(x)))

# Identify and drop samples with fewer than 1,000 total reads (too shallow for 
# reliable community estimation)
# I'm not gonna do this for now. 
# rarefy_col_sum<-colSums(seqtab_12SV5_transpose[c(2:10)])
# rarefy_col_sum<-data.frame(rarefy_col_sum)
# rarefy_col_sum<-setDT(rarefy_col_sum, keep.rownames = TRUE)[]
# rarefy_col_sum<-data.frame(rarefy_col_sum)
# rarefy_col_sum<-rarefy_col_sum[!(rarefy_col_sum$rarefy_col_sum>=1000),]
# drop<-rarefy_col_sum$rn
# seqtab_12SV5_transpose<-seqtab_12SV5_transpose[,!(names(seqtab_12SV5_transpose) %in% drop)]

#-------------------------------------------------------------------------
# Quick initial peek into our BLAST output
#-------------------------------------------------------------------------

blast_12SV5 <- `12SV5_blast_full_allASVs`

colnames(blast_12SV5) <- c(
  "fasta_id",    # 1  - ASV hash
  "sseqid",    # 2  - Subject accession
  "pident",    # 3  - % identity
  "length",    # 4  - Alignment length
  "mismatch",  # 5  - Mismatches
  "gapopen",   # 6  - Gap openings
  "qstart",    # 7  - Query alignment start
  "qend",      # 8  - Query alignment end
  "sstart",    # 9  - Subject alignment start
  "send",      # 10 - Subject alignment end
  "evalue",    # 11 - E-value
  "bitscore",  # 12 - Bit score
  "qlen",      # 13 - Query length
  "slen",      # 14 - Subject length
  "staxids",   # 15 - Taxonomy ID
  "sscinames", # 16 - Scientific name (genus + species)
  "scomnames", # 17 - Common name / subspecies ("aimara")
  "sskingdoms" # 18 - Kingdom (Eukaryota)
)

blast_12SV5 <- blast_12SV5 %>%
  mutate(genus = word(sscinames, 1))

# Run taxonkit locally in computer to get all the taxonomic details for your species ids
# Extract unique taxids from column 15
# cut -d',' -f15 /Volumes/ROSALIND/s25_blast.txt | sort -u > /Volumes/ROSALIND/taxonkit/taxids.txt

# Get full lineage
# taxonkit lineage /Volumes/ROSALIND/taxonkit/taxids.txt \
# --data-dir /Volumes/ROSALIND/taxonkit/ \
# --o /Volumes/ROSALIND/taxonkit/taxids_lineage.txt

# Reformat into ranked columns
# taxonkit reformat /Volumes/ROSALIND/taxonkit/taxids_lineage.txt \
# --data-dir /Volumes/ROSALIND/taxonkit/ \
# --format "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" \
# --fill-miss-rank \
# --o /Volumes/ROSALIND/taxonkit/taxids_reformatted.txt

# taxonkit lineage taxids.txt \
# --data-dir /Volumes/ROSALIND/taxonkit/ \
# --o taxids_lineage.txt

# taxonkit reformat taxids_lineage.txt \
# --data-dir /Volumes/ROSALIND/taxonkit/ \
# --format "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" \
# --fill-miss-rank \
# --o taxids_reformatted.txt

# Filter on thresholds, and then remove human contamination.
# Should the BLAST Thresholding be done first or last?

#--------------START JACK'S EARLY PIPELINE----------------------------------

filtered_12SV5 <- blast_12SV5 %>%
  filter(blast_12SV5$evalue <= 1e-3) %>%      # Keep only e-values <= 1e-3
  filter(blast_12SV5$pident >= 98) %>% # Keep only matches with >= 95% identity
  filter(sscinames != "Homo sapiens") # Remove human contam

sum(blast_12SV5$sscinames == "Homo sapiens")
sum(filtered_12SV5$sscinames == "Homo sapiens")
species_list <- unique(filtered_12SV5$sscinames)

species_counts_12SV5 <- filtered_12SV5 %>%
  count(sscinames, sort = TRUE) %>%
  rename(Species = sscinames, Count = n)

best_hits_12SV5 <- filtered_12SV5 %>%
  group_by(fasta_id) %>%
  slice_min(evalue, n = 5, with_ties = TRUE)

#-----------END JACK'S EARLY PIPELINE----------------------------------------
#---------------------------------------------------------------------------
# SECTION 4: SEQTAB IMPORT AND CLEANING
# Join BLAST taxonomy onto the abundance table, using the raw sequence string as the key;
# this replaces literal sequences with species assignments
# Row names (ASV sequences) are moved into a column called 'Seq' so the
# data frame can be manipulated without relying on implicit row names.
# Stray quotation marks from file export are stripped via gsub().
# A helper function (header.true) promotes the first row of a data frame
# to column names and drops it; this pattern recurs throughout the script
# wherever a transpose leaves the header in row 1.
#-------------------------------------------------------------------------------
      
seqtab_12SV5_tax <- best_hits_12SV5 %>% full_join(seqtab_12SV5, by = "fasta_id")
seqtab_12SV5_tax$rn <- NULL
seqtab_12SV5_tax$sseqid <- NULL
seqtab_12SV5_tax$qstart <- NULL
seqtab_12SV5_tax$qend <- NULL
seqtab_12SV5_tax$sstart <- NULL
seqtab_12SV5_tax$send <- NULL
seqtab_12SV5_tax$qlen <- NULL
seqtab_12SV5_tax$slen <- NULL
seqtab_12SV5_tax<-na.omit(seqtab_12SV5_tax)

i <- c(2:8, 13:219) 
seqtab_12SV5_tax[ , i] <- apply(seqtab_12SV5_tax[ , i], 2, function(x) as.numeric(as.character(x)))

# Calculate a per-sample 0.05% read threshold: ASVs below this proportion will be treated as noise
seqtab_12SV5_tax_col_sum<-colSums(seqtab_12SV5_tax[c(13:219)])
seqtab_12SV5_tax_col_sum<-data.frame(seqtab_12SV5_tax_col_sum)
seqtab_12SV5_tax_col_sum$read_filter<-seqtab_12SV5_tax_col_sum$seqtab_12SV5_tax_col_sum *0.0005
seqtab_12SV5_tax_col_sum<- seqtab_12SV5_tax_col_sum %>% rownames_to_column("ID")

# Remove ASVs with no BLAST species assignment, but we don't want to do this
# seqtab_12SV5_tax<-seqtab_12SV5_tax %>% drop_na(sscinames)

# =============================================================================
# SECTION 5: LONG FORMAT, METADATA JOIN & READ-DEPTH TRACKING PLOT
# The wide ASV table (rows=ASVs, cols=samples) is pivoted to long format
# (one row per ASV x sample observation) using tidyr::gather(). This format
# is required for ggplot2 stacked bar charts and for applying the per-sample
# read threshold row-by-row.
# Sample metadata (site codes, dates, seasons, env. covariates) is read in
# and joined. A corrmorant pairwise correlation matrix is generated for the
# WP2A environmental covariates to identify collinear predictors before
# building PERMANOVA and GAM models.
# The 0.05% read threshold is joined and applied: any ASV-sample
# row where reads <= threshold is removed.
# The filtered long table is written to disk as a checkpoint, then converted
# back to wide format for normalisation.
# A read-depth tracking plot is generated showing reads retained at each
# of the filter stages,
# faceted by sample with bars coloured by filter step.
# =============================================================================
#CONVERTING ASV TABLE INTO LONG FORMAT

# Reshape to long format (one row per ASV × sample combination) for filtering and metadata merging
seqtab_12SV5_tax_long <- seqtab_12SV5_tax %>%
  pivot_longer(
    -c(fasta_id, pident, length, mismatch, gapopen, evalue, 
       bitscore, staxids, sscinames, scomnames, sskingdoms, genus),
    names_to = "Sample",
    values_to = "reads"
  ) %>%
  mutate(reads = as.numeric(reads))

# Load project metadata and merge onto the long-format table

metadata_full <- `metadata_full`

metadata_summary <- metadata_full %>%
  select(sample_id, Month.Collected, Village, SubLocation, 
         GPS.Location, Liters.Pumped, pH..1..3.and.5.min., 
         Temp..1..3.and.5.min., EC..1..3.and.5.min., Region)

full_metadata_summary<-data.frame(metadata_summary)
write.csv(full_metadata_summary,"metadata_summary.csv")

metadata_full<-metadata_full %>% rename(sample_id = sample.id)
seqtab_12SV5_tax_long<-seqtab_12SV5_tax_long %>% rename(sample_id = Sample)
seqtab_12SV5_tax_col_sum<-seqtab_12SV5_tax_col_sum %>% rename(sample_id = ID)
seqtab_12SV5_tax_long$sample_id <- gsub("\\.", "-", seqtab_12SV5_tax_long$sample_id)
seqtab_12SV5_tax_col_sum$sample_id <- gsub("\\.", "-", seqtab_12SV5_tax_col_sum$sample_id)

seqtab_12SV5_tax_long_meta <- merge(metadata_full[, c("sample_id", "Month.Collected", 
"Village", "SubLocation", "GPS.Location", "Liters.Pumped", "pH..1..3.and.5.min.",
"Temp..1..3.and.5.min.", "EC..1..3.and.5.min.", "Control.")],seqtab_12SV5_tax_long, by="sample_id")
sum_reads<- sum(seqtab_12SV5_tax_long_meta$reads)

# Attach per-sample read thresholds
seqtab_12SV5_tax_long_meta_clean_filtered <- merge(seqtab_12SV5_tax_col_sum[, c("sample_id", "read_filter")],seqtab_12SV5_tax_long_meta, by="sample_id")

# Apply read depth filters: remove rows below the 0.05% per-sample threshold and below the 20-read hard floor
seqtab_12SV5_tax_long_meta_clean_filtered<-seqtab_12SV5_tax_long_meta_clean_filtered[!(seqtab_12SV5_tax_long_meta_clean_filtered$reads<=seqtab_12SV5_tax_long_meta_clean_filtered$read_filter),]
seqtab_12SV5_tax_long_meta_clean_filtered<-seqtab_12SV5_tax_long_meta_clean_filtered[!(seqtab_12SV5_tax_long_meta_clean_filtered$reads<20),]

write.csv(seqtab_12SV5_tax_long_meta_clean_filtered,"full_ASV_table_long_filtered.csv")
seqtab_12SV5_tax_long_meta_clean_filtered <- read_csv("full_ASV_table_long_filtered.csv")

# Pivot to wide format (species × sample) for normalisation
seqtab_12SV5_tax_wide_filtered<-seqtab_12SV5_tax_long_meta_clean_filtered
seqtab_12SV5_tax_wide_filtered$Order <- NULL 
seqtab_12SV5_tax_wide_filtered$read_filter <- NULL 
seqtab_12SV5_tax_wide_filtered$Seq_length <- NULL 

seqtab_12SV5_tax_wide_filtered<-seqtab_12SV5_tax_wide_filtered %>% 
  group_by(sample_id,sscinames) %>% 
  summarize(reads=sum(reads)) %>%
  rename(sample_id=sample_id,sscinames=sscinames, reads=reads,)

seqtab_12SV5_tax_wide_filtered <- spread(seqtab_12SV5_tax_wide_filtered,sample_id,reads)
seqtab_12SV5_tax_wide_filtered[is.na(seqtab_12SV5_tax_wide_filtered)] <- 0
seqtab_12SV5_tax_wide_filtered <- data.frame(t(seqtab_12SV5_tax_wide_filtered))

names(seqtab_12SV5_tax_wide_filtered) <- as.matrix(seqtab_12SV5_tax_wide_filtered[1, ])
seqtab_12SV5_tax_wide_filtered <- seqtab_12SV5_tax_wide_filtered[-1, ]
seqtab_12SV5_tax_wide_filtered[] <- lapply(seqtab_12SV5_tax_wide_filtered, function(x) type.convert(as.character(x)))

# Normalise: divide each ASV's read count by the sample total to get relative proportions.
# NOTE: The column index range (1:68) may need to be updated if the number of species changes.
seqtab_12SV5_tax_col_sum2<-rowSums(seqtab_12SV5_tax_wide_filtered[c(1:500)])
seqtab_12SV5_tax_col_sum2<-data.frame(seqtab_12SV5_tax_col_sum2)
seqtab_12SV5_tax_col_sum2$sample_id <- rownames(seqtab_12SV5_tax_col_sum2)

seqtab_12SV5_tax_normalised<- merge(seqtab_12SV5_tax_col_sum2[, c("sample_id", "seqtab_12SV5_tax_col_sum2")],seqtab_12SV5_tax_long_meta_clean_filtered, by="sample_id")
     
seqtab_12SV5_tax_normalised <- transform(seqtab_12SV5_tax_normalised, normalised_reads = reads / seqtab_12SV5_tax_col_sum2)

ggplot(seqtab_12SV5_tax_normalised , aes(x = sequence_id, fill = species, y = normalised_reads)) + 
  geom_bar(stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))+ theme(legend.position = "none")
                                      
#################################################################################
####### COMMUNITY COMP CHECK ###################################################
###############################################################################

# Quick visualisation of normalised community composition across all samples
top_n <- 25
top_species <- seqtab_12SV5_tax_normalised %>%
  group_by(sscinames) %>%
  summarise(total = sum(normalised_reads, na.rm = TRUE)) %>%
  slice_max(total, n = top_n) %>%
  pull(sscinames)

seqtab_12SV5_tax_normalised <- seqtab_12SV5_tax_normalised %>%
  mutate(species_plot = ifelse(sscinames %in% top_species, sscinames, "Other"))


ggplot(seqtab_12SV5_tax_normalised , aes(x = sample_id, fill = species_plot, y = normalised_reads)) + 
  geom_bar(stat = "identity", colour = "white")+
  theme(axis.text.x = element_text(angle = 90),
        legend.text = element_text(face = "italic", size = 9))

####################################################################################################

#COVARIATES IN METADATA PLOT FOR WP2A
library(corrmorant)
lofresh_metadata_WP2A<-lofresh_metadata_WP2_3 %>% filter(WP == "WP2A")
lofresh_metadata_WP2A<-lofresh_metadata_WP2A %>%
  select("season","Days","Dist_from_lake","pH",
"conductivity_uS_cm","gran_alk_ueq_L","daily_flow_by_catchment","monthly_flow_by_catchment","seasonal_flow_by_catchment",
"week_before_flow_by_catchment","river_width_m","river_gradient_percent","rainfall_monthly_average",
"temp_monthly_average")
corrmorant(lofresh_metadata_WP2A, style = "blue_red")
                                           
# Pivot normalised data to wide format and save as the master normalised community matrix
seqtab_12SV5_tax_normalised_wide<-seqtab_12SV5_tax_normalised
seqtab_12SV5_tax_normalised_wide$read_filter <- NULL 
seqtab_12SV5_tax_normalised_wide$seqtab_12SV5_tax_col_sum2 <- NULL 
seqtab_12SV5_tax_normalised_wide$seqtab_12SV5_tax_col_sum1 <- NULL 
seqtab_12SV5_tax_normalised_wide$Month.Collected <- NULL 
seqtab_12SV5_tax_normalised_wide$seqtab_12SV5_tax <- NULL 
seqtab_12SV5_tax_normalised_wide$GPS.Location <- NULL 
seqtab_12SV5_tax_normalised_wide$Liters.Pumped <- NULL 
seqtab_12SV5_tax_normalised_wide$Control. <- NULL 
seqtab_12SV5_tax_normalised_wide$pindent <- NULL 
seqtab_12SV5_tax_normalised_wide$length <- NULL 
seqtab_12SV5_tax_normalised_wide$mismatch <- NULL 
seqtab_12SV5_tax_normalised_wide$gapopen <- NULL 
seqtab_12SV5_tax_normalised_wide$evalue <- NULL 
seqtab_12SV5_tax_normalised_wide$bitscore <- NULL
seqtab_12SV5_tax_normalised_wide$species_plot <- NULL 

seqtab_12SV5_tax_normalised_wide<-seqtab_12SV5_tax_normalised_wide %>% 
  group_by(sample_id,sscinames) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(sample_id=sample_id,sscinames=sscinames, normalised_reads= normalised_reads)

seqtab_12SV5_tax_normalised_wide <- spread(seqtab_12SV5_tax_normalised_wide,sample_id,normalised_reads)
seqtab_12SV5_tax_normalised_wide[is.na(seqtab_12SV5_tax_normalised_wide)] <- 0
seqtab_12SV5_tax_normalised_wide <- data.frame(t(seqtab_12SV5_tax_normalised_wide))

names(seqtab_12SV5_tax_normalised_wide) <- as.matrix(seqtab_12SV5_tax_normalised_wide[1, ])
seqtab_12SV5_tax_normalised_wide <- seqtab_12SV5_tax_normalised_wide[-1, ]
seqtab_12SV5_tax_normalised_wide[] <- lapply(seqtab_12SV5_tax_normalised_wide, function(x) type.convert(as.character(x)))

write.csv(seqtab_12SV5_tax_normalised_wide,"FULL_ASV_table_wide_filtered.csv")

# Split by work package and add further metadata fields; save each sub-dataset as its own CSV

# Extract Set designation from sample_id
guyana_long <- seqtab_12SV5_tax_long_meta_clean_filtered %>%
  left_join(
    metadata_full %>% select(sample_id, Month.Collected, Village, SubLocation,
                             GPS.Location, Liters.Pumped, pH..1..3.and.5.min.,
                             Temp..1..3.and.5.min., EC..1..3.and.5.min., Region),
    by = "sample_id"
  ) %>%
  mutate(Set = sub(".*-(Set\\d+)$", "\\1", sample_id))  # extracts "Set1", "Set2", "Set3"

# Split by Village > SubLocation > Set and save each as its own CSV
guyana_long <- guyana_long %>%
  select(-ends_with(".y"), -`...1`) %>%
  rename_with(~ sub("\\.x$", "", .), ends_with(".x"))

# =============================================================================
# WORK PACKAGE SPLITTING
# The full filtered long dataset is divided into sub-datasets by work package
# and for Mawrow Creek, further by geographic region / river
# region. Each sub-dataset has the same set of metadata columns  attached via
# sequential merge calls, then saved as a standalone CSV file.
# These per-village CSVs are the starting point for all downstream analyses
# onwards.
# =============================================================================
                                             
guyana_edna_long_meta_clean <- read_csv("FULL_ASV_table_wide_filtered.csv")

# Surama
guyana_Surama <- seqtab_12SV5_tax_long_meta_clean_filtered %>%
  filter(Village == "Surama", Month.Collected == "April") %>%
  merge(metadata_full[, c("sample_id", "Control")], by = "sample_id") %>%
  merge(metadata_full[, c("sample_id", "SubLocation")], by = "sample_id") %>%
  merge(metadata_full[, c("sample_id", "pH")], by = "sample_id") %>%                                             
  merge(metadata_full[, c("sample_id", "Temp")], by = "sample_id") %>%
  merge(metadata_full[, c("sample_id", "EC")], by = "sample_id") %>%  
  merge(metadata_full[, c("sample_id", "Liters_Pumped")], by = "sample_id") %>%                                           
  merge(metadata_full[, c("sample_id", "GPS.Location")], by = "sample_id")
write.csv(guyana_Surama_April, "Surama_ASV_table_long_filtered.csv", row.names = FALSE)
guyana_Surama_April <- read_csv("Surama_ASV_table_long_filtered.csv")

# Kamwatta
guyana_Kamwatta_April <- seqtab_12SV5_tax_long_meta_clean_filtered %>%
  filter(Village == "Kamwatta", Month.Collected == "April") %>%
  merge(metadata_full[, c("sample_id", "Control")], by = "sample_id") %>%
  merge(metadata_full[, c("sample_id", "SubLocation")], by = "sample_id") %>%
  merge(metadata_full[, c("sample_id", "pH")], by = "sample_id") %>%                                             
  merge(metadata_full[, c("sample_id", "Temp")], by = "sample_id") %>%
  merge(metadata_full[, c("sample_id", "EC")], by = "sample_id") %>%  
  merge(metadata_full[, c("sample_id", "Liters_Pumped")], by = "sample_id") %>%                                           
  merge(metadata_full[, c("sample_id", "GPS.Location")], by = "sample_id")
write.csv(guyana_Kamwatta_April, "Kamwatta_April_ASV_table_long_filtered.csv", row.names = FALSE)
guyana_Kamwatta_April <- read_csv("Kamwatta_April_ASV_table_long_filtered.csv")

# Mainstay
guyana_Mainstay <- seqtab_12SV5_tax_long_meta_clean_filtered %>%
  filter(Village == "Mainstay", Month.Collected == "April") %>%
  merge(metadata_full[, c("sample_id", "Control")], by = "sample_id") %>%
  merge(metadata_full[, c("sample_id", "SubLocation")], by = "sample_id") %>%
  merge(metadata_full[, c("sample_id", "pH")], by = "sample_id") %>%                                             
  merge(metadata_full[, c("sample_id", "Temp")], by = "sample_id") %>%
  merge(metadata_full[, c("sample_id", "EC")], by = "sample_id") %>%  
  merge(metadata_full[, c("sample_id", "Liters_Pumped")], by = "sample_id") %>%                                           
  merge(metadata_full[, c("sample_id", "GPS.Location")], by = "sample_id")
write.csv(guyana_Mainstay, "Mainstay_April_ASV_table_long_filtered.csv", row.names = FALSE)
guyana_Mainstay <- read_csv("Mainstay_April_ASV_table_long_filtered.csv")                                           

# Wowetta Rainy & Dry
guyana_Wowetta <- seqtab_12SV5_tax_long_meta_clean_filtered %>%
  filter(Village == "Wowetta") %>%
  merge(metadata_full[, c("sample_id", "Control")], by = "sample_id") %>%
  merge(metadata_full[, c("sample_id", "SubLocation")], by = "sample_id") %>%
  merge(metadata_full[, c("sample_id", "pH")], by = "sample_id") %>%                                             
  merge(metadata_full[, c("sample_id", "Temp")], by = "sample_id") %>%
  merge(metadata_full[, c("sample_id", "EC")], by = "sample_id") %>%  
  merge(metadata_full[, c("sample_id", "Liters_Pumped")], by = "sample_id") %>%                                           
  merge(metadata_full[, c("sample_id", "GPS.Location")], by = "sample_id") %>%
  merge(metadata_full[, c("sample_id", "Month.Collected")], by = "sample_id")                                           
write.csv(guyana_Wowetta, "Wowetta_RainPlusDry_ASV_table_long_filtered.csv", row.names = FALSE)
guyana_Wowetta <- read_csv("Wowetta_RainPlusDry_ASV_table_long_filtered.csv")

# MowrowCreek
guyana_MowrowCreek <- seqtab_12SV5_tax_long_meta_clean_filtered %>%
  filter(Village == "MowrowCreek", Month.Collected == "July") %>%
  merge(metadata_full[, c("sample_id", "Control")], by = "sample_id") %>%
  merge(metadata_full[, c("sample_id", "SubLocation")], by = "sample_id") %>%
  merge(metadata_full[, c("sample_id", "pH")], by = "sample_id") %>%                                             
  merge(metadata_full[, c("sample_id", "Temp")], by = "sample_id") %>%
  merge(metadata_full[, c("sample_id", "EC")], by = "sample_id") %>%  
  merge(metadata_full[, c("sample_id", "Liters_Pumped")], by = "sample_id") %>%                                           
  merge(metadata_full[, c("sample_id", "GPS.Location")], by = "sample_id")
write.csv(guyana_MowrowCreek, "MowrowCreek_July_ASV_table_long_filtered.csv", row.names = FALSE)
guyana_MowrowCreek <- read_csv("MowrowCreek_July_ASV_table_long_filtered.csv")


# Whitewater
guyana_Whitewater <- seqtab_12SV5_tax_long_meta_clean_filtered %>%
  filter(Village == "Whitewater", Month.Collected == "April") %>%
  merge(metadata_full[, c("sample_id", "Control")], by = "sample_id") %>%
  merge(metadata_full[, c("sample_id", "SubLocation")], by = "sample_id") %>%
  merge(metadata_full[, c("sample_id", "pH")], by = "sample_id") %>%                                             
  merge(metadata_full[, c("sample_id", "Temp")], by = "sample_id") %>%
  merge(metadata_full[, c("sample_id", "EC")], by = "sample_id") %>%  
  merge(metadata_full[, c("sample_id", "Liters_Pumped")], by = "sample_id") %>%                                           
  merge(metadata_full[, c("sample_id", "GPS.Location")], by = "sample_id")
write.csv(guyana_Whitewater, "Whitewater_April_ASV_table_long_filtered.csv", row.names = FALSE)
guyana_Whitewater <- read_csv("Whitewater_April_ASV_table_long_filtered.csv")

# Mawrow Creek Upper, Middle, and Lower regions can be accounted for in with the SubRegions variable in the Metadata

# Village Rainy vs Dry can be compared in the Wowetta dataset, which holds both seasons. The Kamwatta and Whitewater datasets 
# could potentially be used as the dry season analogs for Mowrow Creek, which was a Rainy season dataset. However, I also like
# the idea of keeping the Mowrow Creek dataset separate from rainy vs dry analyses.





############################################################
# Extras ---------------------
############################################################                                             
                                             
                        
# Convert QIIME FASTA to txt
header_idx <- grep("^>", lines)
headers <- sub("^>", "", lines[header_idx])  # remove the ">" prefix

sequences <- sapply(seq_along(header_idx), function(i) {
  start <- header_idx[i] + 1
  end <- ifelse(i < length(header_idx), header_idx[i + 1] - 1, length(lines))
  paste(lines[start:end], collapse = "")
})

# ---- Build data frame and write to .tsv ----
df <- data.frame(header = headers, sequence = sequences, stringsAsFactors = FALSE)
write.table(df, "repseqs_12SV5.txt", sep = "\t", row.names = FALSE, quote = FALSE)

repseqs_12SV5 <- read_table("repseqs_12SV5.txt")
colnames(repseqs_12SV5)[colnames(repseqs_12SV5) == "header"] <- "fasta_id"
