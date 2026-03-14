library(tidyverse)
library(RColorBrewer)
library(biomformat)
library(readr)
library(data.table)
library(vegan)
library(plotly)
#library(microDecon)
library(ggpubr)
library(effectsize)
library(see)
library(lme4)
library(lmerTest)
library(emmeans)
library(lsmeans)
library(readr)
library(gratia)
setwd("/Users/jacksadler/Desktop/project_thesis/taxonomy/12SV5_BLAST")

#________________________________________________________
# Removing 12S sequences in the taxtable that are the wrong size, then BLAST'ing
# with a broad 12S database.
#________________________________________________________

# -----------Read and parse FASTA file
fasta_path <- "repseqs_12SV5.fasta"
lines <- readLines(fasta_path)

#------------Extract headers and sequences

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

# Strip pre-existing taxonomy columns — only the raw sequences are needed for length filtering
repseqs_12SV5$Kingdom <-NULL
repseqs_12SV5$Phylum <-NULL
repseqs_12SV5$Class <-NULL
repseqs_12SV5$Order <-NULL
repseqs_12SV5$Family <-NULL
repseqs_12SV5$Genus <-NULL
repseqs_12SV5$Species <-NULL

# Calculate the character length of each unique ASV sequence
repseqs_12SV5_length<-repseqs_12SV5
repseqs_12SV5_length$seq_length <- nchar(repseqs_12SV5$sequence)
hist(repseqs_12SV5_length$seq_length)

# Remove ASVs >20% longer than the expected Riaz 12SV5 amplicon size (110 bp → upper limit 132 bp)
repseqs_12SV5_length<-repseqs_12SV5_length[!(repseqs_12SV5_length$seq_length>132),]

# Remove ASVs >20% shorter than the expected Riaz 12SV5 amplicon size (110 bp → lower limit 88 bp)
repseqs_12SV5_length<-repseqs_12SV5_length[!(repseqs_12SV5_length$seq_length<88),]
hist(repseqs_12SV5_length$seq_length)

#We should be reading directly from a .fasta file, so FASTA id's don't need to be created.
# Assign sequential FASTA-style IDs (>SEQ1, >SEQ2, ...) to the length-filtered sequences
# taxTab_WP2_12S_length$Seq_number <- 1:nrow(taxTab_WP2_12S_length) 
# taxTab_WP2_12S_length$Seq<-"SEQ"
# taxTab_WP2_12S_length$div<-">"

# taxTab_WP2_12S_length$Seq_number <-paste(taxTab_WP2_12S_length$div,taxTab_WP2_12S_length$Seq,taxTab_WP2_12S_length$Seq_number)
# taxTab_WP2_12S_length$Seq_number = as.character(gsub(" ", "", taxTab_WP2_12S_length$Seq_number))

# taxTab_WP2_12S_length$Seq <-NULL
# taxTab_WP2_12S_length$div <-NULL

# NOTE: We may need this step

# Save lookup table: links SEQ IDs back to their original sequence lengths (needed for joining later)
write.table(repseqs_12SV5_length,"length_filtered_SEQ_to_sequence_lookup_table.txt",quote = FALSE)

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
write.csv(repseqs_12SV5_length_min,"length_filtered_repseqs_12SV5_BLAST_INPUT.fasta", row.names=FALSE,quote = FALSE)

# --> Now BLAST this on the cluster using the "12SV5 gene" database:
#     -perc_identity 70 -qcov_hsp_perc 80 -evalue 1

# INSERT HERE Sophie's Classifier. May need to install on cluster in new conda

# ---- Extract and read table.qza ----
qza_path <- "table.qza"
tmp_dir <- tempdir()

unzip(qza_path, exdir = tmp_dir)

# Find the BIOM file inside
biom_path <- list.files(tmp_dir, pattern = "\\.biom$", 
                        recursive = TRUE, full.names = TRUE)[1]

# Read and convert to matrix
biom_data <- read_biom(biom_path)
seqtab_12SV5 <- as.data.frame(as.matrix(biom_data(biom_data)))

# BIOM format is ASVs as rows, samples as columns — transpose to match your original format
seqtab_12SV5 <- t(seqtab_12SV5)

# Make it a proper data frame with sample names retained
seqtab_12SV5 <- as.data.frame(seqtab_12SV5)

# Load the DADA2 ASV abundance table if you have it as a .txt file (samples as rows, sequences as columns)
# seqtab_12SV5 <- read_table("seqtab_12SV5.txt")

# Transpose so that sequences become rows and samples become columns, then tidy rownames
seqtab_12SV5_transpose<-t(seqtab_12SV5)
seqtab_12SV5_transpose<-as.data.frame(seqtab_12SV5_transpose)
seqtab_12SV5_transpose<-setDT(seqtab_12SV5_transpose, keep.rownames = TRUE)
names(seqtab_12SV5_transpose)[names(seqtab_12SV5_transpose) == "rn"] <- "fasta_entry"
seqtab_12SV5_transpose[] <- lapply(seqtab_12SV5_transpose, gsub, pattern='"', replacement='')

# Convert read counts to numeric
seqtab_12SV5_transpose<-data.frame(seqtab_12SV5_transpose)
i <- c(2:10) 
seqtab_12SV5_transpose[ , i] <- apply(seqtab_12SV5_transpose[ , i], 2,  
                                           function(x) as.numeric(as.character(x)))

# Identify and drop samples with fewer than 1,000 total reads (too shallow for reliable community estimation)
rarefy_col_sum<-colSums(seqtab_12SV5_transpose[c(2:10)])
rarefy_col_sum<-data.frame(rarefy_col_sum)
rarefy_col_sum<-setDT(rarefy_col_sum, keep.rownames = TRUE)[]
rarefy_col_sum<-data.frame(rarefy_col_sum)

rarefy_col_sum<-rarefy_col_sum[!(rarefy_col_sum$rarefy_col_sum>=1000),]
drop<-rarefy_col_sum$rn

seqtab_12SV5_transpose<-seqtab_12SV5_transpose[,!(names(seqtab_12SV5_transpose) %in% drop)]

#-------------------------------------------------------------------------
# Quick initial peek into our BLAST output
#-------------------------------------------------------------------------

blast_12SV5 <- s25_blast

colnames(blast_12SV5) <- c(
  "qseqid",    # 1  - ASV hash
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

# Filter on thresholds, and then remove human contamination.

filtered_12SV5 <- blast_12SV5 %>%
  filter(blast_12SV5$evalue <= 1e-3) %>%      # Keep only e-values <= 1e-3
  filter(blast_12SV5$pident >= 98) %>% # Keep only matches with >= 95% identity
  filter(sscinames != "Homo sapiens") # Remove human contam

sum(filtered_12SV5$sscinames == "Homo sapiens")
unique(filtered_12SV5$sscinames)

species_counts_12SV5 <- filtered_12SV5 %>%
  count(sscinames, sort = TRUE) %>%
  rename(Species = sscinames, Count = n)

best_hits_12SV5 <- filtered_12SV5 %>%
  group_by(qseqid) %>%
  slice_min(evalue, n = 5, with_ties = TRUE)
ungroup()

# Join BLAST taxonomy onto the abundance table, using the raw sequence string as the key;
# this replaces literal sequences with species assignments
seqtab_12SV5_tax <- best_hits_12SV5 %>% full_join(seqtab_12SV5_transpose, by = c("qseqid" = "fasta_entry"))
seqtab_12SV5_tax$rn <- NULL
seqtab_12SV5_tax$sseqid <- NULL
seqtab_12SV5_tax$qstart <- NULL
seqtab_12SV5_tax$qend <- NULL
seqtab_12SV5_tax$sstart <- NULL
seqtab_12SV5_tax$send <- NULL
seqtab_12SV5_tax$qlen <- NULL
seqtab_12SV5_tax$slen <- NULL
seqtab_12SV5_tax<-na.omit(seqtab_12SV5_tax)

i <- c(2:8, 12:20) 
seqtab_12SV5_tax[ , i] <- apply(seqtab_12SV5_tax[ , i], 2,  
                                     function(x) as.numeric(as.character(x)))

# Calculate a per-sample 0.05% read threshold: ASVs below this proportion will be treated as noise
seqtab_12SV5_tax_col_sum<-colSums(seqtab_12SV5_tax[c(12:20)])
seqtab_12SV5_tax_col_sum<-data.frame(seqtab_12SV5_tax_col_sum)
seqtab_12SV5_tax_col_sum$read_filter<-seqtab_12SV5_tax_col_sum$seqtab_12SV5_tax_col_sum *0.0005
seqtab_12SV5_tax_col_sum<- seqtab_12SV5_tax_col_sum %>% rownames_to_column("ID")

# Remove ASVs with no BLAST species assignment
seqtab_12SV5_tax<-seqtab_12SV5_tax %>% drop_na(sscinames)

# Reshape to long format (one row per ASV × sample combination) for filtering and metadata merging
seqtab_12SV5_tax_long <- seqtab_12SV5_tax %>%
  pivot_longer(
    cols = starts_with("S_SC_"),
    names_to = "Sample",
    values_to = "reads"
  ) %>%
  mutate(reads = as.numeric(reads))

# Load project metadata and merge onto the long-format table
full_metadata_summary <- table(SC_Subset_metadata$Month.Collected,
                               SC_Subset_metadata$Village,
                               SC_Subset_metadata$SubLocation,
                               SC_Subset_metadata$GPS.Location,
                               SC_Subset_metadata$Region)
full_metadata_summary<-data.frame(full_metadata_summary)
write.csv(full_metadata_summary,"full_metadata_summary.csv")

SC_Subset_metadata<-SC_Subset_metadata %>% rename(sample_id = ID)
seqtab_12SV5_tax_long<-seqtab_12SV5_tax_long %>% rename(sample_id = ID)
seqtab_12SV5_tax_col_sum<-seqtab_12SV5_tax_col_sum %>% rename(sample_id = ID)
seqtab_12SV5_tax_long$sample_id <- gsub("\\.", "-", seqtab_12SV5_tax_long$sample_id)
seqtab_12SV5_tax_col_sum$sample_id <- gsub("\\.", "-", seqtab_12SV5_tax_col_sum$sample_id)

seqtab_12SV5_tax_long_meta <- merge(SC_Subset_metadata[, c("sample_id", "Month.Collected", "GPS.Location", "Liters.Pumped", "Control.")],seqtab_12SV5_tax_long, by="sample_id")
sum_reads<- sum(seqtab_12SV5_tax_long_meta$reads)

# Attach per-sample read thresholds
seqtab_12SV5_tax_long_meta_clean_filtered <- merge(seqtab_12SV5_tax_col_sum[, c("sample_id", "read_filter")],seqtab_12SV5_tax_long_meta, by="sample_id")

# Apply read depth filters: remove rows below the 0.05% per-sample threshold and below the 20-read hard floor
seqtab_12SV5_tax_long_meta_clean_filtered<-seqtab_12SV5_tax_long_meta_clean_filtered[!(seqtab_12SV5_tax_long_meta_clean_filtered$reads<=seqtab_12SV5_tax_long_meta_clean_filtered$read_filter),]
seqtab_12SV5_tax_long_meta_clean_filtered<-seqtab_12SV5_tax_long_meta_clean_filtered[!(seqtab_12SV5_tax_long_meta_clean_filtered$reads<20),]

write.csv(seqtab_12SV5_tax_long_meta_clean_filtered,"s25_ASV_table_long_filtered.csv")
seqtab_12SV5_tax_long_meta_clean_filtered <- read_csv("s25_ASV_table_long_filtered.csv")

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
seqtab_12SV5_tax_col_sum2<-rowSums(seqtab_12SV5_tax_wide_filtered[c(1:88)])
seqtab_12SV5_tax_col_sum2<-data.frame(seqtab_12SV5_tax_col_sum2)
seqtab_12SV5_tax_col_sum2$sample_id <- rownames(seqtab_12SV5_tax_col_sum2)

seqtab_12SV5_tax_normalised<- merge(seqtab_12SV5_tax_col_sum2[, c("sample_id", "seqtab_12SV5_tax_col_sum2")],seqtab_12SV5_tax_long_meta_clean_filtered, by="sample_id")

seqtab_12SV5_tax_normalised <- transform(seqtab_12SV5_tax_normalised, normalised_reads = reads / seqtab_12SV5_tax_col_sum2)

# Quick visualisation of normalised community composition across all samples
top_n <- 10
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

write.csv(seqtab_12SV5_tax_normalised_wide,"S25_ASV_table_wide_filtered.csv")

# Split by work package and add further metadata fields; save each sub-dataset as its own CSV
### WP2A  Splitting work packages and adding metadata ###
lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
seqtabNoC_WP2A_12S_tax_long_meta_clean<-seqtab_12SV5_tax_long_meta_clean_filtered %>% filter(WP == "WP2A")
seqtabNoC_WP2A_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "SampleSite_time")],seqtabNoC_WP2A_12S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2A_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "SampleSite_code")],seqtabNoC_WP2A_12S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2A_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "Date_sampled")],seqtabNoC_WP2A_12S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2A_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "Days")],seqtabNoC_WP2A_12S_tax_long_meta_clean, by="ID")

write.csv(seqtabNoC_WP2A_12S_tax_long_meta_clean,"WP2A_ASV_table_long_filtered_family.csv")
seqtabNoC_WP2A_12S_tax_long_meta_clean <- read_csv("WP2A_ASV_table_long_filtered_family.csv")

### WP2C Splitting work packages and adding metadata ###
seqtabNoC_WP2C_12S_tax_long_meta_clean<-seqtab_12SV5_tax_long_meta_clean_filtered %>% filter(WP == "WP2C")
seqtabNoC_WP2C_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "SampleSite_time")],seqtabNoC_WP2C_12S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2C_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "SampleSite_code")],seqtabNoC_WP2C_12S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2C_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "Date_sampled")],seqtabNoC_WP2C_12S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2C_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "Days")],seqtabNoC_WP2C_12S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2C_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "Region")],seqtabNoC_WP2C_12S_tax_long_meta_clean, by="ID")

# Further split WP2C by geographic region (river system)
seqtabNoC_WP2C_12S_tax_long_meta_clean_Gwash<-seqtabNoC_WP2C_12S_tax_long_meta_clean %>% filter(Region == "Leicestershire")
seqtabNoC_WP2C_12S_tax_long_meta_clean_Towy<-seqtabNoC_WP2C_12S_tax_long_meta_clean %>% filter(Region == "Carmarthenshire")
seqtabNoC_WP2C_12S_tax_long_meta_clean_Skan<-seqtabNoC_WP2C_12S_tax_long_meta_clean %>% filter(Region == "New_York")
seqtabNoC_WP2C_12S_tax_long_meta_clean_Glatt<-seqtabNoC_WP2C_12S_tax_long_meta_clean %>% filter(Region == "Switzerland")

write.csv(seqtabNoC_WP2C_12S_tax_long_meta_clean_Gwash,"WP2C_Gwash_ASV_table_long_filtered_family.csv")
write.csv(seqtabNoC_WP2C_12S_tax_long_meta_clean_Towy,"WP2C_Towy_ASV_table_long_filtered_family.csv")
write.csv(seqtabNoC_WP2C_12S_tax_long_meta_clean_Skan,"WP2C_Skan_ASV_table_long_filtered_family.csv")
write.csv(seqtabNoC_WP2C_12S_tax_long_meta_clean_Glatt,"WP2C_Glatt_ASV_table_long_filtered_family.csv")
