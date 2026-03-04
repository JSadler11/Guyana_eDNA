# Analysing BLAST Output — WP2 12S eDNA Pipeline

This document walks through the full analysis pipeline for the LOFRESH WP2 12S MiFish eDNA metabarcoding data. Code is presented in the order of the workflow, grouped by primary purpose. No code has been altered; headings and explanatory comments have been added throughout.

---

## Table of Contents

1. [Libraries & Setup](#1-libraries--setup)
2. [Stage 1 — ASV Length Filtering & BLAST Input Preparation](#2-stage-1--asv-length-filtering--blast-input-preparation)
3. [Stage 2 — Two-Pass BLAST: Non-fish Removal & Fish-Only Re-BLAST](#3-stage-2--two-pass-blast-non-fish-removal--fish-only-re-blast)
4. [Stage 3 — ASV Table Construction, Decontamination & Normalisation](#4-stage-3--asv-table-construction-decontamination--normalisation)
5. [Stage 4 — Per-River Biological Curation & Community Visualisation](#5-stage-4--per-river-biological-curation--community-visualisation)
   - [WP2A (Conwy)](#wp2a-conwy)
   - [WP2C Glatt (Switzerland)](#wp2c-glatt-switzerland)
   - [WP2C Towy (Carmarthenshire)](#wp2c-towy-carmarthenshire)
   - [WP2C Gwash (Leicestershire)](#wp2c-gwash-leicestershire)
   - [WP2C Skaneateles (New York)](#wp2c-skaneateles-new-york)
6. [Stage 5 — Alpha Diversity: Metrics, Modelling & Cross-River Comparison](#6-stage-5--alpha-diversity-metrics-modelling--cross-river-comparison)
7. [Stage 6 — Beta Diversity: PERMANOVA, NMDS & Multivariate Analysis](#7-stage-6--beta-diversity-permanova-nmds--multivariate-analysis)

---

## 1. Libraries & Setup

All required packages are loaded here. Key dependencies: `tidyverse` and `data.table` for data wrangling; `vegan` for ecological diversity metrics and ordination; `microDecon` for decontamination; `lme4`/`lmerTest`/`emmeans` for mixed-effects modelling; `mgcv`/`gratia` for GAMs; `plotly` for interactive plots.

```r
library(tidyverse)
library(readr)
library(data.table)
library(vegan)
library(plotly)
library(microDecon)
library(ggpubr)
library(effectsize)
library(see)
library(lme4)
library(lmerTest)
library(emmeans)
library(lsmeans)
library(gratia)
setwd("WP2_analysis/12S_BLAST")
```

---

## 2. Stage 1 — ASV Length Filtering & BLAST Input Preparation

### Purpose
The Riaz 12SV5 primer set targets a ~110 bp region of the 12S rRNA gene. ASVs that fall outside ±20% of this expected length (i.e. shorter than 88 bp or longer than 132 bp) are assumed to be artefacts and are removed before BLASTing. The surviving sequences are renamed with FASTA-style headers and exported as a `.fasta`-formatted file ready for submission to NCBI BLAST against the broad "12S rRNA gene" database.

> **Cluster BLAST parameters (Pass 1):** `-perc_identity 70 -qcov_hsp_perc 80 -evalue 1`
>
> These deliberately permissive settings are used to cast a wide net and identify any non-fish sequences before the stricter fish-only pass for the Perry paper. However, we are looking for ASVs that correspond to non-fish species and will only exclude non-native species, or human/human commensal ASVs.

```r
#________________________________________________________
# Removing 12S sequences in the taxtable that are the wrong size, then BLAST'ing
# with a broad 12S database.
#________________________________________________________

taxTab_WP2_12S <- read_table2("taxTab_WP2_12S.txt")

# Strip pre-existing taxonomy columns — only the raw sequences are needed for length filtering
repseqs_12SV5$Kingdom <-NULL
repseqs_12SV5$Phylum <-NULL
repseqs_12SV5$Class <-NULL
repseqs_12SV5$Order <-NULL
repseqs_12SV5$Family <-NULL
repseqs_12SV5$Genus <-NULL
repseqs_12SV5$Species <-NULL

# Calculate the character length of each unique ASV sequence
repseqs_12SV5_length<-aggregate(Seq~Seq_length, transform(repseqs_12SV5, Seq_length=Seq),
                                        FUN=function(x) nchar(unique(x)))
hist(repseqs_12SV5_length$Seq)

# Remove ASVs >20% longer than the expected MiFish amplicon size (110 bp → upper limit 132 bp)
repseqs_12SV5_length<-repseqs_12SV5_length[!(repseqs_12SV5_length$Seq>132),]

# Remove ASVs >20% shorter than the expected MiFish amplicon size (110 bp → lower limit 88 bp)
repseqs_12SV5_length<-repseqs_12SV5_length[!(repseqs_12SV5_length$Seq<88),]
hist(repseqs_12SV5_length$Seq)

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
repseqs_12SV5_length$Seq <-paste(repseqs_12SV5_length$Seq_number,repseqs_12SV5_length$Seq_length)
repseqs_12SV5_length$Seq <- as.character(gsub(" ", "\n", repseqs_12SV5_length$Seq))

# Retain a standalone lookup dataframe (SEQ ID ↔ sequence length) for use after BLAST
length_adjust_Seq_number<-repseqs_12SV5_length$Seq_number
length_adjust_Seq_number<-data.frame(length_adjust_Seq_number)
length_adjust_Seq_number$Seq_length<-repseqs_12SV5_length$Seq_length

repseqs_12SV5_length$Seq_number<-NULL
repseqs_12SV5_length$Seq_length<-NULL

# Export the length-filtered sequences as a FASTA-formatted CSV: this is the BLAST input file
write.csv(repseqs_12SV5_length,"length_filtered_taxTab_WP2_12S_BLAST_INPUT.fasta", row.names=FALSE,quote = FALSE)

# --> Now BLAST this on the cluster using the "12SV5 gene" database:
#     -perc_identity 70 -qcov_hsp_perc 80 -evalue 1
```

---

## 3. Stage 2 — Two-Pass BLAST: Non-fish Removal & Fish-Only Re-BLAST

### Purpose
The Perry paper handled BLAST results in two sequential passes. The first pass used a broad 12S database to identify and discard any ASVs assigned to non-fish taxa (e.g. mammals). Only sequences classified as ray-finned fish (`Actinopteri`) were retained. These fish-only sequences were then formatted for a second, stricter BLAST against a curated mitogenome database to achieve species-level taxonomic resolution.

For our purposes, we want to start with a specific BLAST analysis, and will then pass our fish ASVs through Sophie's classifier.

---

### Pass 1 — Read back broad BLAST results and filter to fish only

```r
#________________________________________________________
# Read in the length-adjusted BLAST output file from the first BLAST of NCBI.
# Requires objects created in Stage 1 to be present in the environment.
#________________________________________________________

length_filtered_BLAST_OUTPUT_12S <- read_table2("length_filtered_BLAST_OUTPUT_12S.txt", 
                                                col_names = FALSE)

# Bind the taxonomy (retrieved separately from NCBI) onto the BLAST hit table
TAXONOMY_filtered_BLAST_OUTPUT_12S <- read_csv("TAXONOMY_filtered_BLAST_OUTPUT_12S.csv")
length_filtered_BLAST_OUTPUT_12S<-cbind(length_filtered_BLAST_OUTPUT_12S,TAXONOMY_filtered_BLAST_OUTPUT_12S)

# Clean up the SEQ ID format in the lookup table to match the BLAST output column
length_adjust_Seq_number$length_adjust_Seq_number<-as.character(gsub(">", "", length_adjust_Seq_number$length_adjust_Seq_number))
length_adjust_Seq_number<-length_adjust_Seq_number %>% rename(Seq = "length_adjust_Seq_number")
length_filtered_BLAST_OUTPUT_12S<-length_filtered_BLAST_OUTPUT_12S %>% rename(Seq = "X1")

# Join BLAST results back onto the SEQ ID lookup table
length_adjust_Seq_number <- length_adjust_Seq_number %>% full_join(length_filtered_BLAST_OUTPUT_12S, by = c("Seq"))

# Remove ASVs that could not be assigned a taxon by the broad BLAST
length_adjust_Seq_number<-na.omit(length_adjust_Seq_number)

# Split into fish (Actinopteri) and mammal subsets; non-fish ASVs are discarded from the main pipeline
length_adjust_Seq_number_fish<-length_adjust_Seq_number %>% filter(class == "Actinopteri")
length_adjust_Seq_number_mammals<-length_adjust_Seq_number %>% filter(class == "Mammalia")
```

---

### Pass 1 → Pass 2 — Format fish-only sequences for Sophie's Classifier

```r
# Reformat the fish-only ASVs into FASTA structure for the second BLAST run
BLAST_2_fish<-length_adjust_Seq_number_fish$Seq
BLAST_2_fish<-data.frame(BLAST_2_fish)
BLAST_2_fish$Seq<-length_adjust_Seq_number_fish$Seq_length

BLAST_2_fish$Seq_BLAST <-paste(">",BLAST_2_fish$BLAST_2_fish,BLAST_2_fish$Seq)
BLAST_2_fish$Seq_BLAST<- gsub("> ", ">", BLAST_2_fish$Seq_BLAST)
BLAST_2_fish$Seq_BLAST<- as.character(gsub(" ", "\n", BLAST_2_fish$Seq_BLAST))

BLAST_2_fish$Seq<-NULL
BLAST_2_fish$BLAST_2_fish<-NULL

# Export fish-only sequences as input for Sophie's QIIME2 Classifier
write.csv(BLAST_2_fish,"length_filtered_repseqs_12SV5_BLAST_INPUT_FISH_ONLY.fasta", row.names=FALSE,quote = FALSE)
remove(BLAST_2_fish, length_adjust_Seq_number)

# --> Now run this on the cluster in QIIME2 2023:
#     -perc_identity 90 -qcov_hsp_perc 90 -evalue 0.001
#     (stricter thresholds appropriate for species-level assignment)
```

---

### Pass 2 — Read back Sophie's Classifier results and attach taxonomy

```r
#________________________________________________________
# Read in the fish-only Classifier output generated against the general BLAST database
#________________________________________________________

BLAST_OUTPUT_12S <- read_table2("BLAST_OUTPUT_FISH_ONLY_12S.txt", col_names = FALSE)
BLAST_OUTPUT_12S<-BLAST_OUTPUT_12S %>% rename(rn = X2 )

# I don't think we'll need this, but keep it here in case useful
# Strip formatting artefacts from NCBI reference accession numbers (e.g. "ref|NC_012345|" → "NC_012345")
# mitogenome_annotations <- read_csv("mitogenome_annotations.csv")
# BLAST_OUTPUT_12S$rn <- gsub("ref|", "", BLAST_OUTPUT_12S$rn)
# BLAST_OUTPUT_12S$rn <- gsub("\\|", "", BLAST_OUTPUT_12S$rn)

# Join species-level taxonomy from the mitogenome annotation table using accession number as the key
BLAST_OUTPUT_12S <- BLAST_OUTPUT_12S %>% full_join(mitogenome_annotations, by = c("rn"))
BLAST_OUTPUT_12S<-na.omit(BLAST_OUTPUT_12S)

# Rebuild a lookup linking each SEQ ID to its actual sequence (needed to join onto the abundance table)
SEQ_lookup<-length_adjust_Seq_number_fish$Seq
SEQ_lookup<-data.frame(SEQ_lookup)
SEQ_lookup$Seq<-length_adjust_Seq_number_fish$Seq_length
SEQ_lookup<-SEQ_lookup %>% rename(X1 = SEQ_lookup )

BLAST_OUTPUT_12S <- BLAST_OUTPUT_12S %>% full_join(SEQ_lookup, by = c("X1"))
BLAST_OUTPUT_12S<-na.omit(BLAST_OUTPUT_12S)

# Save the final annotated BLAST output: SEQ IDs with species taxonomy attached
write.table(BLAST_OUTPUT_12S,"BLAST_OUTPUT_12S_IDs_tax.txt",row.names = FALSE,col.names=FALSE,quote = FALSE)

# Clean up intermediate objects and free memory
remove(mitogenome_annotations,length_adjust_Seq_number_fish,length_filtered_BLAST_OUTPUT_12S,SEQ_lookup,
       TAXONOMY_filtered_BLAST_OUTPUT_12S,taxTab_WP2_12S,taxTab_WP2_12S_length)
gc()
```

---

## 4. Stage 3 — ASV Table Construction, Decontamination & Normalisation

### Purpose
The raw DADA2 sequence abundance table (`seqtabNoC`) is transposed, low-depth samples are dropped, and it is joined with the BLAST taxonomy to replace raw sequences with species names. Read counts are filtered using a per-sample 0.05% threshold and a hard floor of 20 reads to remove low-confidence detections. Reads are then normalised to relative proportions within each sample. The combined dataset is split by work package (WP2A and WP2C sub-rivers) and saved for all downstream analyses.

```r
# Load the DADA2 ASV abundance table (samples as rows, sequences as columns)
seqtabNoC_WP2_12S <- read_table2("seqtabNoC_WP2_12S.txt")

# Transpose so that sequences become rows and samples become columns, then tidy rownames
seqtabNoC_WP2_12S_transpose<-t(seqtabNoC_WP2_12S)
seqtabNoC_WP2_12S_transpose<-as.data.frame(seqtabNoC_WP2_12S_transpose)
seqtabNoC_WP2_12S_transpose<-setDT(seqtabNoC_WP2_12S_transpose, keep.rownames = TRUE)
names(seqtabNoC_WP2_12S_transpose)[names(seqtabNoC_WP2_12S_transpose) == "rn"] <- "Seq"
seqtabNoC_WP2_12S_transpose[] <- lapply(seqtabNoC_WP2_12S_transpose, gsub, pattern='"', replacement='')

# Promote the first row (sample names) to column headers
seqtabNoC_WP2_12S_transpose <- na.omit(transform(seqtabNoC_WP2_12S_transpose, Seq = c("Seq", Seq[-nrow(seqtabNoC_WP2_12S_transpose)])))
header.true <- function(seqtabNoC_WP2_12S_transpose) {
  names(seqtabNoC_WP2_12S_transpose) <- as.character(unlist(seqtabNoC_WP2_12S_transpose[1,]))
  seqtabNoC_WP2_12S_transpose[-1,]
}
seqtabNoC_WP2_12S_transpose<-header.true(seqtabNoC_WP2_12S_transpose)
remove(seqtabNoC_WP2_12S)

# Convert read counts to numeric
seqtabNoC_WP2_12S_transpose<-data.frame(seqtabNoC_WP2_12S_transpose)
i <- c(2:1224) 
seqtabNoC_WP2_12S_transpose[ , i] <- apply(seqtabNoC_WP2_12S_transpose[ , i], 2,  
                                           function(x) as.numeric(as.character(x)))

# Identify and drop samples with fewer than 1,000 total reads (too shallow for reliable community estimation)
rarefy_col_sum<-colSums(seqtabNoC_WP2_12S_transpose[c(2:1224)])
rarefy_col_sum<-data.frame(rarefy_col_sum)
rarefy_col_sum<-setDT(rarefy_col_sum, keep.rownames = TRUE)[]
rarefy_col_sum<-data.frame(rarefy_col_sum)

rarefy_col_sum<-rarefy_col_sum[!(rarefy_col_sum$rarefy_col_sum>=1000),]
drop<-rarefy_col_sum$rn

seqtabNoC_WP2_12S_transpose<-seqtabNoC_WP2_12S_transpose[,!(names(seqtabNoC_WP2_12S_transpose) %in% drop)]

# Join BLAST taxonomy onto the abundance table, using the raw sequence string as the key;
# this replaces literal sequences with species assignments
seqtabNoC_WP2_12S_tax <- BLAST_OUTPUT_12S %>% full_join(seqtabNoC_WP2_12S_transpose, by = c("Seq"))
seqtabNoC_WP2_12S_tax$rn <- NULL
seqtabNoC_WP2_12S_tax$X1 <- NULL
seqtabNoC_WP2_12S_tax$X3 <- NULL
seqtabNoC_WP2_12S_tax$X4 <- NULL
seqtabNoC_WP2_12S_tax$X5 <- NULL
seqtabNoC_WP2_12S_tax$Seq<-NULL
seqtabNoC_WP2_12S_tax<-na.omit(seqtabNoC_WP2_12S_tax)

i <- c(2:992) 
seqtabNoC_WP2_12S_tax[ , i] <- apply(seqtabNoC_WP2_12S_tax[ , i], 2,  
                                     function(x) as.numeric(as.character(x)))

# Calculate a per-sample 0.05% read threshold: ASVs below this proportion will be treated as noise
seqtabNoC_WP2_12S_tax_col_sum<-colSums(seqtabNoC_WP2_12S_tax[c(2:992)])
seqtabNoC_WP2_12S_tax_col_sum<-data.frame(seqtabNoC_WP2_12S_tax_col_sum)
seqtabNoC_WP2_12S_tax_col_sum$read_filter<-seqtabNoC_WP2_12S_tax_col_sum$seqtabNoC_WP2_12S_tax_col_sum *0.0005
seqtabNoC_WP2_12S_tax_col_sum<- seqtabNoC_WP2_12S_tax_col_sum %>% rownames_to_column("ID")

# Remove ASVs with no BLAST species assignment
seqtabNoC_WP2_12S_tax<-seqtabNoC_WP2_12S_tax %>% drop_na(species)
seqtabNoC_WP2_12S_tax$X1<-NULL

# Reshape to long format (one row per ASV × sample combination) for filtering and metadata merging
seqtabNoC_WP2_12S_tax_long <- gather(seqtabNoC_WP2_12S_tax, Sample, reads, 2:992, factor_key=TRUE)
seqtabNoC_WP2_12S_tax_long$reads<-   as.numeric(seqtabNoC_WP2_12S_tax_long$reads)

# Load project metadata and merge onto the long-format table
lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3_summary<-table(lofresh_metadata_WP2_3$Days,lofresh_metadata_WP2_3$SampleSite_code)
lofresh_metadata_WP2_3_summary<-data.frame(lofresh_metadata_WP2_3_summary)
write.csv(lofresh_metadata_WP2_3_summary,"lofresh_metadata_WP2_3_summary.csv")
seqtabNoC_WP2_12S_tax_long<-seqtabNoC_WP2_12S_tax_long %>% rename(ID=Sample)

seqtabNoC_WP2_12S_tax_long_meta <- merge(lofresh_metadata_WP2_3[, c("ID", "WP")],seqtabNoC_WP2_12S_tax_long, by="ID")
sum_reads<- sum(seqtabNoC_WP2_12S_tax_long_meta$reads)

# Attach per-sample read thresholds
seqtabNoC_WP2_12S_tax_long_meta_clean_filtered <- merge(seqtabNoC_WP2_12S_tax_col_sum[, c("ID", "read_filter")],seqtabNoC_WP2_12S_tax_long_meta, by="ID")

# Apply read depth filters: remove rows below the 0.05% per-sample threshold and below the 20-read hard floor
seqtabNoC_WP2_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2_12S_tax_long_meta_clean_filtered$reads<=seqtabNoC_WP2_12S_tax_long_meta_clean_filtered$read_filter),]
seqtabNoC_WP2_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2_12S_tax_long_meta_clean_filtered$reads<20),]

write.csv(seqtabNoC_WP2_12S_tax_long_meta_clean_filtered,"WP2_ASV_table_long_filtered_family.csv")
seqtabNoC_WP2_12S_tax_long_meta_clean_filtered <- read_csv("WP2_ASV_table_long_filtered_family.csv")

# Pivot to wide format (species × sample) for normalisation
seqtabNoC_WP2_12S_tax_wide_filtered<-seqtabNoC_WP2_12S_tax_long_meta_clean_filtered
seqtabNoC_WP2_12S_tax_wide_filtered$Order <- NULL 
seqtabNoC_WP2_12S_tax_wide_filtered$read_filter <- NULL 
seqtabNoC_WP2_12S_tax_wide_filtered$Seq_length <- NULL 
seqtabNoC_WP2_12S_tax_wide_filtered$WP <- NULL 

seqtabNoC_WP2_12S_tax_wide_filtered<-seqtabNoC_WP2_12S_tax_wide_filtered %>% 
  group_by(ID,species) %>% 
  summarize(reads=sum(reads)) %>%
  rename(ID=ID,species=species, reads= reads,)

seqtabNoC_WP2_12S_tax_wide_filtered <- spread(seqtabNoC_WP2_12S_tax_wide_filtered,ID,reads)
seqtabNoC_WP2_12S_tax_wide_filtered[is.na(seqtabNoC_WP2_12S_tax_wide_filtered)] <- 0
seqtabNoC_WP2_12S_tax_wide_filtered <- data.frame(t(seqtabNoC_WP2_12S_tax_wide_filtered))

names(seqtabNoC_WP2_12S_tax_wide_filtered) <- as.matrix(seqtabNoC_WP2_12S_tax_wide_filtered[1, ])
seqtabNoC_WP2_12S_tax_wide_filtered <- seqtabNoC_WP2_12S_tax_wide_filtered[-1, ]
seqtabNoC_WP2_12S_tax_wide_filtered[] <- lapply(seqtabNoC_WP2_12S_tax_wide_filtered, function(x) type.convert(as.character(x)))

# Normalise: divide each ASV's read count by the sample total to get relative proportions.
# NOTE: The column index range (1:68) may need to be updated if the number of species changes.
seqtabNoC_WP2_12S_tax_col_sum2<-rowSums(seqtabNoC_WP2_12S_tax_wide_filtered[c(1:68)])
seqtabNoC_WP2_12S_tax_col_sum2<-data.frame(seqtabNoC_WP2_12S_tax_col_sum2)
seqtabNoC_WP2_12S_tax_col_sum2$ID <- rownames(seqtabNoC_WP2_12S_tax_col_sum2)

seqtabNoC_WP2_12S_tax_normalised<- merge(seqtabNoC_WP2_12S_tax_col_sum2[, c("ID", "seqtabNoC_WP2_12S_tax_col_sum2")],seqtabNoC_WP2_12S_tax_long_meta_clean_filtered, by="ID")

seqtabNoC_WP2_12S_tax_normalised <- transform(seqtabNoC_WP2_12S_tax_normalised, normalised_reads = reads / seqtabNoC_WP2_12S_tax_col_sum2)

# Quick visualisation of normalised community composition across all samples
ggplot(seqtabNoC_WP2_12S_tax_normalised , aes(x = ID, fill = species, y = normalised_reads)) + 
  geom_bar(stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))+ theme(legend.position = "none")

# Pivot normalised data to wide format and save as the master normalised community matrix
seqtabNoC_WP2_12S_tax_normalised_wide<-seqtabNoC_WP2_12S_tax_normalised
seqtabNoC_WP2_12S_tax_normalised_wide$read_filter <- NULL 
seqtabNoC_WP2_12S_tax_normalised_wide$Days <- NULL 
seqtabNoC_WP2_12S_tax_normalised_wide$Date_sampled <- NULL 
seqtabNoC_WP2_12S_tax_normalised_wide$SampleSite_code <- NULL 
seqtabNoC_WP2_12S_tax_normalised_wide$SampleSite_time <- NULL 
seqtabNoC_WP2_12S_tax_normalised_wide$WP <- NULL 
seqtabNoC_WP2_12S_tax_normalised_wide$Seq_length <- NULL 
seqtabNoC_WP2_12S_tax_normalised_wide$reads <- NULL 
seqtabNoC_WP2_12S_tax_normalised_wide$seqtabNoC_WP2_12S_tax_col_sum2 <- NULL 

seqtabNoC_WP2_12S_tax_normalised_wide<-seqtabNoC_WP2_12S_tax_normalised_wide %>% 
  group_by(ID,species) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(ID=ID,species=species, normalised_reads= normalised_reads,)

seqtabNoC_WP2_12S_tax_normalised_wide <- spread(seqtabNoC_WP2_12S_tax_normalised_wide,ID,normalised_reads)
seqtabNoC_WP2_12S_tax_normalised_wide[is.na(seqtabNoC_WP2_12S_tax_normalised_wide)] <- 0
seqtabNoC_WP2_12S_tax_normalised_wide <- data.frame(t(seqtabNoC_WP2_12S_tax_normalised_wide))

names(seqtabNoC_WP2_12S_tax_normalised_wide) <- as.matrix(seqtabNoC_WP2_12S_tax_normalised_wide[1, ])
seqtabNoC_WP2_12S_tax_normalised_wide <- seqtabNoC_WP2_12S_tax_normalised_wide[-1, ]
seqtabNoC_WP2_12S_tax_normalised_wide[] <- lapply(seqtabNoC_WP2_12S_tax_normalised_wide, function(x) type.convert(as.character(x)))

write.csv(seqtabNoC_WP2_12S_tax_normalised_wide,"NORM_WP2_ASV_table_wide_filtered_family.csv")

# Split by work package and add further metadata fields; save each sub-dataset as its own CSV
### WP2A  Splitting work packages and adding metadata ###
lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
seqtabNoC_WP2A_12S_tax_long_meta_clean<-seqtabNoC_WP2_12S_tax_long_meta_clean_filtered %>% filter(WP == "WP2A")
seqtabNoC_WP2A_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "SampleSite_time")],seqtabNoC_WP2A_12S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2A_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "SampleSite_code")],seqtabNoC_WP2A_12S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2A_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "Date_sampled")],seqtabNoC_WP2A_12S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2A_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "Days")],seqtabNoC_WP2A_12S_tax_long_meta_clean, by="ID")

write.csv(seqtabNoC_WP2A_12S_tax_long_meta_clean,"WP2A_ASV_table_long_filtered_family.csv")
seqtabNoC_WP2A_12S_tax_long_meta_clean <- read_csv("WP2A_ASV_table_long_filtered_family.csv")

### WP2C Splitting work packages and adding metadata ###
seqtabNoC_WP2C_12S_tax_long_meta_clean<-seqtabNoC_WP2_12S_tax_long_meta_clean_filtered %>% filter(WP == "WP2C")
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
```

---

## 5. Stage 4 — Per-River Biological Curation & Community Visualisation

### Purpose
Each river dataset undergoes an identical block of steps:
1. **Replicate averaging** — technical replicates are averaged by `group_by(SampleSite_time, species)`.
2. **Species curation** — clearly non-native or human-food contamination species (e.g. Pacific cod, Atlantic herring, positive control mackerel) are removed by name. Some species are harmonised to preferred synonyms.
3. **Negative control removal** — known negative control samples are excluded by sample ID.
4. **Wide-format conversion and normalisation** — long-format data is pivoted to a species × site matrix; reads are normalised to proportions.
5. **Visualisation** — stacked bar charts of relative and absolute amplicon abundance; top-taxon plots ordered by cumulative reads.
6. **Alpha diversity** — Shannon index, Simpson index, and species richness are calculated and saved.
7. **PERMANOVA** — Bray-Curtis PERMANOVA tests which environmental variables explain community composition.
8. **NMDS** — Ordination of Bray-Curtis distances; scores merged with metadata for plotting.

The blocks below are shown for each river in full. The logic is identical across all five; river-specific species exclusions and column index ranges are noted where they differ.

---

### WP2A (Conwy)

```r
#________________________________________________________
########  WP2A BLAST TAXONOMY ANALYSIS  ########
#________________________________________________________

library(tidyverse)
library(readr)
library(data.table)
library(vegan)
library(plotly)
library(microDecon)
library(RColorBrewer)
library(emmeans)
library(mgcv)
library(ggplot2)
library(ggpubr)
setwd("WP2_analysis/12S_BLAST")
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered <- read_csv("WP2A_ASV_table_long_filtered_family.csv")

# Average across technical replicates at the SampleSite_time × species level
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% 
 group_by(SampleSite_time,species) %>% 
  summarize(reads=mean(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads)

seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% 
group_by(SampleSite_time,species) %>% 
summarize(reads=mean(reads)) %>%
rename(SampleSite_time=SampleSite_time,species=species, reads= reads)

# Remove non-native/human-food species unlikely to represent wild freshwater fish communities;
# these are probable contaminants from processed fish products or positive controls
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Gadus_macrocephalus"),]#Pacific cod
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Hippoglossus_stenolepis"),]#Pacific halibut
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Ammodytes_personatus"),]#Pacific Sand Lance 
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Clupea_harengus"),]#Atlantic herring 
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Clupea_pallasii"),]#Pacific herring 
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Copadichromis_virginalis"),]#Haplochromine cichlid  
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Favonigobius_gymnauchen"),]#Sand goby   
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Pholis_crassispina"),]#Rock gunnel  
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Melanogrammus_aeglefinus"),]#Haddock
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Platichthys_stellatus"),]#Starry flounder
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Scomber_scombrus"),]#Atlantic mackerel - POSITIVE CONTROL

# Harmonise species names: reassign ambiguous/American species to their European equivalents
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species<-sub('Salmo_obtusirostris','Salmo_trutta',seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species)
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species<-sub('Anguilla_rostrata','Anguilla_anguilla',seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species)
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species<-sub('Salvelinus_fontinalis_x_Salvelinus_malma','Salvelinus_malma',seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species)

# Remove known negative control samples (blank extractions that should contain no fish DNA)
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$SampleSite_time=="ENC_T10"),]
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$SampleSite_time=="ENC_T17"),]
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$SampleSite_time=="ENC_T18"),]
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$SampleSite_time=="ENC_T19"),]

# Pivot to wide format (species × site) for normalisation and ordination
seqtabNoC_WP2A_12S_tax_wide_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered
seqtabNoC_WP2A_12S_tax_wide_filtered$read_filter <- NULL 
seqtabNoC_WP2A_12S_tax_wide_filtered$Days <- NULL 
seqtabNoC_WP2A_12S_tax_wide_filtered$Date_sampled <- NULL 
seqtabNoC_WP2A_12S_tax_wide_filtered$SampleSite_code <- NULL 
seqtabNoC_WP2A_12S_tax_wide_filtered$WP <- NULL 
seqtabNoC_WP2A_12S_tax_wide_filtered$Seq_length <- NULL 

seqtabNoC_WP2A_12S_tax_wide_filtered<-seqtabNoC_WP2A_12S_tax_wide_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=sum(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads,)

seqtabNoC_WP2A_12S_tax_wide_filtered <- spread(seqtabNoC_WP2A_12S_tax_wide_filtered,SampleSite_time,reads)
seqtabNoC_WP2A_12S_tax_wide_filtered[is.na(seqtabNoC_WP2A_12S_tax_wide_filtered)] <- 0
seqtabNoC_WP2A_12S_tax_wide_filtered <- data.frame(t(seqtabNoC_WP2A_12S_tax_wide_filtered))

names(seqtabNoC_WP2A_12S_tax_wide_filtered) <- as.matrix(seqtabNoC_WP2A_12S_tax_wide_filtered[1, ])
seqtabNoC_WP2A_12S_tax_wide_filtered <- seqtabNoC_WP2A_12S_tax_wide_filtered[-1, ]
seqtabNoC_WP2A_12S_tax_wide_filtered[] <- lapply(seqtabNoC_WP2A_12S_tax_wide_filtered, function(x) type.convert(as.character(x)))

write.csv(seqtabNoC_WP2A_12S_tax_wide_filtered,"WP2A_ASV_table_wide_filtered.csv")

# Relative amplicon frequency abundance plot (proportional stacked bar per site)
ggplot(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered, aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))

# Normalise reads by sample total; also calculate per-species total reads across all samples (for ordering plots)
# NOTE: column range (1:16) reflects the number of species in WP2A — update if this changes
seqtabNoC_WP2A_12S_tax_col_sum2<-rowSums(seqtabNoC_WP2A_12S_tax_wide_filtered[c(1:16)])
WP2A_12S_sps_read_summary<-colSums(seqtabNoC_WP2A_12S_tax_wide_filtered[c(1:16)])
WP2A_12S_sps_read_summary<-data.frame(WP2A_12S_sps_read_summary)

seqtabNoC_WP2A_12S_tax_normalised<- merge(seqtabNoC_WP2A_12S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2A_12S_tax_col_sum2")],seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered, by="SampleSite_time")
seqtabNoC_WP2A_12S_tax_top<- merge(seqtabNoC_WP2A_12S_tax_col_sum3[, c("species", "seqtabNoC_WP2A_12S_tax_col_sum3")],seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered, by="species")

seqtabNoC_WP2A_12S_tax_normalised <- transform(seqtabNoC_WP2A_12S_tax_normalised, normalised_reads = reads / seqtabNoC_WP2A_12S_tax_col_sum2)

# Normalised stacked bar plot
ggplot(seqtabNoC_WP2A_12S_tax_normalised , aes(x = SampleSite_time, fill = species, y = normalised_reads)) + 
  geom_bar(stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")

# Add metadata fields (Days, SampleSite_code) for seasonal and spatial breakdowns
lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3_summary = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

seqtabNoC_WP2A_12S_tax_normalised <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "Days")],seqtabNoC_WP2A_12S_tax_normalised, by="SampleSite_time")
seqtabNoC_WP2A_12S_tax_normalised <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "SampleSite_code")],seqtabNoC_WP2A_12S_tax_normalised, by="SampleSite_time")

# Count per-season detections for each species and plot as a faceted line chart
sps_seasonal_detection<-unique(seqtabNoC_WP2A_12S_tax_normalised[c("SampleSite_time", "species")])
sps_seasonal_detection <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "season")],sps_seasonal_detection, by="SampleSite_time")

sps_seasonal_detection_summary<-count(sps_seasonal_detection, season, species)
sps_seasonal_detection_summary$season <- factor(sps_seasonal_detection_summary$season, levels = c("Winter", "Spring", "Summer", "Autumn"))

ggplot(data=sps_seasonal_detection_summary, aes(x=season, y=n)) +
  geom_line()+
  geom_point()+facet_grid((.~species))

# Top-taxon stacked bar plot ordered by cumulative read abundance, faceted by sample site
seqtabNoC_WP2A_12S_tax_top <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "SampleSite_code")],seqtabNoC_WP2A_12S_tax_top, by="SampleSite_time")
seqtabNoC_WP2A_12S_tax_top %>%
  mutate(species = fct_reorder(species, desc(seqtabNoC_WP2A_12S_tax_col_sum3))) %>%
  ggplot(aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  theme(axis.text.x = element_text(angle = 90))+scale_fill_brewer(palette="Paired")+facet_wrap(~ SampleSite_code,scales="free",ncol = 15)

write.csv(seqtabNoC_WP2A_12S_tax_normalised,"NORM_WP2A_ASV_table_LONG_filtered.csv")

# Pivot normalised data to wide format and save
seqtabNoC_WP2A_12S_tax_normalised_wide<-seqtabNoC_WP2A_12S_tax_normalised
seqtabNoC_WP2A_12S_tax_normalised_wide$read_filter <- NULL 
seqtabNoC_WP2A_12S_tax_normalised_wide$Days <- NULL 
seqtabNoC_WP2A_12S_tax_normalised_wide$Date_sampled <- NULL 
seqtabNoC_WP2A_12S_tax_normalised_wide$SampleSite_code <- NULL 
seqtabNoC_WP2A_12S_tax_normalised_wide$WP <- NULL 
seqtabNoC_WP2A_12S_tax_normalised_wide$Seq_length <- NULL 
seqtabNoC_WP2A_12S_tax_normalised_wide$reads <- NULL 
seqtabNoC_WP2A_12S_tax_normalised_wide$seqtabNoC_WP2_12S_tax_col_sum2 <- NULL 

seqtabNoC_WP2A_12S_tax_normalised_wide<-seqtabNoC_WP2A_12S_tax_normalised_wide %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, normalised_reads= normalised_reads,)

seqtabNoC_WP2A_12S_tax_normalised_wide <- spread(seqtabNoC_WP2A_12S_tax_normalised_wide,SampleSite_time,normalised_reads)
seqtabNoC_WP2A_12S_tax_normalised_wide[is.na(seqtabNoC_WP2A_12S_tax_normalised_wide)] <- 0
seqtabNoC_WP2A_12S_tax_normalised_wide <- data.frame(t(seqtabNoC_WP2A_12S_tax_normalised_wide))

names(seqtabNoC_WP2A_12S_tax_normalised_wide) <- as.matrix(seqtabNoC_WP2A_12S_tax_normalised_wide[1, ])
seqtabNoC_WP2A_12S_tax_normalised_wide <- seqtabNoC_WP2A_12S_tax_normalised_wide[-1, ]
seqtabNoC_WP2A_12S_tax_normalised_wide[] <- lapply(seqtabNoC_WP2A_12S_tax_normalised_wide, function(x) type.convert(as.character(x)))

write.csv(seqtabNoC_WP2A_12S_tax_normalised_wide,"NORM_WP2A_ASV_table_wide_filtered.csv")

# PERMANOVA: test which environmental variables explain community composition (Bray-Curtis, 999 permutations)
NORM_WP2A_ASV_table_wide_filtered <- read_csv("NORM_WP2A_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2A_ASV_table_wide_filtered
adonis_data_samples<-adonis_data$...1
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(ID = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = ...1)

lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$ID),]

adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","river_section","season","Days","Dist_from_lake","pH",
                                               "conductivity_uS_cm","gran_alk_ueq_L","daily_flow_by_catchment","monthly_flow_by_catchment","seasonal_flow_by_catchment",
                                               "week_before_flow_by_catchment","river_gradient_percent","rainfall_monthly_average",
                                               "temp_monthly_average"
                                           )],adonis_data, by="SampleSite_time")

adonis_data<-na.omit(adonis_data)
adonis_data_meta<- adonis_data[c(1:16)]
adonis_data<-adonis_data[,-c(1:16)]

adonis_data_meta$Days<-as.numeric(adonis_data_meta$Days)

# Model: season (nested within days) + spatial/environmental predictors, tested by margin
adon.results_WP2<-adonis2(adonis_data ~ adonis_data_meta$season/adonis_data_meta$Days+
                            adonis_data_meta$Dist_from_lake+
                            adonis_data_meta$daily_flow_by_catchment+
                            adonis_data_meta$monthly_flow_by_catchment+
                            adonis_data_meta$seasonal_flow_by_catchment+
                            adonis_data_meta$week_before_flow_by_catchment+
                            adonis_data_meta$river_gradient_percent+
                            adonis_data_meta$rainfall_monthly_average+
                            adonis_data_meta$temp_monthly_average+
                            adonis_data_meta$gran_alk_ueq_L+
                            adonis_data_meta$conductivity_uS_cm+
                            adonis_data_meta$pH,
                          method="bray",perm=999,by="margin")
print(adon.results_WP2)

# NMDS ordination (3 dimensions, Bray-Curtis, up to 1000 iterations)
abund_table_12S_ID <- read.csv("NORM_WP2A_ASV_table_wide_filtered.csv", header = TRUE)
abund_table_12S <- abund_table_12S_ID[,-1]
com <- abund_table_12S[,col(abund_table_12S)]
m_com <- as.matrix(com)
set.seed(123)
nmds = metaMDS(m_com, distance = "bray", try = 1000, trymax = 1000, k = 3)

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds))
#add columns to data frame 
data.scores$SampleSite_time = abund_table_12S_ID[,1]

lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]
data.scores_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","SampleSite_code","Days","Dist_from_lake","season")],data.scores, by="SampleSite_time")

# Plot NMDS: points shaped by season, coloured by distance from lake; red centroids with spider lines
gg <- merge(data.scores_meta_data,aggregate(cbind(NMDS1.x=NMDS1,NMDS2.y=NMDS2)~season,data.scores_meta_data,mean),by="season")
ggplot(gg, aes(x = NMDS1, y = NMDS2,colour = Dist_from_lake)) + scale_shape_manual(values=c(15,16,17,8,18))+
  geom_point(size = 4, aes(shape=season,colour = Dist_from_lake))+
  theme_bw()+  theme(axis.title.x=element_blank())+  theme(axis.title.y=element_blank())+ 
  geom_segment(aes(x=NMDS1.x, y=NMDS2.y, xend=NMDS1, yend=NMDS2))+
  geom_point(aes(x=NMDS1.x,y=NMDS2.y),size=3,colour="red")
```

---

### Alpha diversity modelling — WP2A

```r
########  ALPHA DIVERSITY METRICS  ########

# Calculate Shannon entropy, Simpson index, and observed species richness per sample
shannon_index<-diversity(seqtabNoC_WP2A_12S_tax_normalised_wide, index = "shannon")
shannon_index<-data.frame(shannon_index)
ID <- rownames(shannon_index)
shannon_index <- cbind(shannon_index,ID)

simpson_index<-diversity(seqtabNoC_WP2A_12S_tax_normalised_wide, index = "simpson")
simpson_index<-data.frame(simpson_index)
ID <- rownames(simpson_index)
simpson_index <- cbind(simpson_index,ID)

species_count <- rowSums(seqtabNoC_WP2A_12S_tax_normalised_wide >0)
species_count<-data.frame(species_count)
species_count$ID <- rownames(species_count)

alpha_12S_WP2A<-merge(shannon_index,simpson_index,by="ID")
alpha_12S_WP2A<-merge(alpha_12S_WP2A,species_count,by="ID")

write.csv(alpha_12S_WP2A,"alpha_12S_WP2A.csv")

########  WP2A ALPHA DIVERSITY MODELLING  ##########
library(mgcv)

# Merge species richness with a full suite of environmental covariates for modelling
species_count_dataset<-species_count %>% rename(SampleSite_time = ID)
species_count_dataset<- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time","SampleSite_code",
                                               "season",
                                               "Days",
                                               "Dist_from_lake",
                                               "temp_monthly_average",
                                               "conductivity_uS_cm",
                                               "gran_alk_ueq_L","daily_flow_by_catchment",
                                               "monthly_flow_by_catchment",
                                               "seasonal_flow_by_catchment",
                                               "week_before_flow_by_catchment",
                                               "river_gradient_percent",
                                               "rainfall_monthly_average",
                                               "pH","river_section")],species_count_dataset, by="SampleSite_time")

species_count_dataset$Days<-as.numeric(species_count_dataset$Days)
species_count_dataset<-species_count_dataset[species_count_dataset$species_count != 0, ] 
species_count_dataset<-na.omit(species_count_dataset)

# Exploratory plots: species richness vs distance-from-lake (by season) and vs rainfall
ggplot(species_count_dataset, aes(Dist_from_lake, species_count,colour=season)) + 
  geom_point() + 
  geom_smooth(method = lm, se = TRUE)+theme_bw()
ggplot(species_count_dataset, aes(rainfall_monthly_average, species_count)) + 
  geom_point() + 
  geom_smooth(se = TRUE)+theme_bw()
ggplot(shannon_index_dataset, aes(Dist_from_lake, shannon_index,colour=season)) + 
  geom_point() + 
  geom_smooth(method = lm, se = TRUE)+theme_bw()

hist(species_count_dataset$species_count)

# Linear model: species richness ~ season × distance + environmental predictors
lm_12S_species_count = lm(species_count~season*
                            Dist_from_lake+
                            temp_monthly_average+
                            conductivity_uS_cm+
                            daily_flow_by_catchment+
                            monthly_flow_by_catchment+
                            week_before_flow_by_catchment+
                            river_gradient_percent+
                            rainfall_monthly_average+
                            pH, data = species_count_dataset)
anova(lm_12S_species_count)
```

---

### Sparse partial least squares analysis (sPLS) — WP2A

```r
######## SPARSE PARTIAL LEAST SQUARES ANALYSIS ##########
# sPLS finds latent components that maximally co-vary between species composition (X) and
# environmental variables (Y), identifying which species are most associated with each predictor.
library(mixOmics)

NORM_WP2A_ASV_table_wide_filtered <- read_csv("NORM_WP2A_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2A_ASV_table_wide_filtered

adonis_data_samples<-adonis_data$X1
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(ID = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = X1)

lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time",
                                               "Days",
                                               "Dist_from_lake",
                                               "temp_monthly_average",
                                               "conductivity_uS_cm",
                                               "gran_alk_ueq_L","daily_flow_by_catchment",
                                               "monthly_flow_by_catchment",
                                               "seasonal_flow_by_catchment",
                                               "week_before_flow_by_catchment",
                                               "river_gradient_percent",
                                               "rainfall_monthly_average",
                                               "pH")],adonis_data, by="SampleSite_time")

adonis_data<-na.omit(adonis_data)
adonis_data_meta<- adonis_data[c(1:29)]
adonis_data<-adonis_data[,-c(1:29)]
adonis_data_meta$Days<-as.numeric(adonis_data_meta$Days)
adonis_data_meta$SampleSite_time<-NULL

ncomp = 1
result.spls <- spls(adonis_data, adonis_data_meta, ncomp = ncomp,  mode = 'canonical')

# Cross-validate model performance with 5-fold CV repeated 99 times
design <- data.frame(samp = adonis_data_meta$sample)
tune.spls <- perf(result.spls, validation = 'Mfold', folds = 5,
                  criterion = 'all',nrepeat = 99, progressBar = TRUE)

tune.spls$R2
plot(tune.spls$Q2.total)
abline(h = 0.0975)

plot(tune.spls, criterion="R2", type = 'l')
plot(tune.spls, criterion="Q2", type = 'l')
plot(tune.spls)

# Clustered image map: visualises species–environment co-inertia structure
cim_plot<-cim(result.spls,
              comp = 1:1,
              margins = c(15, 15))
```

---

### Mantel test — WP2A

```r
######## MANTEL TEST ##########
# Tests whether pairwise community dissimilarity correlates with pairwise environmental distance
# (Spearman-based, 9999 permutations)
adonis_data_meta$Days<-NULL
adonis_data_meta$Dist_from_lake<-NULL

dist.abund = vegdist(adonis_data, method = "bray")
dist.meta = dist(adonis_data_meta, method = "euclidean")

abund_meta = mantel(dist.abund, dist.meta, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_meta

# Scatter plot of pairwise Bray-Curtis dissimilarity vs. pairwise environmental distance
aa = as.vector(dist.abund)
mm = as.vector(dist.meta)
mat = data.frame(aa,mm)

ggplot(mat, aes(y = aa, x = mm)) + 
  geom_point(size = 1) + 
  theme_bw()+ geom_smooth(method = lm)+ scale_y_continuous(limits = c(0, 1))
```

---

## 6. Stage 5 — Alpha Diversity: Cross-River Comparison

### Purpose
Alpha diversity outputs from all rivers are combined and filtered to the shared sampling window (the Conwy-only dates not coordinated with the other rivers are removed). Species richness is plotted against distance from lake coloured by river, and a linear model tests whether the richness–distance relationship differs among rivers.

```r
#________________________________________________________
########  COMPILING WP2 CURATED NORMALISED OTU TABLES FOR ALPHA AND BETA DIVERSITY  #####
# (WP2B excluded)
#________________________________________________________

# Load per-river alpha diversity tables (note: Skan commented out for combined analysis)
alpha_12S_WP2A <- read_csv("alpha_12S_WP2A.csv")
alpha_12S_WP2C_Glatt <- read_csv("alpha_12S_WP2C_Glatt.csv")
alpha_12S_WP2C_Gwash <- read_csv("alpha_12S_WP2C_Gwash.csv")
alpha_12S_WP2C_Towy <- read_csv("alpha_12S_WP2C_Towy.csv")
#alpha_12S_WP2C_Skan <- read_csv("alpha_12S_WP2C_Skan.csv")

# Combine all rivers and remove samples with zero Shannon diversity (likely empty/failed samples)
alpha_12S_WP2<-rbind(alpha_12S_WP2A,alpha_12S_WP2C_Glatt,alpha_12S_WP2C_Gwash,alpha_12S_WP2C_Towy)#,alpha_12S_WP2C_Skan)
alpha_12S_WP2<-alpha_12S_WP2[!(alpha_12S_WP2$shannon_index==0),]

# Add spatial and environmental metadata
lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3_summary = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]  

alpha_12S_WP2<-alpha_12S_WP2 %>% rename(SampleSite_time = ID)

alpha_12S_WP2 <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "Dist_from_lake")],alpha_12S_WP2, by="SampleSite_time")
alpha_12S_WP2 <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "River")],alpha_12S_WP2, by="SampleSite_time")
alpha_12S_WP2 <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "Date_sampled")],alpha_12S_WP2, by="SampleSite_time")
alpha_12S_WP2 <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "pH")],alpha_12S_WP2, by="SampleSite_time")
alpha_12S_WP2 <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "conductivity_uS_cm")],alpha_12S_WP2, by="SampleSite_time")

# Remove Conwy-only sampling dates that were not coordinated with the comparative river survey
# (keeping only the shared summer sampling window used by all rivers)
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "01.09.17" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "02.11.17" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "04.01.18" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "11.04.18" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "12.06.17" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "12.08.17" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "12.10.17" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "14.02.18" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "14.12.17" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "18.04.18" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "18.05.17" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "22.09.17" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "23.11.17" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "25.01.18" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "27.04.17" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "28.03.18" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "29.06.17" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "29.06.17" , ] 

# Plot species richness vs distance-from-lake, one trend line per river
ggplot(alpha_12S_WP2, aes(x=Dist_from_lake, y=species_count, color=River)) +
  geom_point() + 
  geom_smooth(method = "lm",fullrange = T)+theme_bw()+scale_color_manual(values=c("#56B4E9","#F0E442","#E69F00","#0072B2","#D55E00","#999999"))

hist(alpha_12S_WP2$shannon_index)

# Linear model: does the richness–distance relationship differ by river?
lm_12S_WP2.1 = lm(species_count~Dist_from_lake*River, data = alpha_12S_WP2)
anova(lm_12S_WP2.1)
```

---

## 7. Stage 6 — Beta Diversity: Cross-River PERMANOVA and NMDS

### Purpose
The normalised ASV tables from all five rivers are pooled and pivoted to wide format. After applying the same date filter used in the alpha diversity comparison, a Bray-Curtis PERMANOVA tests the effects of River and distance-from-lake on community composition. NMDS ordination is then plotted with centroid spider lines coloured by distance from lake and shaped by river.

```r
#BETA DIVERSITY - PERMANOVA ON TAXON
setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/12S_BLAST")

# Load per-river normalised long-format tables
NORM_WP2A_ASV_table_LONG_filtered <- read_csv("NORM_WP2A_ASV_table_LONG_filtered.csv")
NORM_WP2C_Glatt_ASV_table_LONG_filtered <- read_csv("NORM_WP2C_Glatt_ASV_table_LONG_filtered.csv")
NORM_WP2C_Gwash_ASV_table_LONG_filtered <- read_csv("NORM_WP2C_Gwash_ASV_table_LONG_filtered.csv")
NORM_WP2C_Towy_ASV_table_LONG_filtered <- read_csv("NORM_WP2C_Towy_ASV_table_LONG_filtered.csv")
NORM_WP2C_Skan_ASV_table_LONG_filtered <- read_csv("NORM_WP2C_Skan_ASV_table_LONG_filtered.csv")

# Drop extraneous columns before binding
NORM_WP2A_ASV_table_LONG_filtered$Days<-NULL
NORM_WP2A_ASV_table_LONG_filtered$SampleSite_code<-NULL

NORM_WP2A_ASV_table_LONG_filtered[,3] <- NULL
NORM_WP2C_Glatt_ASV_table_LONG_filtered[,3] <- NULL
NORM_WP2C_Gwash_ASV_table_LONG_filtered[,3] <- NULL
NORM_WP2C_Towy_ASV_table_LONG_filtered[,3] <- NULL
NORM_WP2C_Skan_ASV_table_LONG_filtered[,3] <- NULL

# Combine all rivers into a single long-format table
NORM_WP2_ASV_table_LONG_filtered<-rbind(NORM_WP2A_ASV_table_LONG_filtered,NORM_WP2C_Glatt_ASV_table_LONG_filtered,
                                        NORM_WP2C_Gwash_ASV_table_LONG_filtered ,NORM_WP2C_Towy_ASV_table_LONG_filtered,NORM_WP2C_Skan_ASV_table_LONG_filtered)

# Pivot to wide format (species × sample) for ordination
seqtabNoC_WP2A_12S_tax_normalised_wide<-NORM_WP2_ASV_table_LONG_filtered
seqtabNoC_WP2A_12S_tax_normalised_wide$X1<- NULL 
seqtabNoC_WP2A_12S_tax_normalised_wide$reads <- NULL 

seqtabNoC_WP2A_12S_tax_normalised_wide<-seqtabNoC_WP2A_12S_tax_normalised_wide %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, normalised_reads= normalised_reads,)

seqtabNoC_WP2A_12S_tax_normalised_wide <- spread(seqtabNoC_WP2A_12S_tax_normalised_wide,SampleSite_time,normalised_reads)
seqtabNoC_WP2A_12S_tax_normalised_wide[is.na(seqtabNoC_WP2A_12S_tax_normalised_wide)] <- 0
seqtabNoC_WP2A_12S_tax_normalised_wide <- data.frame(t(seqtabNoC_WP2A_12S_tax_normalised_wide))

names(seqtabNoC_WP2A_12S_tax_normalised_wide) <- as.matrix(seqtabNoC_WP2A_12S_tax_normalised_wide[1, ])
seqtabNoC_WP2A_12S_tax_normalised_wide <- seqtabNoC_WP2A_12S_tax_normalised_wide[-1, ]
seqtabNoC_WP2A_12S_tax_normalised_wide[] <- lapply(seqtabNoC_WP2A_12S_tax_normalised_wide, function(x) type.convert(as.character(x)))

adonis_data<-seqtabNoC_WP2A_12S_tax_normalised_wide
adonis_data<-setDT(adonis_data, keep.rownames = TRUE)[]
adonis_data<-data.frame(adonis_data)
adonis_data_samples<-adonis_data$rn
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(ID = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = rn)

lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]
adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","River","Date_sampled","Dist_from_lake")],adonis_data, by="SampleSite_time")

# Apply the same date filter as the alpha diversity comparison to restrict to the shared sampling window
adonis_data<-adonis_data[adonis_data$Date_sampled != "01.09.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "02.11.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "04.01.18" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "11.04.18" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "12.06.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "12.08.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "12.10.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "14.02.18" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "14.12.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "18.04.18" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "18.05.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "22.09.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "23.11.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "25.01.18" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "27.04.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "28.03.18" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "29.06.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "29.06.17" , ] 

adonis_data<-na.omit(adonis_data)
adonis_data_meta<- adonis_data[c(1:5)]
adonis_data<-adonis_data[,-c(1:5)]

# Cross-river PERMANOVA: River × distance-from-lake interaction (Bray-Curtis, 999 permutations)
adon.results_WP2<-adonis2(adonis_data ~ adonis_data_meta$River*adonis_data_meta$Dist_from_lake,
                         method="bray",perm=999)
print(adon.results_WP2)

# Cross-river NMDS ordination (Bray-Curtis, 3D)
abund_table_12S_ID <- read.csv("NORM_WP2C_Skan_ASV_table_wide_filtered.csv", header = TRUE)
abund_table_12S <- abund_table_12S_ID[,-1]

com <- adonis_data[,col(adonis_data)]
m_com <- as.matrix(com)
set.seed(123)
nmds = metaMDS(m_com, distance = "bray", try = 1000, trymax = 1000, k = 3)

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds))
#add columns to data frame 
data.scores$SampleSite_time = adonis_data_meta[,1]

lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]
data.scores_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","Dist_from_lake","River")],data.scores, by="SampleSite_time")

# Plot: points shaped by river, coloured by distance from lake; red centroids per river with spider lines
gg <- merge(data.scores_meta_data,aggregate(cbind(NMDS1.x=NMDS1,NMDS2.y=NMDS2)~River,data.scores_meta_data,mean),by="River")
ggplot(gg, aes(x = NMDS1, y = NMDS2,colour = Dist_from_lake)) + scale_shape_manual(values=c(15,16,17,8,18))+
  geom_point(size = 4, aes(shape=River,colour = Dist_from_lake))+
  theme_bw()+  theme(axis.title.x=element_blank())+  theme(axis.title.y=element_blank())+ 
  geom_point(aes(x=NMDS1.x,y=NMDS2.y),size=3,colour="red")+
  geom_segment(aes(x=NMDS1.x, y=NMDS2.y, xend=NMDS1, yend=NMDS2))
