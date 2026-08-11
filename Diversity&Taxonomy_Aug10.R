# Diversity Analysis MC

# Questions for Sophie
# - Using BLAST for initial annotations to filter out contaminant/commensal species. Should we use the ASV/Bayes first?
# - At which point is it most appropriate to filter out low ASV counts?
# - Advice on comparisons between BLAST and VSEARCH/Bayes? 
# - Best way to identify endangered species

# Required files to run
# - MowrowCreek_ASV_table_long_filtered.csv
# - table_full from DADA2
# - Antonios List
# - Sophies List
# - 12SV5_blast_full_filtered
# - bayesian_taxonomy_Sadler_12SV5_final
# - JackASVs_12SV5_vsearch_bayesian
# - repseqs_12SV5
# - 
 
#====PACKAGES====
library(tidyverse)
library(sessioninfo)
library(dplyr)

# eDNA / metabarcoding
library(ranacapa)
library(phyloseq)
library(microDecon)
library(ngsLCA)
library(iNEXT)

# taxonomy
library(taxonomizr)
library(taxizedb)

# ecology
library(vegan)
library(hillR)
library(MeanRarity)

# statistics
library(gam)
library(limma)

# phylogeny
library(ape)
library(data.tree)
library(itol.toolkit)

# visualization
library(plotly)
library(ggh4x)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(wesanderson)
library(pheatmap)

# overlap plots
library(VennDiagram)
library(eulerr)
library(ggvenn)

# high-performance data handling
library(data.table)
library(reshape2)

#====Workspace Setup====
# Let's make a new folder for this run. And a new workspace
setwd("/Users/jacksadler/Desktop/project_thesis/taxonomy/12SV5_BLAST/Full_Outputs_MowrowCreekBLAST/Condensed_Code_June_Tweaks/")
#Import

# ASV occurrence setup from DADA2 table (converted qza)
seqtab_12SV5 <- table_full
seqtab_12SV5<-as.data.frame(seqtab_12SV5)
seqtab_12SV5<-seqtab_12SV5 %>% rename(fasta_id = OTU.ID)
# Convert read counts to numeric
seqtab_12SV5<-data.frame(seqtab_12SV5)
i <- c(2:208) 
seqtab_12SV5[ , i] <- apply(seqtab_12SV5[ , i], 2, function(x) as.numeric(as.character(x)))

#====SAMPLE METADATA LOOKUPS====
# one row per sample_id
meta_lookup_MC <- seqtab_MC_12SV5_tax_long_meta_clean_filtered %>%
  select(sample_id, SubLocation, Control) %>%
  distinct() %>%
  group_by(sample_id) %>%
  slice(1) %>%
  ungroup()

#====RepSeqs Processing====

##----Calculate the character length of each unique ASV sequence----
repseqs_12SV5_length<-repseqs_12SV5
repseqs_12SV5_length$seq_length <- nchar(repseqs_12SV5$sequence)
hist(repseqs_12SV5_length$seq_length)

##----20% ASV trim----
# ----Remove ASVs >20% longer than the expected Riaz 12SV5 amplicon size (110 bp → upper limit 132 bp)
repseqs_12SV5_length<-repseqs_12SV5_length[!(repseqs_12SV5_length$seq_length>132),]

# ----Remove ASVs >20% shorter than the expected Riaz 12SV5 amplicon size (110 bp → lower limit 88 bp)
repseqs_12SV5_length<-repseqs_12SV5_length[!(repseqs_12SV5_length$seq_length<88),]
hist(repseqs_12SV5_length$seq_length)

##----Save lookup table----
#links SEQ IDs back to their original sequence lengths (needed for joining later)
write.table(repseqs_12SV5_length,"~/full_12SV5_length_filtered_SEQ_to_sequence_lookup_table.txt",quote = FALSE)

##----Format each entry as a two-line FASTA record----
# (header\nsequence) by replacing the space separator with a newline
repseqs_12SV5_length$fasta_entry <- paste(repseqs_12SV5_length$fasta_id, 
                                          repseqs_12SV5_length$sequence, 
                                          sep = "\n")

##----Retain a standalone lookup dataframe----
# (FASTA ID ↔ sequence length) for use after BLAST
length_adjust_fasta_id<-repseqs_12SV5_length$fasta_id
length_adjust_fasta_id<-data.frame(length_adjust_fasta_id)
length_adjust_fasta_id$seq_length<-repseqs_12SV5_length$seq_length

repseqs_12SV5_length_min <- repseqs_12SV5_length
repseqs_12SV5_length_min$fasta_id<-NULL
repseqs_12SV5_length_min$seq_length<-NULL

##----Export the length-filtered sequences as a FASTA-formatted CSV----
# this is the BLAST input file
#write.csv(repseqs_12SV5_length_min$fasta_entry,"R_outputs/full_12SV5_length_filtered_repseqs_12SV5_BLAST_INPUT.fasta", row.names=FALSE,quote = FALSE)
writeLines(
  paste0(">", repseqs_12SV5_length$fasta_id, "\n", repseqs_12SV5_length$sequence),
  con = "~/full_12SV5_length_filtered_repseqs_12SV5_BLAST_INPUT.fasta"
)


#BLAST taxonomy setup
seqtab_MC_12SV5_tax_long_meta_clean_filtered <- MowrowCreek_ASV_table_long_filtered

#====Make taxonomy-independent file for MicroDecon====
# I need to make a df with each FASTA_ID per ASV for each sample site. Then merge
# with Sophie's assignments. Also need to clarify the number of unique ASVs vs read counts
# Use seqtab_12SV5

# merge seqtab_12SV5 with repseqs_12SV5_length. 
# We want the trimmed version that has +/- 20% of our reads outside bp length removed

repseqs_12SV5_length <- repseqs_12SV5_length %>% 
  rename(fasta_id = header)

repseqs_12SV5_filtered <- repseqs_12SV5_length %>%
 select(-seq_length, -fasta_entry, -sequence)

# pivot seqtab of DADA2 outputs to long format
# NOTE: This long version of seqtab_12SV5 is what you'll use for all future analyses too

seqtab_12SV5_long <- seqtab_12SV5 %>%
  pivot_longer(
    -c(fasta_id),
    names_to = "sample_id",
    values_to = "reads"
  ) %>%
  mutate(reads = as.numeric(reads))

seqtab_12SV5_long <- seqtab_12SV5_long %>%
  mutate(sample_id = gsub("\\.", "-", sample_id)) %>%
  group_by(sample_id)

# merge with MC metadata to remove any non-MC reads and remove columns we don't need for MicroDecon

seqtab_12SV5_long_MC <- left_join(seqtab_12SV5_long, meta_lookup_MC, by = "sample_id")
seqtab_12SV5_long_MC <- seqtab_12SV5_long_MC %>%
  drop_na(Control)
seqtab_12SV5_long_MC <- seqtab_12SV5_long_MC %>%
  select(-c(Village, SubLocation, Site, Replicate, pH, Temp, EC))

# Now merge with repseqs_12SV5_filtered to remove any ASVs that we excluded based on size constraints

seqtab_12SV5_long_MC_filtered <-seqtab_12SV5_long_MC %>%
  filter(fasta_id %in% repseqs_12SV5_filtered$fasta_id)

# Double check number of unique ASVs to make sure we've removed all the trimmed ones

length(unique(seqtab_12SV5_long_MC$fasta_id)) #Should be 6343
nrow(repseqs_12SV5_filtered) # Should be 6193

length(unique(seqtab_12SV5_long_MC$fasta_id)) # 6343
length(unique(seqtab_12SV5_long_MC_filtered$fasta_id)) # 6193

# List unique sample IDs to make sure we only have MC samples
# Use long file to assign control values

unique(seqtab_12SV5_long_MC_filtered$sample_id)

# seqtab_MC_12SV5_wide_filtered <- left_join(repseqs_12SV5_filtered, seqtab_12SV5, by = "fasta_id")
# seqtab_MC_12SV5_wide_filtered <- left_join(seqtab_MC_12SV5_wide_filtered, meta_lookup_MC, by = "fasta_id")

# The meta lookup can be used for both BLAST and Bayesian assignments
meta_lookup_MC <- seqtab_MC_12SV5_tax_long_meta_clean_filtered_blast %>%
  select(sample_id, Village, SubLocation, Site, Replicate, pH, Temp, EC, Control) %>%
  distinct()

#====Aggregate technical replicates====
# seqtab_MC_12SV5_tax_long_meta_clean_filtered <- seqtab_MC_12SV5_tax_long_meta_clean_filtered %>% 
#   group_by(SubLocation, Site, sscinames, Control, sample_id, fasta_id) %>% 
#   summarize(reads = mean(reads), .groups = "drop")

# Ok so this is the tricky part. 
# To summarize what I changed from the original MicroDecon run. 
  # First, I wanted to make sure we were not losing ASVs from taxonomic assignments, so I 
  # separated ASVs from taxa entirely. Ran MicroDecon fully before merging with BLAST assignments
  # or Sophie's bayes classifier. 
  # I did this becuase I wanted to look at ASV counts postdecon, but did not want taxonomy to interfere.
  # Realized that the way I did so before did not guarantee that assessment would be correct. 

#====MicroDecon====
# Run for the whole dataset, given that all locations share control replicates
# Now we convert our long table back to wide

seqtab_12SV5_wide_MC_filtered <- seqtab_12SV5_long_MC_filtered %>%
  select(-c(Control)) %>%
  pivot_wider(
    id_cols = fasta_id,
    names_from = sample_id,
    values_from = reads
  )

# Identify samples from LONG table
control_samples_MC <- seqtab_12SV5_long_MC_filtered %>%
  filter(Control == "Y") %>%
  distinct(sample_id) %>%
  pull(sample_id)

real_samples_MC <- seqtab_12SV5_long_MC_filtered %>%
  filter(Control == "N") %>%
  distinct(sample_id) %>%
  pull(sample_id)

##====Build decon table====
# Decon is run ONLY on fasta_id values. Species assignments have nothing to do with it.
seqtab_DECON_MC <- seqtab_12SV5_wide_MC_filtered %>%
  select(fasta_id, all_of(control_samples_MC), all_of(real_samples_MC)) %>%
  relocate(fasta_id, all_of(control_samples_MC)) %>%
  mutate(fasta_id = trimws(as.character(fasta_id))) %>%
  filter(!is.na(fasta_id), fasta_id != "") %>%
  group_by(fasta_id) %>%
  summarise(across(everything(), ~ sum(as.numeric(.), na.rm = TRUE)), .groups = "drop") %>%
  filter(rowSums(across(-fasta_id)) > 0)

numb.blanks_MC <- length(control_samples_MC)
numb.ind_MC <- c(length(real_samples_MC))

##====MicroDecon Diagnostics====
cat("=== MOWROW CREEK DIAGNOSTICS ===\n")
cat("numb.blanks:", numb.blanks_MC, "\n")
cat("numb.ind:", numb.ind_MC, "\n")
cat("ncol seqtab_DECON_MC:", ncol(seqtab_DECON_MC), "\n")
cat("expected ncol:", 1 + numb.blanks_MC + sum(numb.ind_MC), "\n")
cat("duplicated fasta_id rows:", sum(duplicated(seqtab_DECON_MC$fasta_id)), "\n")
cat("Dims:", dim(seqtab_DECON_MC), "\n")
cat("Any duplicated column names:", any(duplicated(colnames(seqtab_DECON_MC))), "\n")
cat("Any NA fasta_id:", any(is.na(seqtab_DECON_MC$fasta_id)), "\n")

stopifnot(
  ncol(seqtab_DECON_MC) == 1 + numb.blanks_MC + sum(numb.ind_MC),
  sum(duplicated(seqtab_DECON_MC$fasta_id)) == 0
)

original_colnames_MC <- colnames(seqtab_DECON_MC)
colnames(seqtab_DECON_MC) <- make.names(colnames(seqtab_DECON_MC), unique = TRUE)

##====Run MicroDecon====
decontaminated_MC <- decon(
  data        = as.data.frame(seqtab_DECON_MC),
  numb.blanks = numb.blanks_MC,
  numb.ind    = numb.ind_MC,
  taxa        = FALSE
)

##====Remove NULL values====
seqtab_12SV5_MC_wide_postdecon <- decontaminated_MC$decon.table
seqtab_12SV5_MC_wide_postdecon$Mean.blank <- NULL

if (!is.null(decontaminated_MC$decon.table) &&
    ncol(decontaminated_MC$decon.table) == length(original_colnames_MC)) {
  colnames(decontaminated_MC$decon.table) <- original_colnames_MC
}

##====Save MicroDecon Outputs====

# Make a table without taxa assignments
seqtab_12SV5_MC_wide_postdecon_filtered <- decontaminated_MC$decon.table %>%
  select(-Mean.blank)

# Make a df for both BLAST and Bayesian
seqtab_12SV5_MC_wide_postdecon_blast <- decontaminated_MC$decon.table %>%
  select(-Mean.blank) %>%
  left_join(tax_lookup_MC_BLAST, by = "fasta_id") %>%
  filter(rowSums(across(-c(fasta_id, sscinames))) > 0)

seqtab_12SV5_MC_wide_postdecon_bayes <- decontaminated_MC$decon.table %>%
  select(-Mean.blank) %>%
  left_join(tax_lookup_MC_bayes, by = "fasta_id") %>%
  filter(rowSums(across(-c(fasta_id, sscinames))) > 0)

reads_removed_MC <- decontaminated_MC$reads.removed
seqtab_12SV5_wide_MC_filtered_postdecon <- decontaminated_MC$decon.table
saveRDS(decontaminated_MC, "/Users/jacksadler/Desktop/project_thesis/taxonomy/12SV5_BLAST/Full_Outputs_MowrowCreekBLAST/Condensed_Code_June_Tweaks/decontaminated_MC.rds")

##====MICRODECON ADJUSTMENTS====

# Yes, the following two steps are a tad redundant.
###====Import/rename the decontaminated object====
#decontaminated_MC <- readRDS("/Users/jacksadler/Desktop/project_thesis/taxonomy/12SV5_BLAST/Full_Outputs_MowrowCreekBLAST/Condensed_Code_June_Tweaks/decontaminated_MC.rds")

###====Extract the decon table====
seqtab_12SV5_wide_MC_filtered_postdecon <- decontaminated_MC$decon.table

###====Move fasta_id out of the body and assign as a row name====
rownames(seqtab_12SV5_MC_wide_postdecon_filtered) <-
  seqtab_12SV5_MC_wide_postdecon_filtered$fasta_id
seqtab_12SV5_MC_wide_postdecon_filtered$fasta_id <- NULL

###====TRANSPOSE====
seqtab_12SV5_MC_long_postdecon_filtered <-
  t(seqtab_12SV5_MC_wide_postdecon_filtered)

seqtab_12SV5_MC_long_postdecon_filtered <-
  data.frame(seqtab_12SV5_MC_long_postdecon_filtered, check.names = FALSE)

###====CONVERT ALL COLUMNS TO NUMERIC====
i <- c(1:ncol(seqtab_12SV5_MC_long_postdecon_filtered))
seqtab_12SV5_MC_long_postdecon_filtered[, i] <- apply(
  seqtab_12SV5_MC_long_postdecon_filtered[, i],
  2,
  function(x) as.numeric(as.character(x))
)

###====REMOVE SAMPLES THAT NOW HAVE NO READS IN THEM====
seqtab_12SV5_MC_long_postdecon_filtered <-
  seqtab_12SV5_MC_long_postdecon_filtered[
    rowSums(seqtab_12SV5_MC_long_postdecon_filtered) > 0, ]

###====CONVERT POSTDECON WIDE -> LONG====
seqtab_12SV5_MC_long_postdecon_filtered <-
  tibble::rownames_to_column(seqtab_12SV5_MC_long_postdecon_filtered, "sample_id")

seqtab_12SV5_MC_long_postdecon_filtered <-
  gather(
    seqtab_12SV5_MC_long_postdecon_filtered,
    fasta_id,
    reads,
    -sample_id,
    factor_key = TRUE
  )

###====ENSURE READS ARE NUMERIC====
seqtab_12SV5_MC_long_postdecon_filtered$reads <-
  as.numeric(seqtab_12SV5_MC_long_postdecon_filtered$reads)

###====REMOVE ZERO-READ COMBINATIONS====
seqtab_12SV5_MC_long_postdecon_filtered <-
  seqtab_12SV5_MC_long_postdecon_filtered %>%
  filter(reads > 0)

###====Clean up Sample IDs====
# This is the file we will use for ASV diversity assessments, being sure to remove
# all commensal ASVs

seqtab_12SV5_MC_long_postdecon_filtered$sample_id <-
  gsub("\\.Set", "-Set", seqtab_12SV5_MC_long_postdecon_filtered$sample_id)

#====TAXONOMY ASSIGNMENTS====
# Here we will readd taxonomy to postdecon file and filter/adjust 
# classifier outputs

##----Bayes Classifier Output Cleanup----

asv_bayes_vsearch <- JackASVs_12SV5_vsearch_bayesian
vsearch100_tax <- vsearch100_taxonomy_Sadler_12SV5_final
bayesian_tax <- bayesian_taxonomy_Sadler_12SV5_final

tax_ranks <- c("K", "P", "C", "O", "F", "G", "S")

# Strip QIIME-style prefixes, whitespace, and coerce blanks/string "NA" to NA
asv_bayes_vsearch <- asv_bayes_vsearch %>%
  mutate(across(all_of(tax_ranks), ~ gsub("[a-z]__", "", .x))) %>%
  mutate(across(all_of(tax_ranks), ~ trimws(.x))) %>%
  mutate(across(all_of(tax_ranks), ~ na_if(.x, ""))) %>%
  mutate(across(all_of(tax_ranks), ~ na_if(.x, "NA")))
  
asv_bayes_vsearch <- asv_bayes_vsearch %>%
  select(-c(X))

# Standardize inconsistencies

asv_bayes_vsearch <- asv_bayes_vsearch %>%
  mutate(
    C = case_when(
      C == "Actinopterygii" ~ "Actinopteri",  # standardize to shorter form
      TRUE ~ C
    )
  )

# Drop species assignments with no genus

asv_bayes_vsearch <- asv_bayes_vsearch %>%
  mutate(S = if_else(is.na(G), NA_character_, S))

# Combine genus + species into a binomial name 
# e.g. G = "Staphylococcus", S = "aureus" → S = "Staphylococcus aureus"

asv_bayes_vsearch <- asv_bayes_vsearch %>%
  mutate(
    S = case_when(
      !is.na(G) & !is.na(S) ~ paste(G, S),
      TRUE ~ S
    )
  )

# Assign LCA label 
# Most resolved non-NA rank for each ASV
asv_bayes_vsearch <- asv_bayes_vsearch %>%
  mutate(lca_label = coalesce(S, G, F, O, C, P, K)) %>%
  rename(fasta_id = OTUid)

##====BLAST Output Cleanup====

blast_12SV5 <- `12SV5_blast_full_filtered`

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

filtered_12SV5 <- blast_12SV5 %>%
  filter(blast_12SV5$evalue <= 1e-3) %>%      # Keep only e-values <= 1e-3
  filter(blast_12SV5$pident >= 98) # Keep only matches with >= 95% identity

species_list <- unique(filtered_12SV5$sscinames)

species_counts_12SV5 <- filtered_12SV5 %>%
  count(sscinames, sort = TRUE) %>%
  rename(Species = sscinames, Count = n)

best_hits_12SV5 <- filtered_12SV5 %>%
  group_by(fasta_id) %>%
  slice_min(evalue, n = 5, with_ties = TRUE)

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

# Calculate a per-sample 0.01% read threshold: ASVs below this proportion will be treated as noise. Moved down for 12S from 0.05
seqtab_12SV5_tax_col_sum<-colSums(seqtab_12SV5_tax[c(13:219)])
seqtab_12SV5_tax_col_sum<-data.frame(seqtab_12SV5_tax_col_sum)
seqtab_12SV5_tax_col_sum$read_filter<-seqtab_12SV5_tax_col_sum$seqtab_12SV5_tax_col_sum *0.0001
seqtab_12SV5_tax_col_sum<- seqtab_12SV5_tax_col_sum %>% rownames_to_column("ID")

seqtab_12SV5_tax_long <- seqtab_12SV5_tax %>%
  pivot_longer(
    -c(fasta_id, pident, length, mismatch, gapopen, evalue, 
       bitscore, staxids, sscinames, scomnames, sskingdoms, genus),
    names_to = "Sample",
    values_to = "reads"
  ) %>%
  mutate(reads = as.numeric(reads))

metadata <- MC_metadata_12SV5_0329

metadata_summary <- metadata %>%
  select(sample_id, Month_Collected, Control, Village, Site, Replicate, SubLocation, 
         GPS_Location, Liters_Pumped, pH, 
         Temp, EC, Region)

metadata_summary<-data.frame(metadata_summary)
write.csv(metadata_summary,"metadata_summary.csv")

seqtab_12SV5_tax_long<-seqtab_12SV5_tax_long %>% rename(sample_id = Sample)
seqtab_12SV5_tax_col_sum<-seqtab_12SV5_tax_col_sum %>% rename(sample_id = ID)
seqtab_12SV5_tax_long$sample_id <- gsub("\\.", "-", seqtab_12SV5_tax_long$sample_id)
seqtab_12SV5_tax_col_sum$sample_id <- gsub("\\.", "-", seqtab_12SV5_tax_col_sum$sample_id)

seqtab_12SV5_tax_long_meta <- merge(metadata[, c("sample_id", "Month_Collected", 
                                                 "Village", "Site", "Replicate", "SubLocation", "GPS_Location", "Liters_Pumped", "pH",
                                                 "Temp", "EC", "Control")],seqtab_12SV5_tax_long, by="sample_id")
sum_reads<- sum(seqtab_12SV5_tax_long_meta$reads)

seqtab_MC_12SV5_tax_long_meta <- 
  seqtab_12SV5_tax_long_meta %>%
  filter(Village == "MowrowCreek")

seqtab_MC_12SV5_tax_long_meta_clean_filtered <- 
  merge(seqtab_12SV5_tax_col_sum[, c("sample_id", "read_filter")],
        seqtab_12SV5_tax_long_meta_MowrowCreek, by="sample_id")

seqtab_MC_12SV5_tax_long_meta_clean_filtered<-
  seqtab_12SV5_tax_long_meta_clean_filtered_MowrowCreek[!(seqtab_12SV5_tax_long_meta_clean_filtered_MowrowCreek$reads<=seqtab_12SV5_tax_long_meta_clean_filtered_MowrowCreek$read_filter),]

seqtab_MC_12SV5_tax_long_meta_clean_filtered<-
  seqtab_12SV5_tax_long_meta_clean_filtered_MowrowCreek[!(seqtab_12SV5_tax_long_meta_clean_filtered_MowrowCreek$reads<5),]

write.csv(seqtab_12SV5_tax_long_meta_clean_filtered_MowrowCreek,"MowrowCreek_ASV_table_long_filtered.csv")
seqtab_12SV5_tax_long_meta_clean_filtered_MowrowCreek <- read_csv("MowrowCreek_ASV_table_long_filtered.csv")

##====Taxonomizr ident buildout====

# if using EHD
#prepareDatabase('/Volumes/YourDriveName/accessionTaxa.sql')

blastResults<-read.table('XXXX.blast',header=FALSE,stringsAsFactors=FALSE)
#grab the 4th |-separated field from the reference name in the second column
accessions<-sapply(strsplit(blastResults[,2],'\\|'),'[',4)

taxaId<-accessionToTaxa(c("LN847353.1","AL079352.3"),"accessionTaxa.sql")
print(taxaId)

getTaxonomy(taxaId,'accessionTaxa.sql')

##====Tax lookup====
tax_lookup_MC_BLAST <- seqtab_MC_12SV5_tax_long_meta_clean_filtered %>%
  count(fasta_id, sscinames, sort = TRUE) %>%
  group_by(fasta_id) %>%
  slice(1) %>%
  ungroup() %>%
  select(fasta_id, sscinames)

tax_lookup_MC_Bayes <- asv_bayes_vsearch %>%
  rename(species = S)

# USE tax_lookup for fasta and taxid filtering with species lists

##====JOIN TAXONOMY BACK====
# Make sure our sample id's in postdecon match metadata. 
# Maybe make one BLAST version and one Bayesian 

# Use the following for ASV, read counts, and classifier assignments
# at the ASV level: seqtab_12SV5_MC_long_postdecon_filtered
# Then use the BLAST and BAYES versions for classification assessment

seqtab_12SV5_MC_long_postdecon_filtered$sample_id <-
  gsub("\\.Set", "-Set", seqtab_12SV5_MC_long_postdecon_filtered$sample_id)

#seqtab_12SV5_MC_long_postdecon_filtered_BLAST <-
#  left_join(
#    seqtab_12SV5_MC_long_postdecon_filtered,
#    tax_lookup_MC_BLAST,
#    by = "fasta_id"
#  )

seqtab_12SV5_MC_long_postdecon_filtered_BLAST <-
  left_join(
  seqtab_12SV5_MC_long_postdecon_filtered,
  tax_lookup_MC_BLAST,
  by = "fasta_id"
 )

seqtab_12SV5_MC_long_postdecon_filtered_Bayes <-
  left_join(
    seqtab_12SV5_MC_long_postdecon_filtered,
    tax_lookup_MC_Bayes,
    by = "fasta_id"
  )

##====JOIN SAMPLE METADATA BACK====

seqtab_12SV5_MC_long_postdecon_filtered_Bayes <-
  left_join(
    seqtab_12SV5_MC_long_postdecon_filtered_Bayes,
    meta_lookup_MC,
    by = "sample_id"
  )

# Count number of unique ASVs for each df.
# They should be the same

n_distinct(seqtab_12SV5_MC_long_postdecon_filtered_Bayes$fasta_id) 
# 2,590

seqtab_12SV5_MC_long_postdecon_filtered_BLAST <-
  left_join(
    seqtab_12SV5_MC_long_postdecon_filtered_BLAST,
    meta_lookup_MC,
    by = "sample_id"
  )

n_distinct(seqtab_12SV5_MC_long_postdecon_filtered_BLAST$fasta_id) 
# 2,590

##====NA order assignments====
# We're starting with NAs at the order level
###====Bayes====
# New data frame with NA assignments
bayes_NA_postdecon <- seqtab_12SV5_MC_long_postdecon_filtered_Bayes %>% 
  filter(is.na(O))
n_distinct(bayes_NA_postdecon$fasta_id)
# 183 at the order level
# 442 at the family level
# 646 at the genus level
# 1036 at the species level

# New data frame with species assignments
bayes_assignments_postdecon <- seqtab_12SV5_MC_long_postdecon_filtered_Bayes %>% 
  filter(!(is.na(O)))
n_distinct(bayes_assignments_postdecon$fasta_id)
# 2,407 at the order level
# 2,148 at the family level
# 1,944 at the genus level
# 1,554 at the species level

ASV_tax_na_losses_bayes <- data.frame(
  Rank = c("None", "Order", "Family", "Genus", "Species"),
  ASV_Count = c(2590, 2407, 2148, 1944, 1554)
)

###====BLAST====

# First we need to expand BLAST output with Taxonkit
#taxonkit_out <- readLines("taxonkit_out.csv")
#taxonkit_clean <- str_remove_all(taxonkit_out, '^"|"$')
#taxonkit_clean <- taxonkit_clean[-1]

#taxonkit_df <- read.delim(
#  text = paste(taxonkit_clean, collapse = "\n"),
#  sep = "\t",
#  header = FALSE,
#  quote = "",
#  stringsAsFactors = FALSE
#)
#colnames(taxonkit_df) <- c("fasta_id", "staxids", "lineage")

#taxonkit_df <- taxonkit_df %>%
#  filter(fasta_id != "fasta_id")
#taxonkit_df <- taxonkit_df %>%
#  separate(
#    lineage,
#    into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
#    sep = ";",
#    fill = "right",
#    extra = "drop"
#  )
#taxonkit_df <- taxonkit_df[-c(1), ]
#taxonkit_df[] <- lapply(taxonkit_df, function(col) gsub('^"|"$', '', col))

# Merge taxonkit with filtered BLAST output

# Now here we need to run LCA to condense BLAST outputs

taxonkit_df_lca <- taxonkit_df %>%
  group_by(fasta_id) %>%
  summarise(
    across(c(kingdom, phylum, class, order, family, genus, species),
           ~ if (n_distinct(.) == 1) first(.) else NA_character_),
    staxids = first(staxids),  
    .groups = "drop"
  )

seqtab_12SV5_MC_long_postdecon_filtered_BLAST <- left_join(
  seqtab_12SV5_MC_long_postdecon_filtered_BLAST,
  taxonkit_df_lca,
  by = "fasta_id",
)

seqtab_12SV5_MC_long_postdecon_filtered_BLAST <- left_join(
  seqtab_12SV5_MC_long_postdecon_filtered_BLAST,
  tax_lookup_MC_BLAST,
  by = "fasta_id"
)

#====WORK FROM HERE====
# Now we want to resolve any differences between taxonkit and BLAST. If at all?

mismatches <- seqtab_12SV5_MC_long_postdecon_filtered_BLAST %>%
  filter(!is.na(species) & !is.na(sscinames) & species != sscinames)

print(mismatches)
nrow(mismatches)

seqtab_12SV5_MC_long_postdecon_filtered_BLAST <- seqtab_12SV5_MC_long_postdecon_filtered_BLAST %>%
  mutate(species_sscinames = coalesce(species, sscinames))

tax_lookup_MC_BLAST <- rows_update(
  tax_lookup_MC_BLAST,
  tax_reference %>% select(sscinames, class, order, family),
  by = "sscinames",
  unmatched = "ignore"
)

# New data frame with NA assignments
BLAST_NA_postdecon <- seqtab_12SV5_MC_long_postdecon_filtered_BLAST %>% 
  filter(is.na(species_sscinames))
n_distinct(BLAST_NA_postdecon$fasta_id)
# 783 at the order level
# 1,745 at the family level
# 2,001 at the genus level
# 2,240 at the species level

# New data frame with assignments up to NA at Order level
BLAST_assignments_postdecon <- seqtab_12SV5_MC_long_postdecon_filtered_BLAST %>% 
  filter(!(is.na(order)))
n_distinct(BLAST_assignments_postdecon$fasta_id)
# 1,807 at the orspecies_sscinames# 1,807 at the order level
# 845 at the family level
# 589 at the genus level
# 350 at the species level

# Ok here is what I think is happening. Taxonkit is dumber 
# than I was expecting. It will not in fact fill out the taxonomy
# as we need it to. I will need to manually merge the family trees
# as I have them recorded from the species assignments. 
# sscinames is the most accurate and correct species detection level
# from BLAST. species is just what taxonkit has access to and 
# they are not the same thing
# I need to fill in tax_lookup_MC_BLAST in order to run my 
# family and order level assignments

ASV_tax_na_losses_BLAST <- data.frame(
  Rank = c("None", "Order", "Family", "Genus", "Species"),
  ASV_Count = c(2590, 1807, 845, 589, 350)
)

#combined_ASV_tax_na_losses <- bind_rows(
#  "Bayes" = ASV_tax_na_losses_bayes, 
#  "BLAST" = ASV_tax_na_losses_BLAST, 
#  .id = "Source"
#)

combined_ASV_tax_na_losses$Rank <- factor(
  combined_ASV_tax_na_losses$Rank,
  levels = c("None", "Order", "Family", "Genus", "Species") 
)

combined_ASV_tax_na_losses %>%
  mutate(Rank = factor(Rank, levels = c("None", "Order", "Family", "Genus", "Species"))) %>%
  ggplot(aes(x = Rank, y = ASV_Count, fill = Source)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Taxonomic Rank", y = "ASV Count",
       title = "Unique ASV Classification Across Taxonomic Ranks",
       fill = "Classifier") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        legend.position = "right",
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(face = "bold"))

##====Commensal species removal====

remove_species_BLAST <- function(df, species_vec, genus_vec) {
  removed_df <- df %>%
    filter(species_sscinames %in% species_vec | genus %in% genus_vec)
  
  filtered_df <- df %>%
    filter(!species_sscinames %in% species_vec & !genus %in% genus_vec)
  
  list(
    filtered = filtered_df,
    removed = removed_df
  )
}

remove_species_bayes <- function(df, species_vec, genus_vec) {
  
  removed_df <- df %>%
    filter(S %in% species_vec | G %in% genus_vec)
  
  filtered_df <- df %>%
    filter(!S %in% species_vec & !G %in% genus_vec)
  
  list(
    filtered = filtered_df,
    removed = removed_df
  )
}

species_to_remove <- c("Homo sapiens", "Pan troglodytes", "Canis lupus", 
                       "Gallus gallus", "Felis cattus", "Bos taurus", 
                       "Mus musculus", "Felis catus", "Sus scrofa", "Rattus rattus")

genera_to_remove <- c()

commensal_filter_BLAST <- remove_species_BLAST(BLAST_assignments_postdecon, 
                                               species_to_remove, genera_to_remove)
commensal_filter_bayes <- remove_species_bayes(bayes_assignments_postdecon, 
                                               species_to_remove, genera_to_remove)

BLAST_assignments_postdecon_filtered <- commensal_filter_BLAST$filtered
n_distinct(BLAST_assignments_species_postdecon_filtered$fasta_id)
# 2,547 ASVs remaining
removed_BLAST_12SV5 <- commensal_filter_BLAST$removed
n_distinct(removed_BLAST_12SV5$fasta_id)
# 242 ASVs removed

bayes_assignments_postdecon_filtered <- commensal_filter_bayes$filtered
n_distinct(bayes_assignments_postdecon_filtered$fasta_id)
# 2,068 ASVs remaining
removed_bayes_12SV5 <- commensal_filter_bayes$removed
n_distinct(removed_bayes_12SV5$fasta_id)
# 339 ASVs removed

# For generating diversity plots, we want to make sure NA asvs are included
# Will need to filter out ASVs from commensal

###====Generate a list of ASV fasta_ids for commensal removals====
# Merge removed_BLAST_12SV5 and removed_bayes_12SV5
# Generate a list of distinct ASVs and then filter by fasta_id, keep read counts
# Plot diversity metrics with all ASVs, even if NA past order level.

removed_commensal_ASVs_BLASTandbayes <- bind_rows(removed_bayes_12SV5,
                                                  removed_BLAST_12SV5) %>%
  distinct(fasta_id, .keep_all = TRUE)

commensal_fasta_id <- c(removed_commensal_ASVs_BLASTandbayes$fasta_id)

# 364 commensal ASVs were identified and removed from our ASV list for diversity analysis

diversity_matrix_MC_ASV <- seqtab_12SV5_MC_long_postdecon_filtered %>%
  filter(!fasta_id %in% commensal_fasta_id)

# This is where we want it all to match, and if not need to look into redoing 
# beta div assessment

n_distinct(diversity_matrix_MC_ASV$fasta_id) # 2226 unique fasta_ids 

#====Species Lists for each Classifier====

n_distinct(BLAST_assignments_postdecon_filtered$species_sscinames) # 53 (remove NA)
species_list_BLAST <- 
  as.list(unique(BLAST_assignments_postdecon_filtered$species_sscinames))

n_distinct(bayes_assignments_postdecon_filtered$S) # 110 (remove NA)
species_list_bayes <- 
  as.list(unique(bayes_assignments_postdecon_filtered$S))

n_distinct(BLAST_assignments_postdecon_filtered$genus) # 27
n_distinct(bayes_assignments_postdecon_filtered$G) # 110

species_bayes <- as.data.frame(do.call(rbind, species_list_bayes)) %>%
  rename(species = V1) %>% 
  filter(!is.na(species)) %>%
  left_join(tax_lookup_MC_Bayes, by = "species") %>%
  select(-fasta_id) %>%
  group_by(species) %>%
  summarise(
    across(everything(), ~ if (n_distinct(na.omit(.x)) > 1) {
      paste(na.omit(unique(.x)), collapse = " | ")
    } else {
      first(na.omit(.x))
    }),
    .groups = "drop"
  ) %>%
  rename(kingdom = K, phylum = P, class = C, order = O, family = F, genus = G)

species_BLAST <- as.data.frame(do.call(rbind, species_list_BLAST)) %>%
  rename(species = V1) %>% 
  filter(!is.na(species)) %>%
  left_join(taxonkit_df_lca, by = "species") %>%
  select(-fasta_id) %>%
  mutate(across(everything(), as.character)) %>%
  group_by(species) %>%
  summarise(
    across(everything(), ~ {
      vals <- na.omit(unique(.x))
      if (length(vals) == 0) NA_character_
      else if (length(vals) > 1) paste(vals, collapse = " | ")
      else vals
    }),
    .groups = "drop"
  )

# Note that taxonkit was not very useful here, only illuminating three out of our
# 21 species

BLAST_bayes_overlap <- as.data.frame(intersect(species_BLAST$species, species_bayes$species)) %>% 
  rename(species = "intersect(species_BLAST$species, species_bayes$species)")

shared_cols_tax <- c("species", "kingdom", "family", "class", "order", "genus")
BLAST_sub <- species_BLAST %>% select(all_of(shared_cols_tax))
bayes_sub <- species_bayes %>% select(all_of(shared_cols_tax))

BLAST_bayes_merged <- bind_rows(BLAST_sub, bayes_sub) %>%
  group_by(species) %>%
  summarise(
    across(everything(), ~ {
      vals <- na.omit(unique(.x))
      if (length(vals) == 0) NA_character_
      else if (length(vals) > 1) paste(vals, collapse = " | ")
      else vals
    }),
    .groups = "drop"
  )

#====Color Palettes====
##====Species color palette====

species_levels <- sort(unique(BLAST_bayes_merged$species))
species_colors <- setNames(
  colorRampPalette(brewer.pal(12, "Paired"))(length(species_levels)),
  species_levels
)

# It would be good too to have a number and specific ASV per species assignment parameter
# I think that can wait until later though
# By doing so, we could match fasta_id colors to species
# This is an "extra" for the thesis. Not critical but cool if you have the time.

##====Family Color Palette====

family_levels <- sort(unique(BLAST_bayes_merged$family))
family_colors <- setNames(
  colorRampPalette(brewer.pal(12, "Paired"))(length(family_levels)),
  family_levels
)

##====Order Color Palette====

order_levels <- sort(unique(BLAST_bayes_merged$order))
order_colors <- setNames(
  colorRampPalette(brewer.pal(12, "Paired"))(length(order_levels)),
  order_levels
)


##====fasta color palette====

fasta_levels <- sort(unique(diversity_matrix_MC_ASV$fasta_id))
fasta_colors <- setNames(
  colorRampPalette(brewer.pal(11, "Spectral"))(length(fasta_levels)),
  fasta_levels
)

fasta_species_map <- tibble(fasta_id = fasta_levels) %>%
  left_join(all_tax_lookup_resolved %>% select(fasta_id, species), by = "fasta_id") %>%
  left_join(all_tax_lookup_resolved %>% select(fasta_id, G), by = "fasta_id") %>%
  left_join(all_tax_lookup_resolved %>% select(fasta_id, F), by = "fasta_id") %>%
  left_join(all_tax_lookup_resolved %>% select(fasta_id, O), by = "fasta_id")

fasta_species_map <- fasta_species_map %>%
  mutate(
    category = case_when(
      !(fasta_id %in% all_tax_lookup_resolved$fasta_id) ~ "no_taxonomy_below_order",   
      is.na(G) ~ "family_only",
      is.na(F) ~ "order_only",
      is.na(species)                                     ~ "genus_only",   
      TRUE                                                ~ "species_known"
    )
  )

species_levels <- sort(unique(na.omit(fasta_species_map$species)))

species_colors <- setNames(
  colorRampPalette(brewer.pal(11, "Spectral"))(length(species_levels)),
  species_levels
)

species_colors_full <- c(species_colors,
                         "family_only" = "#854D65",
                         "order_only" = "#7B3B8B",
                         "genus_only"  = "#DAB1DA",
                         "no_taxonomy_below_order" = "black")

genera_levels <- sort(unique(na.omit(fasta_species_map$G)))
family_levels <- sort(unique(na.omit(fasta_species_map$F)))
order_levels <- sort(unique(na.omit(fasta_species_map$O)))

order_colors <- setNames(
  colorRampPalette(brewer.pal(9, "PuRd"))(length(order_levels)),
  order_levels
)

fasta_species_map <- fasta_species_map %>%
  mutate(color_key = ifelse(category == "species_known", species, category))

fasta_colors_by_species <- setNames(
  species_colors_full[fasta_species_map$color_key],
  fasta_species_map$fasta_id
)

# Should equal length(fasta_levels) - confirms every fasta_id got a color
length(fasta_colors_by_species)

# Should be 0 - confirms nothing fell through uncolored
sum(is.na(fasta_colors_by_species))

# Quick visual gut-check of category counts
table(fasta_species_map$category)

##====Site color palette====

site_levels <- sort(unique(meta_lookup_MC$Site))
site_colors <- setNames(
  colorRampPalette(brewer.pal(12, "Paired"))(length(site_levels)),
  site_levels
)

#====ASV Abundance Analysis====
# Commensal contaminants were not considered here. 
# We are looking at purely whether the PCR amplified our primers
# equally across replicates

asv_counts <- seqtab_12SV5_MC_long_postdecon_filtered %>%
  group_by(sample_id) %>%
  summarise(n_ASVs = n_distinct(fasta_id))

asv_counts_rep_merged <- asv_counts %>%
  mutate(sample_id = sub("\\.Set", "-Set", sample_id)) %>%
  mutate(sample_id = sub("_J-Set[0-9]+$", "", sample_id)) %>%
  mutate(site = str_extract(sample_id, "^K[0-9]+(?=_)")) %>%
  group_by(sample_id)

asv_counts_sample_merged <- asv_counts_rep_merged %>%
  mutate(sample_id = sub("_R[0-9]]+$", "", sample_id)) %>%
  mutate(site = str_extract(sample_id, "^K[0-9]+(?=_)")) %>%
  group_by(sample_id)

asv_counts_sample_merged <- aggregate(cbind(n_ASVs) ~ site, 
                                   data = asv_counts_rep_merged, 
                                   FUN = sum)
asv_counts_rep_merged <- asv_counts_rep_merged %>%
  mutate(site = factor(site,
                       levels = c("K1", "K2", "K4", "K5",
                                  "K6", "K8", "K9", "K10", "K11")))

ggplot(asv_counts_rep_merged,
       aes(x = sample_id, y = n_ASVs)) +
  geom_bar(stat = "identity") +
  facet_grid(~site, scales = "free_x", space = "free") +
  labs(x = "Sample Replicate", 
       y = "ASV Number",
       fill = "FASTA ID") +
  theme_bw(base_family = "Baskerville") +
  theme(axis.title = element_text(size = 18, face = "plain"),
        axis.text = element_text(size = 16, face = "plain"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        strip.text = element_text(size = 18, face = "plain"),
        legend.position = "right",
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.4, "cm")
  )

ggplot(asv_counts_sample_merged,
       aes(x = site, y = n_ASVs)) +
  geom_bar(stat = "identity") +
  facet_grid(~site, scales = "free_x", space = "free") +
  labs(x = "Sample", y = "ASV Number",
       title = "ASV Abundance - Per Sample",
       fill = "FASTA ID") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(face = "bold")) 

##====Plot ASV Abundances with FASTA====

y_max <- max(
  c(diversity_matrix_MC_ASV_Upper$reads,
    diversity_matrix_MC_ASV_Middle$reads,
    diversity_matrix_MC_ASV_Lower$reads),
  na.rm = TRUE
)

ggplot(diversity_matrix_MC_ASV_Upper,
       aes(x = sample_id, y = reads, fill = fasta_id)) +
  geom_bar(stat = "identity") +
  facet_grid(~SubLocation, scales = "free_x", space = "free") +
  scale_fill_manual(values = fasta_colors_by_species) +
  labs(x = "Sample Replicate", y = "Reads",
       fill = "FASTA ID") +
  theme_bw(base_size = 16, base_family = "Baskerville") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(hjust = 1, size = 14),
        strip.text = element_text(size = 18, face = "plain"),
        axis.title = element_text(size = 20),
        legend.position = "none",
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(face = "bold")) + 
  ylim(0, y_max)+ 
  coord_flip()

ggplot(diversity_matrix_MC_ASV_Middle,
       aes(x = sample_id, y = reads, fill = fasta_id)) +
  geom_bar(stat = "identity") +
  facet_grid(~SubLocation, scales = "free_x", space = "free") +
  scale_fill_manual(values = fasta_colors_by_species) +
  labs(x = "Sample Replicate", y = "Reads",
       fill = "FASTA ID") +
  theme_bw(base_size = 16, base_family = "Baskerville") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(hjust = 1, size = 14),
        strip.text = element_text(size = 18, face = "plain"),
        axis.title = element_text(size = 20),
        legend.position = "none",
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(face = "bold")) + 
  ylim(0, y_max)+ 
  coord_flip()

ggplot(diversity_matrix_MC_ASV_Lower,
       aes(x = sample_id, y = reads, fill = fasta_id)) +
  geom_bar(stat = "identity") +
  facet_grid(~SubLocation, scales = "free_x", space = "free") +
  scale_fill_manual(values = fasta_colors_by_species) +
  labs(x = "Sample Replicate", y = "Reads",
       fill = "FASTA ID") +
  theme_bw(base_size = 16, base_family = "Baskerville") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(hjust = 1, size = 14),
        strip.text = element_text(size = 18, face = "plain"),
        axis.title = element_text(size = 20),
        legend.position = "none",
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(face = "bold")) + 
  ylim(0, y_max) + 
  coord_flip()

##====Plot ASV Abundances with species====

# I need a lookup table to merge with each fasta_id, its corresponding species 
# assignment for both Bayes and BLAST

bayes_spec_asv_pairs <- unique(bayes_assignments_postdecon_filtered[, c("fasta_id", "S")]) %>%
  filter(!is.na(S) %>%
  rename(species = S)

BLAST_spec_asv_pairs <- unique(BLAST_assignments_species_postdecon_filtered[, c("fasta_id", "sscinames")]) %>%
  filter(!is.na(sscinames)) %>%
  rename(species = sscinames))

all_spec_asv_pairs <- full_join(BLAST_spec_asv_pairs, bayes_spec_asv_pairs)

# Check for duplicate assignments. Not that it matters really, but to be tidy
sum(duplicated(all_spec_asv_pairs))

# Repeat for family and order levels

bayes_family_asv_pairs <- unique(bayes_assignments_postdecon_filtered[, c("fasta_id", "F")]) %>%
  filter(!is.na(F)) %>%
  rename(order = F)
bayes_order_asv_pairs <- unique(bayes_assignments_postdecon_filtered[, c("fasta_id", "O")]) %>%
  filter(!is.na(O)) %>%
  rename(order = O)

BLAST_family_asv_pairs <- unique(BLAST_assignments_species_postdecon_filtered[, c("fasta_id", "family")]) %>%
  filter(!is.na(family)) %>%
           rename(species = sscinames))

all_spec_asv_pairs <- full_join(BLAST_spec_asv_pairs, bayes_spec_asv_pairs)



#====Divide div matrix into regions of study====

##====Merge metadata and species assignments back====
diversity_matrix_MC_ASV <- 
  left_join(diversity_matrix_MC_ASV, meta_lookup_MC, by = "sample_id") %>%
  select(-c("Control"))

diversity_matrix_MC_ASV_taxa <- left_join(diversity_matrix_MC_ASV, 
                                          all_spec_asv_pairs,
                                          relationship = "many-to-many")

diversity_matrix_MC_ASV_Upper <- diversity_matrix_MC_ASV %>%
  filter(SubLocation == "Upper")

diversity_matrix_MC_ASV_Middle <- diversity_matrix_MC_ASV %>%
  filter(SubLocation == "Middle")

diversity_matrix_MC_ASV_Lower <- diversity_matrix_MC_ASV %>%
  filter(SubLocation == "Lower")


# Species Assignments

diversity_matrix_MC_taxa_Upper <- diversity_matrix_MC_ASV_taxa %>%
  filter(SubLocation == "Upper") 

diversity_matrix_MC_taxa_Middle <- diversity_matrix_MC_ASV_taxa %>%
  filter(SubLocation == "Middle")

diversity_matrix_MC_taxa_Lower <- diversity_matrix_MC_ASV_taxa %>%
  filter(SubLocation == "Lower") 

ggplot(diversity_matrix_MC_taxa_Upper,
       aes(x = sample_id, y = reads, fill = species)) +
  geom_bar(stat = "identity") +
  facet_grid(~SubLocation, scales = "free_x", space = "free") +
  scale_fill_manual(values = species_colors) +
  labs(x = "Sample", y = "Reads",
       title = "ASV Abundance - Upper Mowrow Creek 12SV5",
       fill = "Species") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(face = "bold")) + 
  ylim(0, 300000)

ggplot(diversity_matrix_MC_taxa_Middle,
       aes(x = sample_id, y = reads, fill = species)) +
  geom_bar(stat = "identity") +
  facet_grid(~SubLocation, scales = "free_x", space = "free") +
  scale_fill_manual(values = species_colors) +
  labs(x = "Sample", y = "Reads",
       title = "ASV Abundance - Middle Mowrow Creek 12SV5",
       fill = "Species") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(face = "bold")) + 
  ylim(0, 300000)

ggplot(diversity_matrix_MC_taxa_Lower,
       aes(x = sample_id, y = reads, fill = species)) +
  geom_bar(stat = "identity") +
  facet_grid(~SubLocation, scales = "free_x", space = "free") +
  scale_fill_manual(values = species_colors) +
  labs(x = "Sample", y = "Reads",
       title = "ASV Abundance - Lower Mowrow Creek 12SV5",
       fill = "Species") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(face = "bold")) + 
  ylim(0, 300000)

# I think we also want to show family and order levels
# These may be better show with circle-plots

all_tax_lookup <- all_spec_asv_pairs %>%
  left_join(
    bayes_assignments_postdecon_filtered %>%
      select(fasta_id, O, F, G, lca_label),
    by = "fasta_id") %>%
  distinct()

all_tax_lookup %>% count(fasta_id) %>% filter(n > 1)

all_tax_lookup_resolved <- all_tax_lookup %>%
  group_by(fasta_id, O, F, G) %>%   
  summarise(
    species = if (n_distinct(species) > 1) NA_character_ else first(species),
    lca_label = first(G), 
    .groups = "drop"
  )

#====NORMALISING READS====

# We will normalize by sample replicate, but confirm whether that is correct to do or not.
# Following is replicate-specific normalization. PCR reps are merged. Discuss!!!

##====CALCULATE TOTAL READS PER BIOLOGICAL REPLICATE====

diversity_matrix_MC_ASV_taxa <- diversity_matrix_MC_ASV_taxa %>%
  mutate(replicate_id = paste(Site, Replicate, sep = "_"))

replicate_totals_MC_ASV_taxa <- 
  diversity_matrix_MC_ASV_taxa %>%
  group_by(replicate_id) %>%
  summarise(replicate_total_reads = sum(reads, na.rm = TRUE), .groups = "drop")

## NORMALIZE LONG TABLE BY replicate_id TOTAL and re-add replicate id column

diversity_matrix_MC_ASV_taxa_normalized <- 
  diversity_matrix_MC_ASV_taxa %>%
  left_join(replicate_totals_MC_ASV_taxa, by = "replicate_id") %>%
  left_join(all_tax_lookup_resolved %>%
              select(fasta_id, O, F, G),
            by = "fasta_id") %>%
  mutate(normalised_reads = reads / replicate_total_reads)

replicate_check <- 
  diversity_matrix_MC_ASV_taxa_normalized %>%
  group_by(replicate_id) %>%
  summarise(total_norm = sum(normalised_reads, na.rm = TRUE), .groups = "drop")

print(replicate_check, n=27)

genus_check <- all_tax_lookup %>%
  group_by(fasta_id) %>%
  summarise(n_rows = n(), n_genus = n_distinct(G), .groups = "drop")

resolvable_ids <- genus_check %>% filter(n_rows > 1, n_genus == 1) %>% pull(fasta_id)

unresolved_ids <- genus_check %>% filter(n_rows > 1, n_genus > 1) %>% pull(fasta_id)

cat("Resolvable by dropping to genus:", length(resolvable_ids), "\n")
cat("Still conflicting even at genus level:", length(unresolved_ids), "\n")

##====species Normalized====
### Collapse normalized reads within technical (PCR) replicates into one value per species

diversity_matrix_MC_ASV_taxa_normalized_species <- 
  diversity_matrix_MC_ASV_taxa_normalized %>%
  group_by(replicate_id, species, Site, SubLocation) %>%
  summarise(normalised_reads = sum(normalised_reads, na.rm = TRUE), .groups = "drop")

##====family Normalized====
diversity_matrix_MC_ASV_taxa_normalized_family <- 
  diversity_matrix_MC_ASV_taxa_normalized %>%
  group_by(replicate_id, F, Site, SubLocation) %>%
  summarise(normalised_reads = sum(normalised_reads, na.rm = TRUE), .groups = "drop")

##====order Normalized====
diversity_matrix_MC_ASV_taxa_normalized_order <- 
  diversity_matrix_MC_ASV_taxa_normalized %>%
  group_by(replicate_id, O, Site, SubLocation) %>%
  summarise(normalised_reads = sum(normalised_reads, na.rm = TRUE), .groups = "drop")

##====fasta Normalized====
### Collapse normalized reads within technical (PCR) replicates into one value per fasta_id

diversity_matrix_MC_ASV_taxa_normalized_fasta <- 
  diversity_matrix_MC_ASV_taxa_normalized %>%
  group_by(replicate_id, fasta_id, Site, SubLocation) %>%
  summarise(normalised_reads = sum(normalised_reads, na.rm = TRUE), .groups = "drop")

# Whole point of normalizing is to have an level mode of comparison across all locations

##====Plot normalized reads====

# Plot species and fasta_id legend

legend_species <- data.frame(
  species = names(species_colors_full),
  fill_color = species_colors_full
)

n_cols <- 3
col_width <- 10

legend_species <- legend_species %>%
  mutate(
    index = row_number() - 1,
    col = index %% n_cols,
    row = index %/% n_cols,
    x_pos = col * col_width  
  )

ggplot(legend_species, aes(y = -row)) +
  geom_tile(aes(x = x_pos, fill = species), width = 0.8, height = 0.6, colour = "white") +
  geom_text(aes(x = x_pos + 0.55, label = species), hjust = 0, size = 3.2, fontface = "italic") +
  scale_fill_manual(values = species_colors_full) +
  xlim(-0.5, n_cols * col_width + 1) +
  theme_void() +
  theme(legend.position = "none")

ggplot(diversity_matrix_MC_ASV_taxa_normalized_fasta, 
    aes(x = replicate_id, fill = fasta_id, y = normalised_reads)) + 
  scale_fill_manual(values = fasta_colors_by_species) +
  geom_bar(stat = "identity", colour = "white")+
  labs(
    x = "Replicate ID",
    y = "Normalized ASV Abundances"
  ) +
  theme_bw(base_size = 16, base_family = "Baskerville") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, 
                                   vjust = 1, size = 14),
        legend.position = "none",
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 18, face = "plain"),
        axis.title = element_text(size = 20)) + 
  facet_wrap(~ SubLocation, scales = "free_x") 

ggplot(diversity_matrix_MC_ASV_taxa_normalized_species, 
       aes(x = replicate_id, fill = species, y = normalised_reads)) + 
  scale_fill_manual(values = species_colors_full) +
  geom_bar(stat = "identity", colour = "white")+
  labs(
    x = "Replicate ID",
    y = "Normalized ASV Abundances"
  ) +
  theme_bw(base_size = 16, base_family = "Baskerville") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, 
                                   vjust = 1, size = 14),
        legend.position = "none",
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 18, face = "plain"),
        axis.title = element_text(size = 20)) + 
  facet_wrap(~ SubLocation, scales = "free_x") 

ggplot(diversity_matrix_MC_ASV_taxa_normalized_family, 
       aes(x = replicate_id, fill = F, y = normalised_reads)) + 
  scale_fill_manual(values = family_colors) +
  geom_bar(stat = "identity", colour = "white")+
  labs(
    x = "Replicate ID",
    y = "Normalized ASV Abundances"
  ) +
  theme_bw(base_size = 16, base_family = "Baskerville") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, 
                                   vjust = 1, size = 14),
        legend.position = "bottom",
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 18, face = "plain"),
        axis.title = element_text(size = 20)) + 
  facet_wrap(~ SubLocation, scales = "free_x") 

ggplot(diversity_matrix_MC_ASV_taxa_normalized_order, 
       aes(x = replicate_id, fill = O, y = normalised_reads)) + 
  scale_fill_manual(values = order_colors) +
  geom_bar(stat = "identity", colour = "white")+
  labs(
    x = "Replicate ID",
    y = "Normalized ASV Abundances"
  ) +
  theme_bw(base_size = 16, base_family = "Baskerville") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, 
                                   vjust = 1, size = 14),
        legend.position = "bottom",
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 18, face = "plain"),
        axis.title = element_text(size = 20)) + 
  facet_wrap(~ SubLocation, scales = "free_x") 

# write.csv(seqtab_MC_12SV5_tax_normalised_replicate,"NORM_MC_ASV_table_LONG_filtered.csv")

# Now let's make a version that plots according to genus, family, and order
# Perhaps another one of these would be good for taxonomy of fishes, etc

library(treemapify)

order_abundance <- diversity_matrix_MC_ASV_taxa_normalized %>%
  group_by(O) %>%
  mutate(O = ifelse(is.na(O), "Unassigned", O)) %>%
  summarise(normalised_reads = sum(normalised_reads), .groups = "drop")

ggplot(order_abundance, 
       aes(area = normalised_reads, fill = O, label = O)) +
  scale_fill_manual(values = order_colors, na.value = "black") +
  labs(title = "Relative Abundance of Classified Orders Across Mowrow Creek") +
  geom_treemap() +
  geom_treemap_text(colour = "white",
                    place = "centre",
                    reflow = T,
                    fontface = "italic")



##====CONVERTING NORMALISED LONG FORMAT INTO WIDE FORMAT====
# Rows are replicate_id
# Columns are fasta_id or species
# Values are normalized abundance

diversity_matrix_MC_ASV_taxa_normalized_wide <- diversity_matrix_MC_ASV_taxa_normalized %>%
  select(-c(Village, sample_id, SubLocation, Site, SubLocation, pH, EC, Temp, reads, Replicate, replicate_total_reads))

diversity_matrix_MC_ASV_taxa_normalized_wide_fasta_id <- diversity_matrix_MC_ASV_taxa_normalized_wide %>%
  select(-c(species))

diversity_matrix_MC_ASV_taxa_normalized_wide_fasta_id <- 
  diversity_matrix_MC_ASV_taxa_normalized %>%
  pivot_wider(
    names_from = fasta_id,
    values_from = normalised_reads,
    values_fill = 0
  )

replicate_metadata_MC <- 
  MC_metadata_12SV5_0329 %>%
  group_by(replicate_id) %>%
  summarise(
    sample_id = first(sample_id),
    Village = first(Village),
    Site = first(Site),
    Month_Collected = first(Month_Collected),
    pH = first(pH),
    Temp = first(Temp),
    EC = first(EC),
    Region = first(Region),
    .groups = "drop"
  )

# Use wide formats for NMDS/Permanova
# Use long formats for plotting stacked bar plots

# You only need merged data if wide and long dont already have metadata attached

# merged_data_MC_wide <- 
#  left_join(replicate_metadata_MC, seqtab_MC_12SV5_tax_normalised_replicate_wide, by = "replicate_id")

# merged_data_MC_long <- 
#  left_join(replicate_metadata_MC, seqtab_MC_12SV5_tax_normalised_replicate, by = "replicate_id")

# ggplot(seqtab_MC_12SV5_tax_normalised_replicate_wide %>%
#         filter(!grepl("^K_C", replicate_id)),
#       aes(x = replicate_id, fill = species, y = normalised_reads)) +
#  geom_bar(stat = "identity", colour = "white") +
#  labs(title = "Normalized ASV Abundance - Mowrow Creek 12SV5",
#       x = "Site",
#       y = "Relative abundance") +
#  theme(axis.text.x = element_text(angle = 90),
#        legend.position = "none") +
#  facet_wrap(~Site, scales = "free_x")

# Need to check the format of my wide tables from original code to make sure
# we're in the same ballpark.  

#====PERMANOVA (ALL TAXA, THEN PER-PHYLUM SUBSETS)====

# seqtab_MC_12SV5_tax_long_meta_clean_postdecon <- 
#   seqtab_MC_12SV5_tax_long_meta_clean_postdecon %>%
#   mutate(
#     sample_id = gsub("\\.Set", "-Set", sample_id),
#     replicate_id = sub("-Set[0-9]+$", "", sample_id)
#   )

##====1. Row Bind Long, non-normalized but decontaminated datasets====

all_long_MC <- diversity_matrix_MC_ASV_taxa %>%
  mutate(replicate_id = sub("-Set[0-9]+$", "", sample_id))
# adjust name as needed

##====2. Aggregate reads per replicate_id====
#(in case multiple ASVs map to same species)
# Running diversity on an ASV basis, change as needed for on a species level

all_long_agg_MC_fasta <- all_long_MC %>%
  group_by(replicate_id, fasta_id) %>%
  summarise(reads = sum(reads), .groups = "drop")

###====2b. Buid one metadata row per replicate, replacing direct sample id matching====

metadata_rep_MC <- MC_metadata_12SV5_0329 %>%
  mutate(replicate_id = sub("-Set[0-9]+$", "", sample_id)) %>%
  group_by(replicate_id) %>%
  summarise(
    sample_id = first(sample_id),
    Village = first(Village),
    Site = first(Site),
    Replicate = first(Replicate),
    SubLocation = first(SubLocation),
    Month_Collected = first(Month_Collected),
    pH = first(pH),
    Temp = first(Temp),
    EC = first(EC),
    Region = first(Region),
    .groups = "drop"
  )

merged_data_MC <- merge(
  metadata_rep_MC[, c("replicate_id", "sample_id", "Village", "Site", "Replicate",
                      "Month_Collected", "pH", "Temp", "EC", "Region")],
  all_long_agg_MC_fasta,
  by = "replicate_id"
)

##====3. Pivot to wide format====
# (samples x species)
# all_wide_MC_species <- all_long_agg_MC %>%
#  pivot_wider(names_from = species, values_from = reads, values_fill = 0)

all_wide_MC_fasta <- all_long_agg_MC_fasta %>%
  pivot_wider(names_from = fasta_id, values_from = reads, values_fill = 0)

##====4. Merge with metadata====
merged_data_MC_fasta <- merge(
  metadata_rep_MC[, c("replicate_id", "sample_id", "Site", "SubLocation", "Month_Collected",
                      "pH", "Temp", "EC", "Region")],
  all_wide_MC_fasta,
  by = "replicate_id"
)

merged_data_MC_fasta <- na.omit(merged_data_MC_fasta)

##====5. Split into metadata and community matrix====
n_meta_cols <- 9  # adjust to match number of metadata columns selected above

adonis_meta_MC_fasta <- merged_data_MC_fasta[, 1:n_meta_cols]
adonis_comm_MC_fasta <- merged_data_MC_fasta[, (n_meta_cols + 1):ncol(merged_data_MC_fasta)]

##====6. Drop any all-zero fasta/species columns====
adonis_comm_MC_fasta <- adonis_comm_MC_fasta[, colSums(adonis_comm_MC_fasta) > 0]

##====7. Run PERMANOVA====
adon.results_MC_fasta <- adonis2(
  adonis_comm_MC_fasta ~ SubLocation,
  data = adonis_meta_MC_fasta,
  method = "raup",
  perm = 999,
  by = "margin"
)

print(adon.results_MC_fasta)

disp_MC_fasta <- betadisper(vegdist(adonis_comm_MC_fasta, method = "raup"),
                            adonis_meta_MC_fasta$SubLocation)
anova(disp_MC_fasta)

plot(disp_MC_fasta)

#----PLOTTING BETA DIVERSITY NMDS (Use Normalized read counts)----

replicate_totals_MC <- all_long_MC %>%
  group_by(replicate_id) %>%
  summarise(replicate_total_reads = sum(reads, na.rm = TRUE), .groups = "drop")

all_long_norm_MC <- all_long_MC %>%
  left_join(replicate_totals_MC, by = "replicate_id") %>%
  mutate(normalised_reads = reads / replicate_total_reads)

##====1. Aggregate normalized reads per replicate/fasta_id====
all_long_norm_agg_MC_fasta <- all_long_norm_MC %>%
  group_by(replicate_id, fasta_id) %>%
  summarise(normalised_reads = sum(normalised_reads), .groups = "drop")

##====2. Pivot normalized data to wide format====
all_wide_norm_MC_fasta <- all_long_norm_agg_MC_fasta %>%
  pivot_wider(names_from = fasta_id, values_from = normalised_reads, values_fill = 0)

###====OPTIONAL: remove K_C controls if needed====
all_wide_norm_MC_fasta <- all_wide_norm_MC_fasta[!grepl("^K_C", all_wide_norm_MC_fasta$replicate_id), ]

##====3. Merge normalized wide matrix with metadata====
merged_data_MC_norm_wide_fasta <- merge(
  metadata_rep_MC[, c("replicate_id", "sample_id", "Site", "SubLocation",
                      "pH", "Temp", "EC", "Region")],
  all_wide_norm_MC_fasta,
  by = "replicate_id"
)

##====4. Omit NA values====

merged_data_MC_norm_wide_fasta <- na.omit(merged_data_MC_norm_wide_fasta)

##====5. Build NMDS dataframe, env, and set up matrices====

nmds_df <- merged_data_MC_norm_wide_fasta
com <- nmds_df %>%
  select(-replicate_id, -sample_id, -Site, -SubLocation, -pH, -Temp, -EC, -Region)
env <- nmds_df %>%
  select(replicate_id, sample_id, Site, SubLocation, pH, Temp, EC, Region)

m_com <- as.matrix(com)
m_com <- m_com[, colSums(m_com, na.rm = TRUE) > 0, drop = FALSE]
# m_com <- decostand(m_com, method = "hellinger")
dist_mat <- vegdist(m_com, method = "raup")

##====6. Run and plot initial NMDS====

set.seed(123)

nmds = metaMDS(m_com, distance = "raup", k=3)

plot(nmds)
stressplot(nmds)
dimcheck_out <- 
  dimcheckMDS(com,
              distance = "raup",
              k = 6)
print(dimcheck_out)

##====7. extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package====
#add columns to data frame 
data.scores_MC <- as.data.frame(scores(nmds, display = "sites"))
data.scores_MC$replicate_id <- env$replicate_id
data.scores_MC$Site <- env$Site
data.scores_MC$SubLocation <- env$SubLocation

ggplot(data.scores_MC, aes(x = NMDS1, y = NMDS2, color = Site, shape = SubLocation)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = SubLocation, color=SubLocation),
               linewidth = 1, alpha = 0.6) + 
  theme_bw() +
  labs(title = "NMDS - Mowrow Creek 12SV5")

head(data.scores_MC)

##====8. Run envfit====

envfit_data_MC <- env %>%
  select(Site, SubLocation, pH, Temp, EC) %>%
  mutate(
    pH = as.numeric(pH),
    Temp = as.numeric(Temp),
    EC = as.numeric(EC),
    Site = as.factor(Site),
    SubLocation = as.factor(SubLocation)
  )

en <- envfit(nmds, envfit_data_MC, permutations = 999, na.rm = TRUE)
en_coord_cont <- as.data.frame(scores(en, "vectors")) * 2
en_coord_cat  <- as.data.frame(scores(en, "factors")) * 2
colnames(en_coord_cont)[1:2] <- c("NMDS1", "NMDS2")
colnames(en_coord_cat)[1:2]  <- c("NMDS1", "NMDS2")

##===9. Plot finalized NMDS====

ggplot(data = data.scores_MC, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color=Site), size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = SubLocation, color=SubLocation),
               linewidth = 1, alpha = 0.6) +
  geom_segment(data = en_coord_cont,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               linewidth = 0.8, alpha = 0.6, colour = "grey30",
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text(data = en_coord_cont,
            aes(x = NMDS1, y = NMDS2,
                label = row.names(en_coord_cont)),
            colour = "grey30", fontface = "bold", vjust = -0.8) +
  labs(shape = "Site", linetype = "Site") +
  theme_minimal(base_size = 14, base_family = "Baskerville") +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 16, colour = "grey30"),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 20, colour = "grey30"),
    legend.key = element_blank(),
    legend.title = element_text(face = "bold", size = 16, colour = "grey30"),
    legend.text = element_text(size = 15, colour = "grey30")
  )

# geom_point(data = en_coord_cat,
#            aes(x = NMDS1, y = NMDS2),
#            shape = 18, size = 4, alpha = 0.7, colour = "navy") +
# geom_polygon(data = data.scores_MC, 
#             aes(x = NMDS1, y = NMDS2, fill = SubLocation, group = SubLocation), 
#             alpha = 0.30) +

#====Raup-Crick from Chase et al., 2016====
# The way I see it, this example is different from where we used vegan to calc
# NMDS plots of raup and jaccard, as here we are using only species assignments
# whereas in vegan we were using only ASV measures and environmental factors. We
# wanted to do this as we want to account for hits that were logged only at the 
# genus or family levels. Here, instead we will just focus on the assignments
# from bayes and BLAST in a p/a matrix.

# Approaching 1: More dissimilar/ more different
# Zero: Dissimilar
# Approaching -1: Less dissimilar/ more similar
# ... than expected by random chance

##====Dataframe setup====

fasta_reads_asvXsite <- diversity_matrix_MC_ASV_taxa %>%
  distinct() %>%
  mutate(replicate_id = sub("_J-Set[0-9]+$", "", sample_id)) %>%
  select(-c(sample_id, Village, SubLocation, Site, Replicate, pH, Temp, EC, species)) %>%
  pivot_wider(names_from = fasta_id,
              values_from = reads,
              values_fn = sum,
              values_fill = 0) %>%
  column_to_rownames("replicate_id") %>%
  as.data.frame()

species_pa_both_spXsite <- both_accum_data_with_reads_and_sub %>%
  distinct() %>%
  select(-c(fasta_id, SubLocation, Quadrat, reads)) %>%
  mutate(Presence = 1) %>%
  pivot_wider(names_from = Species, 
              values_from = Presence, 
              values_fill = 0,
              values_fn = max) %>%
  select(-`NA`) %>%         
  column_to_rownames("Site") %>%
  as.data.frame() 

order_pa_both_spXsite <- both_accum_data_with_reads_and_sub %>%
  distinct() %>%
  left_join(
    both_data_taxon %>%
      select(Species, order) %>%
      distinct(),
    by = "Species"
  ) %>%
  select(-c(fasta_id, SubLocation, Quadrat, reads, Species)) %>%
  mutate(Presence = 1) %>%
  pivot_wider(
    names_from = order,
    values_from = Presence,
    values_fill = 0,
    values_fn = max
  ) %>%
  select(-`NA`) %>%
  column_to_rownames("Site") %>%
  as.data.frame()

family_pa_both_spXsite <- both_accum_data_with_reads_and_sub %>%
  distinct() %>%
  left_join(
    both_data_taxon %>%
      select(Species, family) %>%
      distinct(),
    by = "Species"
  ) %>%
  select(-c(fasta_id, SubLocation, Quadrat, reads, Species)) %>%
  mutate(Presence = 1) %>%
  pivot_wider(
    names_from = family,
    values_from = Presence,
    values_fill = 0,
    values_fn = max
  ) %>%
  select(-`NA`) %>%
  column_to_rownames("Site") %>%
  as.data.frame()

spXsite <- species_pa_both_spXsite
spXsite$Site <- NULL
spXsite$"NA" <- NULL

asvXsite <- fasta_reads_asvXsite
asvXsite$replicate_id <- NULL

##====Raup-crick for spXsite====

Sites <- rownames(spXsite)

raup_crick=function(spXsite, plot_names_in_col1=FALSE, classic_metric=FALSE, 
                    split_ties=TRUE, reps=9999, set_all_species_equal=FALSE,
                    as.distance.matrix=TRUE, report_similarity=FALSE){
  ##expects a species by site matrix for spXsite, with row names for plots, or 
  # optionally plots named in column 1.  By default calculates a modification 
  # of the Raup-Crick metric (standardizing the metric to range from -1 to 1 
  # instead of 0 to 1). Specifying classic_metric=TRUE instead calculates the 
  # original Raup-Crick metric that ranges from 0 to 1. The option split_ties 
  # (defaults to TRUE) adds half of the number of null observations that are 
  # equal to the observed number of shared species to the calculation- this is 
  # highly recommended.  The argument report_similarity defaults to FALSE so the 
  # function reports a dissimilarity (which is appropriate as a measure of beta diversity).  
  # Setting report_similarity=TRUE returns a measure of similarity, as Raup and 
  # Crick originally specified.  If ties are split (as we recommend) the dissimilarity 
  # (default) and similarity (set report_similarity=TRUE) calculations can be 
  # flipped by multiplying by -1 (for our modification, which ranges from -1 to 1) 
  # or by subtracting the metric from 1 (for the classic metric which ranges from 0 to 1). 
  # If ties are not split (and there are ties between the observed and expected shared number of species) 
  # this conversion will not work. The argument reps specifies the number of randomizations
  # (a minimum of 999 is recommended- default is 9999).  set_all_species_equal 
  # weights all species equally in the null model instead of weighting species by frequency of occupancy.  
  
  # Note that the choice of how many plots (rows) to include has a real impact on 
  # the metric, as species and their occurrence frequencies across the set of plots 
  # is used to determine gamma and the frequency with which each species is drawn from the null model	
  
  # ie, we want to use as many "plots" as possible. We'll run by filter replicate then, as it is our base
  # biological replicate.
  
  ##this section moves plot names in column 1 (if specified as being present) 
  # into the row names of the matrix and drops the column of names
#  if(plot_names_in_col1){
#    row.names(spXsite)<-spXsite[,1]
#    spXsite<-spXsite[,-1]
#  }
  
  ## count number of sites and total species richness across all plots (gamma)
  n_sites<-nrow(spXsite)
  gamma<-ncol(spXsite)
  
  ##make the spXsite matrix into a pres/abs. (overwrites initial spXsite matrix):
  ceiling(spXsite/max(spXsite))->spXsite
  
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur<-apply(spXsite, MARGIN=2, FUN=sum)
  
  ##NOT recommended- this is a non-trivial change to the metric:
  ##sets all species to occur with equal frequency in the null model
  ##e.g.- discards any occupancy frequency information
 # if(set_all_species_equal){
 #   occur<-rep(1,gamma)
 # }
  
  ## determine how many unique species richness values are in the dataset
  ##this is used to limit the number of null communities that have to be calculated
  alpha_levels<-sort(unique(apply(spXsite, MARGIN=1, FUN=sum)))
  
  ##make_null:
  ##alpha_table is used as a lookup to help identify which null distribution to 
  # use for the tests later.  It contains one row for each combination of alpha richness levels. 
  
  alpha_table<-data.frame(c(NA), c(NA))
  names(alpha_table)<-c("smaller_alpha", "bigger_alpha")
  col_count<-1
  
  ##null_array will hold the actual null distribution values.  Each element of 
  # the array corresponds to a null distribution for each combination of alpha values.  
  # The alpha_table is used to point to the correct null distribution- the row 
  # numbers of alpha_table correspond to the [[x]] indices of the null_array.  
  # Later the function will find the row of alpha_table with the right combination 
  # of alpha values.  That row number is used to identify the element of null_array 
  # that contains the correct null distribution for that combination of alpha levels. 
  
  null_array<-list()
  
  ##looping over each combination of alpha levels:
  
  for(a1 in 1:length(alpha_levels)){
    for(a2 in a1:length(alpha_levels)){
      ##build a null distribution of the number of shared species for a pair of alpha values:
      null_shared_spp<-NULL
      for(i in 1:reps){
        ##two empty null communities of size gamma:
        com1<-rep(0,gamma)
        com2<-rep(0,gamma)
        ##add alpha1 number of species to com1, weighting by species occurrence frequencies:
        com1[sample(1:gamma, alpha_levels[a1], replace=FALSE, prob=occur)]<-1
        ##same for com2:
        com2[sample(1:gamma, alpha_levels[a2], replace=FALSE, prob=occur)]<-1
        ##how many species are shared in common?
        null_shared_spp[i]<-sum((com1+com2)>1)
      }
      
      ##store null distribution, record values for alpha 1 and 2 in the alpha_table 
      # to help find the correct null distribution later:
      null_array[[col_count]]<-null_shared_spp
      alpha_table[col_count, which(names(alpha_table)=="smaller_alpha")]<-alpha_levels[a1]
      alpha_table[col_count, which(names(alpha_table)=="bigger_alpha")]<-alpha_levels[a2]
      
      #increment the counter for the columns of the alpha table/ elements of the null array
      col_count<-col_count+1
    }
  }
  
  ##create a new column with both alpha levels to match on:
  alpha_table$matching<-paste(alpha_table[,1], alpha_table[,2], sep="_")
  ##do the test:
  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  results<-matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))
  
  ##for each pair of sites (duplicates effort now to make a full matrix instead 
  # of a half one- but this part should be minimal time as compared to the null model building)
  for(i in 1:n_sites){
    for(j in 1:n_sites){
      ##how many species are shared between the two sites:
      n_shared_obs<-sum((spXsite[i,]+spXsite[j,])>1)
      ## what was the observed richness of each site?
      obs_a1<-sum(spXsite[i,])
      obs_a2<-sum(spXsite[j,])
      ##place these alphas into an object to match against alpha_table (sort so smaller alpha is first)
      obs_a_pair<-sort(c(obs_a1, obs_a2))
      ##match against the alpha table- row index identifies which element of the null array contains the correct null distribution for the observed combination of alpha values:
      null_index<-which(alpha_table$matching==paste(obs_a_pair[1], obs_a_pair[2], sep="_"))
      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null<-sum(null_array[[null_index]]==n_shared_obs)
      ##how many null values are bigger than the observed value?
      num_greater_in_null<-sum(null_array[[null_index]]>n_shared_obs)
      rc<-(num_greater_in_null)/reps
      if(split_ties){
        rc<-((num_greater_in_null+(num_exact_matching_in_null)/2)/reps)
      }
      if(!classic_metric){
        ##our modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
        rc<-(rc-.5)*2
      }
      ## at this point rc represents an index of dissimilarity- multiply by -1 to convert to a similarity as specified in the original 1979 Raup Crick paper
      if(report_similarity & !classic_metric){
        rc<- rc*-1
      }
      ## the switch to similarity is done differently if the original 0 to 1 range of the metric is used:
      if(report_similarity & classic_metric){
        rc<- 1-rc
      }
      ##store the metric in the results matrix:
      results[i,j]<-round(rc, digits=2)
    }
  }

  if(as.distance.matrix){
    results<-as.dist(results)
  }	

  return(results)
}

rc_results <- raup_crick(spXsite, plot_names_in_col1 = FALSE)

rc_matrix <- as.matrix(rc_results)

rc_melt <- melt(rc_matrix, varnames = c("Site1", "Site2"), value.name = "Raup_Crick")

ggplot(rc_melt_asv, aes(x = Site1, y = Site2, fill = Raup_Crick)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "turquoise", high = "purple", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), name = "Raup-Crick") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(hjust = 1, vjust = 1, size = 12),
        panel.grid = element_blank(),
        text = element_text(family = "Baskerville"),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 10)) +
  labs(title = "Raup-Crick Distance Matrix: ASV Presence/Absence", x = "Filter Replicates", y = "Filter Replicates")

ggplot(rc_melt, aes(x = Site1, y = Site2, fill = Raup_Crick)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "turquoise", high = "purple", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), name = "Raup-Crick") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(hjust = 1, vjust = 1, size = 12),
        panel.grid = element_blank(),
        text = element_text(family = "Baskerville"),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 10)) +
  labs(title = "Raup-Crick Distance Matrix: Species Presence/Absence", x = "Filter Replicates", y = "Filter Replicates")


# Reframe to fit metaMDS format
# Approaching 0: Less dissimilar
# Approaching 0.5: Dissimilar (neutral)
# Approaching 1: More dissimilar

rc_dist <- as.dist((rc_matrix + 1) / 2)

nmds_result_rc <- metaMDS(rc_dist, k = 2, trymax = 100)

stressplot(nmds_result_rc_asv)

metadata_rep_MC_nmds <- metadata_rep_MC %>%
  mutate(replicate_id = gsub("_J", "", replicate_id))

metadata_rep_MC_nmds <- metadata_rep_MC_nmds %>%
  select() %>%
  rename(Site = replicate_id)

# Extract scores and plot

env_rc <- nmds_scores_rc %>%
  select(sample_id, replicate_id, Site, SubLocation, pH, Temp, EC, Region)

nmds_scores_rc <- as.data.frame(scores(nmds_result_rc, display = "sites")) 

nmds_scores_rc$Site <- rownames(nmds_scores_rc)

nmds_scores_rc <- nmds_scores_rc %>%
  rename(replicate_id = Site) %>%
  left_join(metadata_rep_MC_nmds, by = "replicate_id")

envfit_data_rc <- env_rc %>%
  select(Site, replicate_id, SubLocation, pH, Temp, EC) %>%
  mutate(
    pH = as.numeric(pH),
    Temp = as.numeric(Temp),
    EC = as.numeric(EC),
    Site = as.factor(Site),
    SubLocation = as.factor(SubLocation)
  )

en_rc <- envfit(nmds_result_rc, envfit_data_rc, permutations = 999, na.rm = TRUE)
en_coord_cont_rc <- as.data.frame(scores(en_rc, "vectors")) * 2
en_coord_cat_rc  <- as.data.frame(scores(en_rc, "factors")) * 2
colnames(en_coord_cont_rc)[1:2] <- c("NMDS1", "NMDS2")
colnames(en_coord_cat_rc)[1:2]  <- c("NMDS1", "NMDS2")


ggplot(nmds_scores_rc, aes(x = NMDS1, y = NMDS2, label = Site)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.8, size = 3) +
  theme_bw() +
  labs(title = "Raup-Crick NMDS",
       subtitle = paste("Stress:", round(nmds_result_rc$stress, 3)))

ggplot(data = nmds_scores_rc, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color=Site), size = 4, alpha = 0.8) +
  scale_color_manual(values = site_colors) + 
  new_scale_color() +
  stat_ellipse(aes(group = SubLocation, color=SubLocation),
               linewidth = 1, alpha = 0.6) +
  scale_color_manual(values = sublocation_colors, name = "SubLocation") +
  geom_segment(data = en_coord_cont,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               linewidth = 0.8, alpha = 0.6, colour = "grey30",
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text(data = en_coord_cont_rc,
            aes(x = NMDS1, y = NMDS2,
                label = row.names(en_coord_cont_rc)),
            colour = "grey30", fontface = "bold", vjust = -0.8) +
  labs(shape = "Site", linetype = "Site", title = "Raup-Crick based NMDS Evaluation of Species Presence/Absence in Mowrow Creek") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(face = "bold", colour = "grey30"),
    legend.key = element_blank(),
    legend.title = element_text(face = "bold", colour = "grey30"),
    legend.text = element_text(colour = "grey30")
  )

##====Raup-crick for asvXsite====

rc_results_asv <- raup_crick(asvXsite, plot_names_in_col1 = FALSE)

rc_matrix_asv <- as.matrix(rc_results_asv)

rc_melt_asv <- melt(rc_matrix_asv, varnames = c("Site1", "Site2"), value.name = "Raup_Crick")

ggplot(rc_melt_asv, aes(x = Site1, y = Site2, fill = Raup_Crick)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "turquoise", high = "purple", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), name = "Raup-Crick") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid = element_blank()) +
  labs(title = "Raup-Crick Distance Matrix: ASV Abundance", 
       x = "Filter Replicates", y = "Filter Replicates")

# Convert to metaMDS format

rc_dist_asv <- as.dist((rc_matrix_asv + 1) / 2)

nmds_result_rc_asv <- metaMDS(rc_dist_asv, k = 3, trymax = 100)

stressplot(nmds_result_rc_asv)

dimcheckMDS(asvXsite, distance = "raup", trymax = 20)

metadata_rep_MC_nmds <- metadata_rep_MC %>%
  mutate(replicate_id = gsub("_J", "", replicate_id))

metadata_rep_MC_nmds <- metadata_rep_MC_nmds %>%
  select() %>%
  rename(Site = replicate_id)

# Extract scores and plot

nmds_scores_rc_asv <- as.data.frame(scores(nmds_result_rc_asv, display = "sites")) 

nmds_scores_rc_asv$replicate_id <- rownames(nmds_scores_rc_asv)

nmds_scores_rc_asv <- nmds_scores_rc_asv %>%
  left_join(metadata_rep_MC_nmds, by = "replicate_id")

env_rc_asv <- nmds_scores_rc_asv %>%
  select(replicate_id, Site, SubLocation, pH, Temp, EC, Region)

envfit_data_rc_asv <- env_rc_asv %>%
  select(replicate_id, Site, SubLocation, pH, Temp, EC) %>%
  mutate(
    pH = as.numeric(pH),
    Temp = as.numeric(Temp),
    EC = as.numeric(EC),
    Site = as.factor(Site),
    SubLocation = as.factor(SubLocation)
  )

en_rc_asv <- envfit(nmds_result_rc_asv, envfit_data_rc_asv, permutations = 999, na.rm = TRUE)
en_coord_cont_rc_asv <- as.data.frame(scores(en_rc_asv, "vectors")) * 2
en_coord_cat_rc_asv  <- as.data.frame(scores(en_rc_asv, "factors")) * 2
colnames(en_coord_cont_rc_asv)[1:2] <- c("NMDS1", "NMDS2")
colnames(en_coord_cat_rc_asv)[1:2]  <- c("NMDS1", "NMDS2")

ggplot(nmds_scores_rc_asv, aes(x = NMDS1, y = NMDS2, label = Site)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.8, size = 3) +
  theme_bw() +
  labs(title = "Raup-Crick NMDS",
       subtitle = paste("Stress:", round(nmds_result_rc_asv$stress, 3)))

library(ggnewscale)

sublocation_colors <- c("Lower" = "#E69F00", "Middle" = "#56B4E9", "Upper" = "#009E73")

ggplot(data = nmds_scores_rc_asv, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color=Site), size = 4, alpha = 0.8) +
  scale_color_manual(values = site_colors) + 
  new_scale_color() +
  stat_ellipse(aes(group = SubLocation, color=SubLocation),
               linewidth = 1, alpha = 0.6) +
  scale_color_manual(values = sublocation_colors, name = "SubLocation") +
  geom_segment(data = en_coord_cont,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               linewidth = 0.8, alpha = 0.6, colour = "grey30",
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text(data = en_coord_cont_rc,
            aes(x = NMDS1, y = NMDS2,
                label = row.names(en_coord_cont_rc)),
            colour = "grey30", fontface = "bold", vjust = -0.8) +
  labs(shape = "Site", linetype = "Site", title = "Raup-Crick based NMDS Evaluation of ASV Abundance in Mowrow Creek") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(face = "bold", colour = "grey30"),
    legend.key = element_blank(),
    legend.title = element_text(face = "bold", colour = "grey30"),
    legend.text = element_text(colour = "grey30")
  )
  

#----MaAsLin and iNEXT-----
library(Maaslin2)
library(iNEXT)

# We are doing an ASV level assessment of read count data

BLAST_bayes_merged <- BLAST_bayes_merged %>%
  rename(Species = species)

both_data_taxon <- both_accum_data_with_reads_and_sub %>%
  left_join(BLAST_bayes_merged, by = "Species")

both_data_order <- both_data_taxon %>%
  select(order, reads, Site) %>%
  na.omit(both_data_order)

both_data_family <- both_data_taxon %>%
  select(family, reads, Site) %>%
  na.omit(both_data_family)

both_data_species <- both_data_taxon %>%
  select(Species, reads, Site) %>%
  na.omit(both_data_species)

# Now we transform each long table into a wide format

both_data_order_wide <- both_data_order %>%
  distinct() %>%
  pivot_wider(names_from = order, 
              values_from = reads, 
              values_fill = 0,
              values_fn = max) %>%
  column_to_rownames("Site") %>%
  as.data.frame()   

both_data_family_wide <- both_data_family %>%
  distinct() %>%
  pivot_wider(names_from = family, 
              values_from = reads, 
              values_fill = 0,
              values_fn = max) %>%
  column_to_rownames("Site") %>%
  as.data.frame()   

both_data_species_wide <- both_data_species %>%
  distinct() %>%
  pivot_wider(names_from = Species, 
              values_from = reads, 
              values_fill = 0,
              values_fn = max) %>%
  column_to_rownames("Site") %>%
  as.data.frame()   



metadata_rep_MC_maaslin <- metadata_rep_MC_nmds %>%
  as.data.frame(metadata_rep_MC_maaslin)

# GOALS:
# We want to break down regional differences in taxonomy
# in order for taxonomy to work, we have to filter out any asvs 
# that do not have an assignment to the given level
# We want to run RC with read-count data on an ASV-level basis
# Dissect with MaAsLin.
# We want our diversity matrix, non-normalized

# join tax lookup for both
# filter one out to order level
# filter another out to the family level for each order
# of the two or three dominant orders, filter out to the genera level

# Incorporate into thesis updated alpha and beta div. plots
# Spec accum. curve in discussion along with heatmap for beta-div
# Discuss differences between using ASV or species-level monitoring approaches
# Dig into taxonomic differences between

rownames(metadata_rep_MC_maaslin) <- metadata_rep_MC_maaslin$Site

metadata_rep_MC_maaslin$Site <- NULL

fit_data = Maaslin2(
  input_data = both_data_species_wide, 
  input_metadata = metadata_rep_MC_maaslin, 
  output = "species_reads_sublocation", 
  fixed_effects = c("SubLocation"),
  reference = c("SubLocation,Middle"))

fit_data = Maaslin2(
  input_data = species_pa_both_spXsite, 
  input_metadata = metadata_rep_MC_maaslin, 
  output = "species_pa_sublocation", 
  fixed_effects = c("SubLocation"),
  reference = c("SubLocation,Middle"))

fit_data = Maaslin2(
  input_data = order_pa_both_spXsite, 
  input_metadata = metadata_rep_MC_maaslin, 
  output = "order_pa_sublocation", 
  fixed_effects = c("SubLocation"),
  reference = c("SubLocation,Middle"))

fit_data = Maaslin2(
  input_data = family_pa_both_spXsite, 
  input_metadata = metadata_rep_MC_maaslin, 
  output = "family_pa_sublocation", 
  fixed_effects = c("SubLocation"),
  reference = c("SubLocation,Middle"))

# Using Maaslin2 for Order or family level-taxonomy would be good to add 
# to discussion!!!

#----ANOSIM for beta-div----
# use same matrices as nmds

nmds_df <- merged_data_MC_norm_wide_fasta
com <- nmds_df %>%
  select(-replicate_id, -sample_id, -Site, -SubLocation, -pH, -Temp, -EC, -Region)
env <- nmds_df %>%
  select(replicate_id, sample_id, Site, SubLocation, pH, Temp, EC, Region)

m_com <- as.matrix(com)
m_com <- m_com[, colSums(m_com, na.rm = TRUE) > 0, drop = FALSE]
# m_com <- decostand(m_com, method = "hellinger")

dist_mat <- vegdist(m_com, method = "raup")

ano = anosim(m_com, nmds_df$SubLocation, distance = "raup", permutations = 9999)
ano

#Call:
#  anosim(x = m_com, grouping = nmds_df$SubLocation, permutations = 9999,      distance = "hellinger") 
# Dissimilarity: hellinger 

#ANOSIM statistic R: 0.431 
#Significance: 1e-04 

#Permutation: free
#Number of permutations: 9999

##----Test for beta-dispersion----

bd <- betadisper(vegdist(m_com, method = "raup", binary = TRUE), nmds_df$SubLocation)
anova(bd)

permutest(bd)

plot(bd)

boxplot(bd)

##########   ALPHA DIVERSITY METRICS ################
# Divide up normalized reads into each upper, middle and lower set
# We will be using merged_data_MC_norm, as this is a replicate-normalized, wide format
# dataset with metadata appended.

merged_data_MC_Upper_norm_fasta <- 
  merged_data_MC_norm_wide_fasta  %>%
  filter(SubLocation == "Upper")

merged_data_MC_Middle_norm_fasta <- 
  merged_data_MC_norm_wide_fasta  %>%
  filter(SubLocation == "Middle")

merged_data_MC_Lower_norm_fasta <- 
  merged_data_MC_norm_wide_fasta  %>%
  filter(SubLocation == "Lower")

##====Shannon Index====
shannon_index_MC_Upper_fasta<-diversity(merged_data_MC_Upper_norm_fasta[,9:2249], index = "shannon")
shannon_index_MC_Middle_fasta<-diversity(merged_data_MC_Middle_norm_fasta[,9:2249], index = "shannon")
shannon_index_MC_Lower_fasta<-diversity(merged_data_MC_Lower_norm_fasta[,9:2249], index = "shannon")

shannon_index_MC_Upper_fasta<-data.frame(
  replicate_id = merged_data_MC_Upper_norm_fasta$replicate_id,
  shannon_index_MC_Upper_fasta = shannon_index_MC_Upper_fasta)

shannon_index_MC_Middle_fasta<-data.frame(
  replicate_id = merged_data_MC_Middle_norm_fasta$replicate_id,
  shannon_index_MC_Middle_fasta = shannon_index_MC_Middle_fasta)

shannon_index_MC_Lower_fasta<-data.frame(
  replicate_id = merged_data_MC_Lower_norm_fasta$replicate_id,
  shannon_index_MC_Lower_fasta = shannon_index_MC_Lower_fasta)

##====Simpson Index====

simpson_index_MC_Upper_fasta<-diversity(merged_data_MC_Upper_norm_fasta[,9:2249], index = "simpson")
simpson_index_MC_Middle_fasta<-diversity(merged_data_MC_Middle_norm_fasta[,9:2249], index = "simpson")
simpson_index_MC_Lower_fasta<-diversity(merged_data_MC_Lower_norm_fasta[,9:2249], index = "simpson")

simpson_index_MC_Upper_fasta<-data.frame(
  replicate_id = merged_data_MC_Upper_norm_fasta$replicate_id,
  simpson_index_MC_Upper_fasta = simpson_index_MC_Upper_fasta)

simpson_index_MC_Middle_fasta<-data.frame(
  replicate_id = merged_data_MC_Middle_norm_fasta$replicate_id,
  simpson_index_MC_Middle_fasta = simpson_index_MC_Middle_fasta)

simpson_index_MC_Lower_fasta<-data.frame(
  replicate_id = merged_data_MC_Lower_norm_fasta$replicate_id,
  simpson_index_MC_Lower_fasta = simpson_index_MC_Lower_fasta)

##====FASTA Counts====
fasta_count_upper <- rowSums(merged_data_MC_Upper_norm_fasta[,9:2249] >0)
fasta_count_middle <- rowSums(merged_data_MC_Middle_norm_fasta[,9:2249] >0)
fasta_count_lower <- rowSums(merged_data_MC_Lower_norm_fasta[,9:2249] >0)

fasta_count_upper<-data.frame(
  replicate_id = merged_data_MC_Upper_norm_fasta$replicate_id,
  fasta_count_upper = fasta_count_upper)

fasta_count_middle<-data.frame(
  replicate_id = merged_data_MC_Middle_norm_fasta$replicate_id,
  fasta_count_middle = fasta_count_middle)

fasta_count_lower<-data.frame(
  replicate_id = merged_data_MC_Lower_norm_fasta$replicate_id,
  fasta_count_lower = fasta_count_lower)

##====Merge Shannon, Simpson, and fasta counts====
alpha_12SV5_MC_Upper_fasta<-merge(shannon_index_MC_Upper_fasta,simpson_index_MC_Upper_fasta,by="replicate_id")
alpha_12SV5_MC_Middle_fasta<-merge(shannon_index_MC_Middle_fasta,simpson_index_MC_Middle_fasta,by="replicate_id")
alpha_12SV5_MC_Lower_fasta<-merge(shannon_index_MC_Lower_fasta,simpson_index_MC_Lower_fasta,by="replicate_id")

alpha_12SV5_MC_Upper_fasta<-merge(alpha_12SV5_MC_Upper_fasta,fasta_count_upper,by="replicate_id")
alpha_12SV5_MC_Middle_fasta<-merge(alpha_12SV5_MC_Middle_fasta,fasta_count_middle,by="replicate_id")
alpha_12SV5_MC_Lower_fasta<-merge(alpha_12SV5_MC_Lower_fasta,fasta_count_lower,by="replicate_id")

write.csv(alpha_12SV5_MC_Upper,"alpha_12SV5_MC_Upper.csv")
write.csv(alpha_12SV5_MC_Middle,"alpha_12SV5_MC_Middle.csv")
write.csv(alpha_12SV5_MC_Lower,"alpha_12SV5_MC_Lower.csv")

##====rename columns====
alpha_12SV5_MC_Lower_fasta <- rename(alpha_12SV5_MC_Lower_fasta,shannon = shannon_index_MC_Lower_fasta,
                                     simpson = simpson_index_MC_Lower_fasta, fasta_count = fasta_count_lower)

alpha_12SV5_MC_Middle_fasta <- rename(alpha_12SV5_MC_Middle_fasta,shannon = shannon_index_MC_Middle_fasta,
                                      simpson = simpson_index_MC_Middle_fasta, fasta_count = fasta_count_middle)

alpha_12SV5_MC_Upper_fasta <- rename(alpha_12SV5_MC_Upper_fasta,shannon = shannon_index_MC_Upper_fasta,
                                     simpson = simpson_index_MC_Upper_fasta, fasta_count = fasta_count_upper)

alpha_12SV5_MC_All_fasta <- rbind(alpha_12SV5_MC_Lower_fasta, alpha_12SV5_MC_Middle_fasta, alpha_12SV5_MC_Upper_fasta)

##====Attach metadata====
alpha_12SV5_MC_All_metadata_fasta <- merge(metadata_rep_MC[, c(
  "replicate_id","Site","Replicate","SubLocation","pH", "Temp", "EC")],
  alpha_12SV5_MC_All_fasta, by="replicate_id")

alpha_12SV5_MC_All_metadata_fasta <- alpha_12SV5_MC_All_metadata_fasta %>%
  mutate(
    richness = fasta_count,
    evenness = ifelse(fasta_count > 0, shannon/log(fasta_count), NA)
  )

#====Alpha-div tests by SubLocation====

##====Kruskal Test by SubLocation====
kruskal.test(shannon ~SubLocation, data = alpha_12SV5_MC_All_metadata_fasta)
kruskal.test(simpson ~SubLocation, data = alpha_12SV5_MC_All_metadata_fasta)
kruskal.test(richness ~SubLocation, data = alpha_12SV5_MC_All_metadata_fasta)
kruskal.test(evenness ~SubLocation, data = alpha_12SV5_MC_All_metadata_fasta)

##====Pairwise Wilcox by SubLocation====
pairwise.wilcox.test(
  alpha_12SV5_MC_All_metadata_fasta$shannon,
  alpha_12SV5_MC_All_metadata_fasta$SubLocation,
  p.adjust.method = "BH"
)

pairwise.wilcox.test(
  alpha_12SV5_MC_All_metadata_fasta$simpson,
  alpha_12SV5_MC_All_metadata_fasta$SubLocation,
  p.adjust.method = "BH"
)

pairwise.wilcox.test(
  alpha_12SV5_MC_All_metadata_fasta$richness,
  alpha_12SV5_MC_All_metadata_fasta$SubLocation,
  p.adjust.method = "BH"
)

pairwise.wilcox.test(
  alpha_12SV5_MC_All_metadata_fasta$evenness,
  alpha_12SV5_MC_All_metadata_fasta$SubLocation,
  p.adjust.method = "BH"
)

##====Plot alpha-div metrics by SubLocation====

ggplot(alpha_12SV5_MC_All_metadata_fasta,
       aes(x = Site, y = shannon, fill = SubLocation))+
  geom_boxplot() + 
  theme_bw()

ggplot(alpha_12SV5_MC_All_metadata_fasta,
       aes(x = Site, y = simpson, fill = SubLocation))+
  geom_boxplot() + 
  theme_bw()

ggplot(alpha_12SV5_MC_All_metadata_fasta,
       aes(x = SubLocation, y = richness, fill = SubLocation))+
  geom_boxplot() + 
  theme_bw()

ggplot(alpha_12SV5_MC_All_metadata_fasta,
       aes(x = SubLocation, y = evenness, fill = SubLocation))+
  geom_boxplot() + 
  theme_bw()

#====Alpha-div tests by Site====

##====Kruskal Test by Site====

kruskal.test(shannon ~Site, data = alpha_12SV5_MC_All_metadata_fasta)
kruskal.test(simpson ~Site, data = alpha_12SV5_MC_All_metadata_fasta)
kruskal.test(richness ~Site, data = alpha_12SV5_MC_All_metadata_fasta)
kruskal.test(evenness ~Site, data = alpha_12SV5_MC_All_metadata_fasta)

##====Pairwise Wilcox Test by Site====

pairwise.wilcox.test(
  alpha_12SV5_MC_All_metadata_fasta$shannon,
  alpha_12SV5_MC_All_metadata_fasta$Site,
  p.adjust.method = "BH"
)

pairwise.wilcox.test(
  alpha_12SV5_MC_All_metadata_fasta$simpson,
  alpha_12SV5_MC_All_metadata_fasta$Site,
  p.adjust.method = "BH"
)

pairwise.wilcox.test(
  alpha_12SV5_MC_All_metadata_fasta$evenness,
  alpha_12SV5_MC_All_metadata_fasta$Site,
  p.adjust.method = "BH"
)

pairwise.wilcox.test(
  alpha_12SV5_MC_All_metadata_fasta$richness,
  alpha_12SV5_MC_All_metadata_fasta$Site,
  p.adjust.method = "BH"
)

##====Plot alpha-div metrics by Site====

ggplot(alpha_12SV5_MC_All_metadata_fasta,
       aes(x = Site, y = shannon, color = Site, fill = Site))+
  geom_boxplot() + 
  theme_bw() + 
  labs(
    title = "Mowrow Creek Shannon Alpha Diversity by Site",
    x = "Site",
    y = "Shannon Diversity Index"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(alpha_12SV5_MC_All_metadata_fasta,
       aes(x = Site, y = simpson, color = Site, fill = Site))+
  geom_boxplot() + 
  theme_bw() + 
  labs(
    title = "Mowrow Creek Simpson Alpha Diversity by Site",
    x = "Site",
    y = "Simpson Diversity Index"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(alpha_12SV5_MC_All_metadata_fasta,
       aes(x = Site, y = richness, fill = Site))+
  geom_boxplot() + 
  theme_bw()

ggplot(alpha_12SV5_MC_All_metadata_fasta,
       aes(x = Site, y = evenness, fill = Site))+
  geom_boxplot() + 
  theme_bw()

#====Big plots for alpha-div====

ggplot(alpha_12SV5_MC_All_metadata_fasta, 
       aes(x=SubLocation, y=shannon, color = Site, fill = Site)) +
  scale_fill_manual(values = site_colors) +
  scale_color_manual(values = site_colors) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  theme_bw()+
  labs(
    x = "SubLocation",
    y = "Shannon Diversity Index"
  ) +
  theme_bw(base_family = "Baskerville") +
  theme(
    axis.text.x = element_text(angle = 45, size = 14, hjust = 1),
    axis.text.y = element_text(size = 16, hjust = 1),
    axis.title = element_text(size = 20, face = "plain"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    legend.key.size = unit(0.5, "cm"),
    panel.grid.major.x = element_blank()
  )

ggplot(alpha_12SV5_MC_All_metadata_fasta, 
       aes(x=SubLocation, y=simpson, color=Site, fill = Site)) +
  scale_fill_manual(values = site_colors) +
  scale_color_manual(values = site_colors) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  theme_bw()+
  labs(
    x = "SubLocation",
    y = "Simpson's Diversity Index"
  ) +
  theme_bw(base_family = "Baskerville") +
  theme(
    axis.text.x = element_text(angle = 45, size = 14, hjust = 1),
    axis.text.y = element_text(size = 16, hjust = 1),
    axis.title = element_text(size = 20, face = "plain"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    legend.key.size = unit(0.5, "cm"),
    panel.grid.major.x = element_blank()
  )

ggplot(alpha_12SV5_MC_All_metadata_fasta, 
       aes(x=SubLocation, y=richness, color=Site, fill = Site)) +
  scale_fill_manual(values = site_colors) +
  scale_color_manual(values = site_colors) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  theme_bw()+
  labs(
    x = "SubLocation",
    y = "Richness (# of Unique ASVs)"
  ) +
  theme_bw(base_family = "Baskerville") +
  theme(
    axis.text.x = element_text(angle = 45, size = 14, hjust = 1),
    axis.text.y = element_text(size = 16, hjust = 1),
    axis.title = element_text(size = 20, face = "plain"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    legend.key.size = unit(0.5, "cm"),
    panel.grid.major.x = element_blank()
  )

ggplot(alpha_12SV5_MC_All_metadata_fasta, 
       aes(x=SubLocation, y=evenness, color=Site, fill = Site)) +
  scale_fill_manual(values = site_colors) +
  scale_color_manual(values = site_colors) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  theme_bw()+
  labs(
    x = "SubLocation",
    y = "Evenness"
  ) +
  theme_bw(base_family = "Baskerville") +
  theme(
    axis.text.x = element_text(angle = 45, size = 14, hjust = 1),
    axis.text.y = element_text(size = 16, hjust = 1),
    axis.title = element_text(size = 20, face = "plain"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    legend.key.size = unit(0.5, "cm"),
    panel.grid.major.x = element_blank()
  )

#====Rarefaction Curve ====

asv_mat_MC <- seqtab_12SV5_MC_long_postdecon_filtered %>%
  select(fasta_id, sample_id, reads) %>%
  tidyr::pivot_wider(
    names_from = fasta_id,
    values_from = reads,
    values_fill = 0
  ) %>%
  column_to_rownames("sample_id") %>%
  as.matrix()

mode(asv_mat_MC) <- "numeric"

par(
  mar = c(5, 5, 2, 2),
  las = 1,
  cex.axis = 1.1,
  cex.lab = 1.2
)

library(purrr)
library(scales)

get_rarefaction_df <- function(counts, sample_name, step = 100) {
  counts <- as.numeric(counts)
  total_reads <- sum(counts, na.rm = TRUE)
  
  if (total_reads <= 0) {
    return(NULL)
  }
  
  if (total_reads < step) {
    depths <- unique(c(1, total_reads))
  } else {
    depths <- unique(c(1, seq(step, total_reads, by = step), total_reads))
  }
  
  richness <- sapply(depths, function(d) {
    vegan::rarefy(matrix(counts, nrow = 1), sample = d)
  })
  
  data.frame(
    sample_id = sample_name,
    reads = depths,
    richness = as.numeric(richness)
  )
}

rare_df_MC <- purrr::map_dfr(
  seq_len(nrow(asv_mat_MC)),
  function(i) {
    get_rarefaction_df(asv_mat_MC[i, ], rownames(asv_mat_MC)[i], step = 100)
  }
)

meta_MC <- diversity_matrix_MC_ASV_taxa_normalized %>%
  select(sample_id, replicate_id, Site) %>%
  distinct(sample_id, .keep_all = TRUE)

rare_df_MC <- rare_df_MC %>%
  left_join(meta_MC, by = "sample_id")

ggplot(rare_df_MC, aes(x = reads, y = richness, group = sample_id, color = Site)) +
  geom_line(alpha = 0.5, linewidth = 0.6) +
  scale_fill_manual(values = site_colors) +
  scale_x_log10(labels = comma) +
  labs(
    title = "Rarefaction Curves - Mowrow Creek 12SV5",
    x = "Sequencing depth (reads, log scale)",
    y = "Observed ASV richness"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )

#----Chao Richness Estimates ----

chao_vals_MC <- t(estimateR(asv_mat_MC)) 
head(chao_vals_MC)

chao_df_MC <- as.data.frame(t(estimateR(asv_mat_MC)))
chao_df_MC$sample_id <- rownames(chao_df_MC)
head(chao_df_MC)

chao_meta_MC <- merge(
  metadata_rep_MC[, c("replicate_id", "sample_id", "Site")],
  chao_df_MC,
  by="sample_id"
)


#----Hill Number Diviersity Estimates----




#====TAXONOMY====

# Here, we can run a BLAST version, and Sophie's version. For now, I will 
# prioritize Sophie's classifier. We can run BLAST version later and stitch it in.

MC_taxa <- tax_lookup_MC

#----Sophie's Classifier----

asv_bayes_vsearch <- JackASVs_12SV5_vsearch_bayesian
vsearch100_tax <- vsearch100_taxonomy_Sadler_12SV5_final
bayesian_tax <- bayesian_taxonomy_Sadler_12SV5_final

##----Remove species----

remove_species_tax_v2 <- function(df, species_vec) {
  
  removed_df_v2 <- df %>%
    filter(S %in% species_vec)
  
  filtered_df_v2 <- df %>%
    filter(!S %in% species_vec)
  
  list(
    filtered = filtered_df_v2,
    removed = removed_df_v2
  )
}

species_to_remove_v2 <- c("s__sapiens", "s__troglodytes", "s__lupus",
                          "s__gallus", "s__cattus", "s__taurus", "s__musculus")

asv_bayes_vsearch_clean<- asv_bayes_vsearch %>%
  remove_species_tax_v2(species_to_remove_v2)

asv_bayes_vsearch_filtered <- asv_bayes_vsearch_clean$filtered

asv_bayes_vsearch_filtered <- asv_bayes_vsearch_filtered %>%
  rename(fasta_id = OTUid)

fasta_id_check<-unique(asv_bayes_vsearch_filtered$fasta_id)

####----Newick Tree edits and redo, v2----

tax_ranks <- c("K", "P", "C", "O", "F", "G", "S")

# ── Step 1: Clean all rank columns ───────────────────────────────────────────
# Strip QIIME-style prefixes, whitespace, and coerce blanks/string "NA" to NA
asv_bayes_vsearch_filtered_MC_meta <- asv_bayes_vsearch_filtered_MC_meta %>%
  mutate(across(all_of(tax_ranks), ~ gsub("[a-z]__", "", .x))) %>%
  mutate(across(all_of(tax_ranks), ~ trimws(.x))) %>%
  mutate(across(all_of(tax_ranks), ~ na_if(.x, ""))) %>%
  mutate(across(all_of(tax_ranks), ~ na_if(.x, "NA")))

# ── Step 1b: Standardize known classifier inconsistencies ────────────────────
asv_bayes_vsearch_filtered_MC_meta <- asv_bayes_vsearch_filtered_MC_meta %>%
  mutate(
    C = case_when(
      C == "Actinopterygii" ~ "Actinopteri",  # standardize to shorter form
      TRUE ~ C
    )
  )

# ── Step 2: Drop species assignments with no genus ────────────────────────────
# Prevents Family → Species jumps in the tree
asv_bayes_vsearch_filtered_MC_meta <- asv_bayes_vsearch_filtered_MC_meta %>%
  mutate(S = if_else(is.na(G), NA_character_, S))

# ── Step 3: Combine genus + species into a binomial name ─────────────────────
# e.g. G = "Staphylococcus", S = "aureus" → S = "Staphylococcus aureus"
asv_bayes_vsearch_filtered_MC_meta <- asv_bayes_vsearch_filtered_MC_meta %>%
  mutate(
    S = case_when(
      !is.na(G) & !is.na(S) ~ paste(G, S),
      TRUE ~ S
    )
  )

# ── Step 3b: Sanitize special characters from tip labels ─────────────────────
# Newick reserved characters: ( ) , ; : [ ] and sometimes spaces
asv_bayes_vsearch_filtered_MC_meta <- asv_bayes_vsearch_filtered_MC_meta %>%
  mutate(
    S = gsub(":", "_", S),      # colons break Newick branch length parsing
    S = gsub(";", "_", S),      # semicolons terminate Newick strings
    S = gsub(",", "_", S),      # commas separate taxa in Newick
    S = gsub("\\(", "_", S),    # parentheses define tree structure
    S = gsub("\\)", "_", S),
    S = gsub("\\[", "_", S),    # square brackets are Newick comments
    S = gsub("\\]", "_", S),
    S = trimws(S)
  )

# ── Step 4: Build pathString ──────────────────────────────────────────────────
# Drops NA/empty ranks per row; paths truncate at deepest known rank
asv_bayes_vsearch_filtered_MC_meta$pathString <- purrr::pmap_chr(
  asv_bayes_vsearch_filtered_MC_meta[, tax_ranks],
  function(...) {
    x <- c(...)
    x <- x[!is.na(x) & x != ""]
    paste(c("Root", x), collapse = "/")
  }
)

# ── Step 4b: Sanitize pathString for Newick compatibility ────────────────────
asv_bayes_vsearch_filtered_MC_meta <- asv_bayes_vsearch_filtered_MC_meta %>%
  mutate(
    pathString = gsub(":", "_", pathString),
    pathString = gsub(";", "_", pathString),
    pathString = gsub(",", "_", pathString),
    pathString = gsub("\\(", "_", pathString),
    pathString = gsub("\\)", "_", pathString),
    pathString = gsub("\\[", "_", pathString),
    pathString = gsub("\\]", "_", pathString)
  )

# ── Step 5: Assign LCA label ──────────────────────────────────────────────────
# Most resolved non-NA rank for each ASV
asv_bayes_vsearch_filtered_MC_meta <- asv_bayes_vsearch_filtered_MC_meta %>%
  mutate(lca_label = coalesce(S, G, F, O, C, P, K))

# ── Step 6: Diagnostic checks — run these before building the tree ────────────

# Are any orders (or other ranks) still duplicated in paths?
asv_bayes_vsearch_filtered_MC_meta %>%
  filter(grepl("Characiformes", pathString)) %>%  # swap in any suspect taxon
  select(pathString) %>%
  distinct()

# Are any full pathStrings duplicated across rows?
asv_bayes_vsearch_filtered_MC_meta %>%
  count(pathString) %>%
  filter(n > 1)

# Quick sanity check on the S column
asv_bayes_vsearch_filtered_MC_meta %>%
  filter(!is.na(S)) %>%
  select(G, S) %>%
  distinct() %>%
  head(20)

# ── Step 7: Build and export tree ─────────────────────────────────────────────
tree_input_MC <- asv_bayes_vsearch_filtered_MC_meta %>%
  distinct(pathString, .keep_all = TRUE) %>%
  mutate(pathString = paste0(pathString, "/", fasta_id))

asv_bayes_vsearch_MC_tree <- as.Node(tree_input_MC)
phy_MC <- as.phylo.Node(asv_bayes_vsearch_MC_tree)

# Relabel tips from ASV IDs back to taxon names

phy_MC$tip.label <- tree_input_MC$lca_label[match(
  gsub("'", "", phy_MC$tip.label),
  tree_input_MC$fasta_id
)]
phy_MC$tip.label <- paste0("'", make.unique(phy_MC$tip.label), "'")

phy_MC$edge.length[is.na(phy_MC$edge.length)] <- 0
write.tree(phy_MC, file = "asv_bayes_vsearch_MC_tree_v1.newick")

# Verify
phy_check_MC <- read.tree("asv_bayes_vsearch_MC_tree_v1.newick")
sum(duplicated(phy_check_MC$tip.label))  # should be 0

#----Number of species MC----

library(ape)
library(treeio)

MC_tree_full <- read.tree("~/Desktop/project_thesis/code_and_notebooks/ASV_Analysis_Perry_Unfiltered/asv_bayes_vsearch_MC_tree_v1.newick")
head(MC_tree_full$tip.label, 20)

library(stringr)

clean_labels <- str_remove_all(MC_tree_full$tip.label, "'")

genera <- str_split_fixed(clean_labels, "_", 2)[, 1]
species <- clean_labels

n_genera <- length(unique(genera))
n_species <- length(unique(clean_labels))

cat("Unique genera:", n_genera, "\n")
cat("Unique species:", n_species, "\n")

genus_summary <- as.data.frame(table(genera))
colnames(genus_summary) <- c("Genus", "Species_Count")
genus_summary[order(-genus_summary$Species_Count), ]  # sorted by richness

#####----Building subtrees----
build_subtree <- function(df, class_filter, terminal_rank, outfile) {
  ranks_to_use <- tax_ranks[1:which(tax_ranks == terminal_rank)]
  
  df %>%
    filter(C == class_filter) %>%
    mutate(pathString = purrr::pmap_chr(
      pick(all_of(ranks_to_use)),
      function(...) {
        x <- c(...)
        x <- x[!is.na(x) & x != ""]
        paste(c("Root", x), collapse = "/")
      }
    )) %>%
    distinct(pathString, .keep_all = TRUE) %>%
    as.Node() %>%
    as.phylo.Node() %>%
    { .$edge.length[is.na(.$edge.length)] <- 0; . } %>%
    write.tree(file = outfile)
}

build_subtree(asv_bayes_vsearch_filtered_MC_meta, 
              class_filter  = "Actinopteri", 
              terminal_rank = "S", 
              outfile       = "fish_species_tree_MC.newick")

build_subtree(asv_bayes_vsearch_filtered_MC_meta, 
              class_filter  = "Actinopteri", 
              terminal_rank = "G", 
              outfile       = "fish_genus_tree_MC.newick")

build_subtree(asv_bayes_vsearch_filtered_MC_meta, 
              class_filter  = "Mammalia", 
              terminal_rank = "S", 
              outfile       = "mammal_species_tree_MC.newick")

build_subtree(asv_bayes_vsearch_filtered_MC_meta, 
              class_filter  = "Mammalia", 
              terminal_rank = "G", 
              outfile       = "mammal_genus_tree_MC.newick")

#----Venn Diagram of Species----

##----MowrowCreek----

asv_both_MC_venn_species <- both_accum_data_with_reads_and_sub

asv_both_MC_venn_species <- na.omit(asv_both_MC_venn_species)

asv_both_MC_venn_species[, c("Site", "fasta_id", "reads", "Quadrat")]<-NULL 

asv_both_MC_venn_species<-
  unique(asv_both_MC_venn_species)

###----Generate Venn Data and Plot----

upper_Venn<-asv_both_MC_venn_species %>% 
  filter(SubLocation == "Upper")
lower_Venn<-asv_both_MC_venn_species %>% 
  filter(SubLocation == "Lower")
middle_Venn<-asv_both_MC_venn_species %>% 
  filter(SubLocation == "Middle")

upper  <- upper_Venn$Species
lower  <- lower_Venn$Species
middle <- middle_Venn$Species

UL  <- intersect(upper, lower)
UM  <- intersect(upper, middle)
LM  <- intersect(lower, middle)
ALL <- Reduce(intersect, list(upper, lower, middle))

length(UL)
length(UM)
length(LM)
length(ALL)

UL_only <- setdiff(UL, middle)
UM_only <- setdiff(UM, lower)
LM_only <- setdiff(LM, upper)

length(UL_only)
length(UM_only)
length(LM_only)

lower_only <- setdiff(lower, union(middle, upper))
middle_only <- setdiff(middle, union(lower, upper))
upper_only <- setdiff(upper, union(lower, middle))

length(lower_only)
length(upper_only)
length(middle_only)

# common_all_MC <- semi_join(upper_Venn, lower_Venn, middle_Venn, by = "Species")
# common_UL_MC <- semi_join(upper_Venn, lower_Venn, by = "Species")
# common_UM_MC <- semi_join(upper_Venn, middle_Venn, by = "Species")
# common_LM_MC <- semi_join(lower_Venn, middle_Venn, by = "Species")

VennDiag <- euler(c("Lower" = length(lower_only),
                    "Middle" = length(middle_only),
                    "Upper" = length(upper_only),
                    "Lower&Middle" = length(LM_only),
                    "Lower&Upper" = length(UL_only), 
                    "Middle&Upper" = length(UM_only),
                    "Lower&Middle&Upper" = length(ALL)))

plot(
  VennDiag,
  quantities = TRUE,
  fills = sublocation_colors,
  alpha = 0.5,
  labels = list(font = 2, cex = 1),
  main = "Unique and Shared Species Counts"),
  quantities_args = list(cex = 1))

library(gridExtra)
library(VennDiagram)

# Euler plot
euler_plot <- plot(
  VennDiag,
  quantities = TRUE,
  fills = sublocation_colors,
  alpha = 0.5,
  labels = list(font = 2, cex = 1),
  quantities_args = list(cex = 1)
)

# Triple venn plot
triple_plot <- draw.triple.venn(
  area1 = length(lower),
  area2 = length(middle),
  area3 = length(upper),
  n12 = length(intersect(lower, middle)),
  n13 = length(intersect(lower, upper)),
  n23 = length(intersect(middle, upper)),
  n123 = length(intersect(intersect(lower, middle), upper)),
  category = c("Lower", "Middle", "Upper"),
  fill = sublocation_colors
)

grid.arrange(
  euler_plot,
  ncol = 1
)

#----Heat Map of Species Across Sites----

species_heatmap_df <- species_heatmap_df %>%
  mutate(
    abundance_bin = case_when(
      Reads < 500 ~ "<500",
      Reads < 10000 ~ "<10,000",
      Reads < 50000 ~ "<50,000",
      Reads < 100,000 ~ "<100,000",
      Reads < 500000 ~ "<500,000",
      Reads < 1000000 ~ "<1,000,000",
      Reads < 1500000 ~ "<1,500,000"
      TRUE ~ ">1,500,000"
    ),
    abundance_bin = factor(
      abundance_bin,
      levels = c("<500", "<10,000", "<50,000",
                 "<100,000", "<500,000", "<1,000,000",
                 "<1,500,000", ">1,500,000")
    )
  )

sample_order <- c("SP-01","SP-02","SP-03","T-04","T-05","T-06","OSP-07","OSP-08","OSP-09","OSP-10","OSP-11")

plot_df <- plot_df %>%
  mutate(
    Sample = factor(Sample, levels = sample_order)
  )

plot_df <- plot_df %>%
  group_by(Class, Order) %>%
  mutate(Taxon = fct_inorder(Taxon)) %>%
  ungroup()

#====Comparing Sophie's and Antonio's list to our Classifier Output====

sophie_list <- Sophies_List
qualifiers <- c("cf\\.", "aff\\.", "sp\\.", "nr\\.")

sophie_list$species <- 
  trimws(gsub("\\s*\\(.*\\)$", "", sophie_list$SPECIES))

sophie_list$species <- sapply(strsplit(trimws(sophie_list$species), "\\s+"), 
  function(x) {
  if (length(x) >= 3 && grepl(paste(qualifiers, collapse = "|"), x[2])) {
    paste(x[1:3], collapse = " ")
  } else {
    paste(x[1:2], collapse = " ")
  }
})

sophie_list <- sophie_list %>%
  rename(class = CLASS, order = ORDER, family = FAMILY, species_name_ref = SPECIES, 
         common_names = Common.Names, reference = Reference, notes = Notes, species = species)


antonio_list <- Antonios_Listv2
antonio_list$species <- 
  trimws(gsub("\\s*\\(.*\\)$", "", antonio_list$SPECIES.NAME))

antonio_list <- antonio_list %>%
  rename(order = ORDER, family = FAMILY, species_name_ref = SPECIES.NAME, 
         common_names = COMMON.NAMES, reference = REFERENCE, species = species)

classifier_list <- `MA_Guyana_DB_summary.(1)` %>%
  select(-c(X)) %>%
  rename(class = C, specimen_count = n)



# First, we want to look at comparisons between the lists from Sophie and Antonio
# to identify any missed regions from the classifier.
# Perhaps also worth going through Taphorn and identifying any species
# with "Orinoco", "Barima", or "waini" in notes

# With species identified by Taphorn we will need to species by species to 
# determine locality

##====Overlap between Antonio and Sophie's lists with classifier====

antonio_sophie_overlap <- as.data.frame(intersect(antonio_list$species, sophie_list$species)) %>% 
  rename(species = "intersect(antonio_list$species, sophie_list$species)")

# 124 shared species between Antonio and Sophie's lists

shared_cols <- c("order", "family", "species_name_ref", "common_names", "reference", "species")
antonio_sub <- antonio_list %>% select(all_of(shared_cols))
sophie_sub <- sophie_list %>% select(all_of(shared_cols))

antonio_sophie_merged <- bind_rows(antonio_sub, sophie_sub) %>%
  group_by(species) %>%
  summarise(
    across(everything(), ~ if (n_distinct(na.omit(.x)) > 1) {
      paste(na.omit(unique(.x)), collapse = " | ")
    } else {
      first(na.omit(.x))
    }),
    .groups = "drop"
  )

# 917 species in the merged list between Antonio and Sophie
# This means 4 more were discarded due to mismatching species names

merged_classifier_overlap <- as.data.frame(
  intersect(antonio_sophie_merged$species, classifier_list$species)
) %>% 
  rename(species = "intersect(antonio_sophie_merged$species, classifier_list$species)")

# 209 species were shared between merged A & S lists, and the classifier
# 2,200 species in the classifier. Many of which are not native to Guyana alone

#====Classifier Efficiency====

# Next we want to look at the efficiency of each classifier. 
# Note lack of mock community trials in study.
# Given initial ASV counts, how many of those were converted into a species assignment?
# A genus assignment? Order? Family?
# make it points-based? Out of all ASVs sampled, each assignement receives
# an increasing number of points as closer to species assignment.

#====Species Accumulation Curve=====

# Let's make both a species accumulation and ASV accumulation curve
# Long format table with species, site, region

##====vegan Spec. Accum.====

# in addition to the bayes only spec accum, i'd also like to run one
# with both Bayes and BLAST assignments. I'll need to make that df 
# anyway to run a Raup-crick NMDS on my species assignments.

###====Bayes Spec Accum====

bayes_accum_data_with_reads <- bayes_assignments_postdecon_filtered %>%
  select(-c(K, P, C, O, F, G, lca_label, SubLocation, Village,
            Replicate, pH, EC, Temp, Control)) %>%
  rename(Species = S, Site = sample_id, Quadrat = Site)

bayes_accum_data_with_reads_and_sub <- bayes_assignments_postdecon_filtered %>%
  select(-c(K, P, C, O, F, G, lca_label, Village,
            Replicate, pH, EC, Temp, Control)) %>%
  rename(Species = S, Site = sample_id, Quadrat = Site)

bayes_accum_data <- bayes_accum_data_with_reads %>%
  select(-c(reads, fasta_id)) %>%
  distinct() %>%
  na.omit(bayes_accum_data)

species_accum_bayes <- bayes_accum_data %>%
  distinct() %>%
  mutate(Presence = 1) %>%
  pivot_wider(names_from = Species, 
              values_from = Presence, 
              values_fill = 0) %>% 
  select(-Site, -Quadrat)

head(species_accum_bayes)

plot_accum_bayes <- specaccum(species_accum_bayes, 
                              method = "random")
plot_accum_bayes

plot(plot_accum_bayes, ci.type = "polygon", ci.col = "lightblue",
     main = "Bayes Species Accumulation Curve",
     xlab = "PCR Replicates", ylab = "Species Richness")

specpool(species_accum_bayes)

        # Species chao    chao.se    jack1  jack1.se    jack2     boot  boot.se  n
# All     112    130.9474 9.540851 138.6667 6.491527 146.7008 124.9791 3.513001 81

###====BLAST Spec. Accum====

BLAST_accum_data_with_reads <- BLAST_assignments_species_postdecon_filtered %>%
  rename(Species = sscinames) %>%
  mutate(Site = gsub("_J-Set[0-9]+$", "", Site)) %>%
  select(-c(SubLocation, Village,
            Replicate, pH, EC, Temp, Control)) %>%
  rename(Site = sample_id, Quadrat = Site)

BLAST_accum_data_with_reads_and_sub <- BLAST_assignments_species_postdecon_filtered %>%
  rename(Species = sscinames) %>%
  mutate(Site = gsub("_J-Set[0-9]+$", "", Site)) %>%
  select(-c(Village,
            Replicate, pH, EC, Temp, Control)) %>%
  rename(Site = sample_id, Quadrat = Site)

BLAST_accum_data <- BLAST_accum_data_with_reads %>%
  select(-c(reads, fasta_id)) %>%
  distinct() %>%
  na.omit(BLAST_accum_data)

species_accum_BLAST <- BLAST_accum_data %>%
  distinct() %>%
  mutate(Presence = 1) %>%
  pivot_wider(names_from = Species, 
              values_from = Presence, 
              values_fill = 0) %>% 
  select(-Site, -Quadrat)

plot_accum_BLAST <- specaccum(species_accum_BLAST, 
                             method = "random")
plot_accum_BLAST

plot(plot_accum_BLAST, ci.type = "polygon", ci.col = "lightblue",
     main = "BLAST Species Accumulation Curve",
     xlab = "PCR Replicates", ylab = "Species Richness")

specpool(species_accum_BLAST)

###==== Both Spec Accum ====

both_accum_data_with_reads <- full_join(bayes_accum_data_with_reads, 
                                        BLAST_accum_data_with_reads) %>%
  distinct()

both_accum_data_with_reads_and_sub <- full_join(bayes_accum_data_with_reads_and_sub, 
                                        BLAST_accum_data_with_reads_and_sub) %>%
  distinct()

both_accum_data <- full_join(bayes_accum_data, 
                             BLAST_accum_data) %>%
  distinct()
  
species_accum_both <- both_accum_data %>%
  distinct() %>%
  mutate(Presence = 1) %>%
  pivot_wider(names_from = Species, 
              values_from = Presence, 
              values_fill = 0) %>% 
  select(-Site, -Quadrat)

plot_accum_both <- specaccum(species_accum_both, 
                              method = "random")
plot_accum_both

plot(plot_accum_both, ci.type = "polygon", ci.col = "lightblue",
     main = "Bayes & BLAST Species Accumulation Curve",
     xlab = "PCR Replicates", ylab = "Species Richness")

specpool(species_accum_bayes)

##====Michaelis-Menten====
# Here this code is interchangable for whichever version you are running

 
plot_accum_bayes_df <- data.frame(sites = plot_accum_bayes$sites, 
                                  richness = plot_accum_bayes$richness)
plot_accum_BLAST_df <- data.frame(sites = plot_accum_BLAST$sites, 
                                  richness = plot_accum_BLAST$richness)
plot_accum_bayes_df <- data.frame(sites = plot_accum_both$sites, 
                                  richness = plot_accum_both$richness)


mm_fit_bayes <- nls(richness ~ (a * sites) / (b + sites),
              data = plot_accum_bayes_df,
              start = list(a = max(plot_accum_bayes_df$richness), b = 5))
mm_fit_BLAST <- nls(richness ~ (a * sites) / (b + sites),
                    data = plot_accum_BLAST_df,
                    start = list(a = max(plot_accum_BLAST_df$richness), b = 5))
mm_fit_both <- nls(richness ~ (a * sites) / (b + sites),
                    data = plot_accum_both_df,
                    start = list(a = max(plot_accum_both_df$richness), b = 5))

# Predictions
sites_seq_bayes <- seq(1, max(plot_accum_bayes_df$sites), by = 0.1)
sites_seq_BLAST <- seq(1, max(plot_accum_BLAST_df$sites), by = 0.1)
sites_seq_both <- seq(1, max(plot_accum_both_df$sites), by = 0.1)

mm_pred_bayes <- data.frame(
  sites = sites_seq_bayes,
  richness = predict(mm_fit_bayes, newdata = data.frame(sites = sites_seq_bayes))
)
mm_pred_BLAST <- data.frame(
  sites = sites_seq_BLAST,
  richness = predict(mm_fit_BLAST, newdata = data.frame(sites = sites_seq_BLAST))
)
mm_pred_both <- data.frame(
  sites = sites_seq_both,
  richness = predict(mm_fit_both, newdata = data.frame(sites = sites_seq_both))
)

# Plot original accumulation curve
plot(plot_accum_BLAST, ci.type = "polygon", ci.col = "lightblue",
     main = "Species Accumulation with Michaelis–Menten Fit",
     xlab = "PCR Replicates", ylab = "Species richness") +
  lines(mm_pred$sites, mm_pred$richness, 
        col = "red", lwd = 2) +
  legend("bottomright", legend = "Michaelis–Menten fit", 
         col = "red", lwd = 2)

###=====Clench fit====

clench_fit_bayes <- nls(richness ~ (a * sites) / (1 + b * sites),
                  data = plot_accum_bayes_df,
                  start = list(a = 1, b = 0.1))
clench_fit_BLAST <- nls(richness ~ (a * sites) / (1 + b * sites),
                        data = plot_accum_BLAST_df,
                        start = list(a = 1, b = 0.1))
clench_fit_both <- nls(richness ~ (a * sites) / (1 + b * sites),
                        data = plot_accum_both_df,
                        start = list(a = 1, b = 0.1))

clench_pred_bayes <- data.frame(
  sites = sites_seq_bayes,
  richness = predict(clench_fit_bayes, newdata = data.frame(sites = sites_seq_bayes))
)
clench_pred_BLAST <- data.frame(
  sites = sites_seq_BLAST,
  richness = predict(clench_fit_BLAST, newdata = data.frame(sites = sites_seq_BLAST))
)
clench_pred_both <- data.frame(
  sites = sites_seq_both,
  richness = predict(clench_fit_both, newdata = data.frame(sites = sites_seq_both))
)

# Plot accumulation with Clench fit
plot(plot_accum_BLAST, ci.type = "polygon", ci.col = "lightblue",
     main = "Species Accumulation with Clench Fit",
     xlab = "Sites", ylab = "Species richness") +
  lines(clench_pred$sites, clench_pred$richness, 
        col = "darkgreen", lwd = 2) + 
  legend("bottomright", legend = "Clench fit", 
         col = "darkgreen", lwd = 2)

chao_est_BLAST <- specpool(species_accum_BLAST)$chao
chao_est_bayes <- specpool(species_accum_bayes)$chao
chao_est_both <- specpool(species_accum_both)$chao

ymax = chao_est + 10

plot(plot_accum_BLAST, ci.type = "polygon", ci.col = "lightblue",
     main = "Species Accumulation and Diversity Estimates with BLAST",
     xlab = "Quadrats", ylab = "Species Richness") + 
     ylim(0, ymax) +
  lines(mm_pred$sites, mm_pred$richness, col = "red", lwd = 4) + 
  lines(clench_pred$sites, clench_pred$richness, 
        col = "darkgreen", lwd = 2) + 
  abline(h = chao_est, col = "blue", lty = 2, lwd = 2) + 
  legend("bottomright",
       legend = c("Michaelis–Menten", "Clench", "Chao1"),
       col = c("red", "darkgreen", "blue"),
       lty = c(1, 1, 2), lwd = 2)

accum_df_BLAST <- data.frame(
  sites    = plot_accum_BLAST$sites,
  richness = plot_accum_BLAST$richness,
  sd       = plot_accum_BLAST$sd
)

accum_df_bayes <- data.frame(
  sites    = plot_accum_bayes$sites,
  richness = plot_accum_bayes$richness,
  sd       = plot_accum_bayes$sd
)

accum_df_both <- data.frame(
  sites    = plot_accum_both$sites,
  richness = plot_accum_both$richness,
  sd       = plot_accum_both$sd
)

##-------------Plotting Spec Accum in Vegan-------------
ggplot() +
  
  # --Ribbons
  geom_ribbon(data = accum_df_both,
              aes(x    = sites,
                  ymin = richness - sd,
                  ymax = richness + sd),
              fill  = "lightblue",
              alpha = 0.9) +
  geom_ribbon(data = accum_df_BLAST,
              aes(x    = sites,
                  ymin = richness - sd,
                  ymax = richness + sd),
              fill  = "lightgreen",
              alpha = 0.9) +
  geom_ribbon(data = accum_df_bayes,
              aes(x    = sites,
                  ymin = richness - sd,
                  ymax = richness + sd),
              fill  = "pink",
              alpha = 0.9) +
  
  # Both classifiers
  geom_line(data = accum_df_both,
            aes(x = sites, y = richness)) +
  geom_line(data = mm_pred_both,
            aes(x = sites, y = richness, color = "MM_both"),
            lwd = 1.2) +
  geom_line(data = clench_pred_both,
            aes(x = sites, y = richness, color = "Clench_both"),
            lwd = 0.8) +
  geom_hline(aes(yintercept = chao_est_both, color = "Chao_both"),
             linetype = "dashed",
             lwd = 0.8) +
  
  # Bayes Classifier
  geom_line(data = accum_df_bayes,
            aes(x = sites, y = richness)) +
  geom_line(data = mm_pred_bayes,
            aes(x = sites, y = richness, color = "MM_bayes"),
            lwd = 1.2) +
  geom_line(data = clench_pred_bayes,
            aes(x = sites, y = richness, color = "Clench_bayes"),
            lwd = 0.8) +
  geom_hline(aes(yintercept = chao_est_bayes, color = "Chao_bayes"),
             linetype = "dashed",
             lwd = 0.8) +
  
  # BLAST Classifier
  geom_line(data = accum_df_BLAST,
            aes(x = sites, y = richness)) +
  geom_line(data = mm_pred_BLAST,
            aes(x = sites, y = richness, color = "MM_BLAST"),
            lwd = 1.2) +
  geom_line(data = clench_pred_BLAST,
            aes(x = sites, y = richness, color = "Clench_BLAST"),
            lwd = 0.8) +
  geom_hline(aes(yintercept = chao_est_BLAST, color = "Chao_BLAST"),
             linetype = "dashed",
             lwd = 0.8) +
  
  # Ribbon text labels (placed at right edge of x-axis) 
  annotate("text", x = max(accum_df_both$sites),
           y = max(accum_df_both$richness + accum_df_both$sd),
           label = "Both", hjust = 1.1, vjust = -0.5,
           color = "#2171B5", size = 4, fontface = "bold") +
  
  annotate("text", x = max(accum_df_BLAST$sites),
           y = max(accum_df_BLAST$richness + accum_df_BLAST$sd),
           label = "BLAST", hjust = 1.1, vjust = -0.5,
           color = "#238B45", size = 4, fontface = "bold") +
  
  annotate("text", x = max(accum_df_bayes$sites),
           y = max(accum_df_bayes$richness + accum_df_bayes$sd),
           label = "Bayes", hjust = 1.1, vjust = -0.5,
           color = "#CB181D", size = 4, fontface = "bold") +
  
  # Color values
scale_color_manual(
  name   = "Model / Classifier",
  values = c(
    "MM_BLAST"     = "#1A9E9E",
    "Clench_BLAST" = "#0D5F5F",
    "MM_bayes"     = "#7B2D8B",
    "Clench_bayes" = "#4A1060",
    "MM_both"      = "#D94F00",
    "Clench_both"  = "#8B2F00",
    "Chao_BLAST"  = "#238B45",
    "Chao_bayes"  = "#CB181D",
    "Chao_both"   = "#2171B5"
  ),
  # Optional: reorder legend entries
  breaks = c("MM_BLAST", "Clench_BLAST", "Chao_BLAST",
             "MM_bayes", "Clench_bayes", "Chao_bayes",
             "MM_both",  "Clench_both",  "Chao_both")
) +
  coord_cartesian(ylim = c(0, chao_est_both + 10)) +
  labs(
    title = "Species Accumulation and Diversity Estimates Across Classifiers",
    x     = "Sample Sites (PCR Replicates)",
    y     = "Species Richness (# of Species Observed)"
  ) +
  theme_classic(18) +
  theme(legend.position = "right",
        legend.key.width = unit(1.5, "cm"),   # makes linetype visible in legend
        plot.title       = element_text(size = 14, face = "bold")
  )


##====iNEXT Spec. Accum.====

# For iNEXT we want to generate two files. One with each site 
# as an assemblage, and another with each region as an assemblage
# Any NAs must be converted to "0".

# Note that antonio's list is in French Guianea
# Also, divide up species list into region of MC and 
# Comment on ecological role of each species
# Add per filter rarefaction and gel plots to appendix
# along with code

###====Format to each filter being a df with all species represented in a row====

bayes_accum_data_with_reads$Site = gsub("_J-Set[0-9]+$", "", 
                                        bayes_accum_data_with_reads$Site)
both_accum_data_with_reads_and_sub$Site = gsub("_J-Set[0-9]+$", "", 
                                        both_accum_data_with_reads$Site)

bayes_accum_data_chao <- bayes_accum_data_with_reads %>%
  select(-c(fasta_id, Quadrat)) %>%
  distinct() %>%
  na.omit(bayes_accum_data) %>%
  group_by(Site)

bayes_accum_data_chao_sublocation <- bayes_accum_data_with_reads_and_sub %>%
  select(-c(fasta_id, Quadrat)) %>%
  distinct() %>%
  na.omit(bayes_accum_data) %>%
  group_by(Site)

both_accum_data_chao_sublocation <- both_accum_data_with_reads_and_sub %>%
  select(-c(fasta_id, Quadrat)) %>%
  distinct() %>%
  na.omit(both_accum_data) %>%
  group_by(Site)

bayes_accum_data_chao_list_original <- split(bayes_accum_data_chao, 
                                      bayes_accum_data_chao$Site)

bayes_accum_data_chao_list_sub <- split(bayes_accum_data_chao_sublocation, 
                                             bayes_accum_data_chao_sublocation$SubLocation)

both_accum_data_chao_list_sub <- split(both_accum_data_chao_sublocation, 
                                        both_accum_data_chao_sublocation$SubLocation)

bayes_accum_data_chao_list_original[[1]] %>%
  filter(Species == "Gymnotus carapo") %>%
  select(Species, reads) %>%
  summarise(total = sum(reads, na.rm = TRUE))

bayes_accum_data_chao_list[[1]] %>%
  filter(Species == "Gymnotus carapo")

sum(bayes_accum_data_chao_list_original[[1]]$reads, na.rm = TRUE)

sum(bayes_accum_data_chao_list[[1]]$reads, na.rm = TRUE)

totals_before <- sapply(bayes_accum_data_chao_list_original, function(df) sum(df$reads, na.rm = TRUE))
totals_after <- sapply(bayes_accum_data_chao_list, function(df) sum(df$reads, na.rm = TRUE))

all.equal(totals_before, totals_after)


bayes_accum_data_chao_list <- lapply(bayes_accum_data_chao_list, function(df) {
  df %>%
    ungroup() %>%
    group_by(Species) %>%
    summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop")
})

bayes_accum_data_chao_list_sub <- lapply(bayes_accum_data_chao_list_sub, function(df) {
  df %>%
    ungroup() %>%
    group_by(Species) %>%
    summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop")
})

both_accum_data_chao_list_sub <- lapply(both_accum_data_chao_list_sub, function(df) {
  df %>%
    ungroup() %>%
    group_by(Species) %>%
    summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop")
})

all_species <- unique(unlist(lapply(bayes_accum_data_chao_list, function(df) df$Species)))

bayes_accum_data_chao_list <- lapply(bayes_accum_data_chao_list, function(df) {
  df <- as.data.frame(df)
  missing_species <- setdiff(all_species, df$Species)
  if (length(missing_species) > 0) {
    zero_rows <- data.frame(
      Species = as.character(missing_species),
      reads = 0
    )
    df <- rbind(df, zero_rows)
  }
  df <- df[order(df$Species), ]
  return(df)
})

bayes_accum_data_chao_list_sub <- lapply(bayes_accum_data_chao_list_sub, function(df) {
  df <- as.data.frame(df)
  missing_species <- setdiff(all_species, df$Species)
  if (length(missing_species) > 0) {
    zero_rows <- data.frame(
      Species = as.character(missing_species),
      reads = 0
    )
    df <- rbind(df, zero_rows)
  }
  df <- df[order(df$Species), ]
  return(df)
})

both_accum_data_chao_list_sub <- lapply(both_accum_data_chao_list_sub, function(df) {
  df <- as.data.frame(df)
  missing_species <- setdiff(all_species, df$Species)
  if (length(missing_species) > 0) {
    zero_rows <- data.frame(
      Species = as.character(missing_species),
      reads = 0
    )
    df <- rbind(df, zero_rows)
  }
  df <- df[order(df$Species), ]
  return(df)
})

nrow(bayes_accum_data_chao_list[[1]])
nrow(bayes_accum_data_chao_list[[27]])  

both_accum_data_chao_list_sub_incidence <- lapply(
  both_accum_data_chao_list_sub,
  function(df) {
    df$pa <- ifelse(df$pa > 0, 1, 0)
    df
  }
)

###====Now we can get going on processing====

both_accum_data_chao_list_sub_incidence <- both_accum_data_with_reads_and_sub %>%
  filter(!is.na(Species), Species != "") %>%
  mutate(pa = ifelse(reads > 0, 1, 0)) %>%
  group_by(SubLocation, Species, Site) %>%
  summarise(pa = max(pa), .groups = "drop") %>%
  group_split(SubLocation) %>%
  setNames(unique(both_accum_data_with_reads_and_sub$SubLocation)) %>%
  map(~ .x %>%
        select(Species, Site, pa) %>%
        pivot_wider(
          names_from = Site,
          values_from = pa,
          values_fill = 0
        ) %>%
        column_to_rownames("Species"))

inext_input_abundance_sub <- lapply(both_accum_data_chao_list_sub, function(df) {
  x <- df$reads
  names(x) <- df$Species
  return(x)
})

DataInfo(inext_input_abundance_sub, datatype = "abundance")
DataInfo(inext_input_incidence, datatype = "incidence")

# RUN THIS OVERNIGHT

max_n <- max(sapply(inext_input_abundance_sub$Lower, sum))

# Running this on incidence

out_sub_incidence <- iNEXT(
  both_accum_data_chao_list_sub_incidence,
  q = c(0, 1, 2),
  datatype = "incidence_raw"
)

# Running on abundance

out_sub_abundance_20mil <- iNEXT(inext_input_abundance_sub,
             q = c(0, 1, 2),
             datatype = "abundance",
             endpoint = 20000000,
             knots = 20)

ggiNEXT(out_sub_abundance_4mil, type=1, facet.var="Assemblage") +
  facet_wrap(~ Assemblage, nrow = 9, ncol = 3) +
  xlab("Individual Read Counts")

ggiNEXT(out_sub_incidence, type=1, facet.var="Assemblage") +
  facet_wrap(~ Assemblage, ncol = 3) +
  xlab("Species Counts") +
  theme(text = element_text(family = "Baskerville"),
        panel.spacing.x = unit(0.5, "cm"),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 18, face = "plain")
        )

  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(hjust = 1, size = 14),
        strip.text = element_text(size = 18, face = "plain"),
        axis.title = element_text(size = 20),
        legend.position = "none",
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(face = "bold")) + 

ggiNEXT(out_sub_abundance, type=1, facet.var="Order.q")
# Use below to get a coverage-based sampling curve
ggiNEXT(out_sub_abundance, type=2)
# Then below for coverage-based R/E curves for each site
ggiNEXT(out, type=3, facet.var="Assemblage")

# For richness and incidence we want to use ASVs
# Maybe it would be easier to just use iNEXT for all beta-div calcs? 
# No for now we're just going to use what we used before. But it's an option

ChaoRichness(inext_input, datatype = "abundance", conf = 0.95)
ChaoShannon(inext_input, datatype = "abundance", transform = FALSE, conf = 0.95, B = 200)
ChaoSimpson(inext_input, datatype = "abundance", transform = FALSE, conf = 0.95, B = 200)

# inexT function can use incidence or abundance
# Extrapolation of Hill Numbers with order q

iNEXT(
  inext_input,
  q = 0,
  datatype = "incidence",
  size = NULL,
  endpoint = NULL,
  knots = 40,
  se = TRUE,
  conf = 0.95,
  nboot = 50
)

# Ok here is what we want
#   - PCR read count abundance accumulation curves for each site, and for each region 
#     with enough buffer for extrapolation
#   - Incidence accumulation curves for each site and for each region
#   - If you have extra time after taxonomy, dig into iNEXT too
  

  