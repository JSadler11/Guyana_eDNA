# =============================================================================
# ANALYSIS -- DECONTAMINATION, ABUNDANCE PLOTS,
#             NORMALISATION & PHYLUM-LEVEL COMPOSITION
#
# This section re-loads the filtered CSV and applies the following
# workflow, which is subsequently repeated for each river in WP2C:
#
#  a) Aggregate technical replicates by taking the mean read count per
#     species x site-time combination (not sum, because replicates measure
#     the same community).
#
#  b) Plot a proportional stacked bar chart (position='fill') as a first
#     look at relative amplicon frequency across sites.
#
#  c) Convert to wide format, then run microDecon::decon() for blank-based
#     statistical decontamination. Blank samples are relocated to the front
#     of the table (required by decon()). The function statistically
#     identifies taxa whose abundance pattern mirrors the field/extraction
#     blanks and subtracts them. The number of blanks (numb.blanks=15)
#     and sample counts per event (numb.ind) must match the study design.
#     Samples that end up with zero reads post-decontamination are dropped.
#
#  d) Convert the decontaminated wide table back to long format.
#
#  e) Rarefaction curve code (commented out) is retained for reference.
#     Venn diagram code (commented out) would show ASV overlap by river
#     section and season using VennDiagram/eulerr.
#
#  f) Normalise reads to relative abundance (same approach as Section 8).
#     Column sums (per-ASV totals across all samples) are also computed
#     for use in the 'top taxa' ranking plots.
#
#  g) Generate phylum-level composition plots by splitting the semicolon-
#     delimited taxonomy string into separate columns; X2 holds the phylum
#     rank in the SILVAngs hierarchy. Reads are aggregated to phylum level
#     and plotted as proportional stacked bars, faceted by site.
#     A fixed factor order and manual colour palette ensure consistency
#     across all river plots. Minor phyla are rendered in white.
#
#  h) Export the normalised wide table as NORM_WP2A_ASV_table_wide_filtered.csv
#     -- the direct input to all PERMANOVA, NMDS, CCA, and alpha-diversity
#     analyses in Sections 11-14.
# =============================================================================

library(tidyverse)
library(readr)
library(data.table)
library(vegan)
library(plotly)
library(microDecon)
library(gam)

setwd("~project_thesis/blast/12SV5_BLAST/Full_Outputs")
seqtab_MawrowCreek_12SV5_tax_long_meta_clean_filtered <- read_csv("R_outputs/MawrowCreek_ASV_table_long_filtered.csv")

#aggregating replicates
seqtab_MawrowCreek_12SV5_tax_long_meta_clean_filtered<-seqtab_MawrowCreek_12SV5_tax_long_meta_clean_filtered %>% 
  group_by(SubLocation,sscinames,Control.,sample_id,fasta_id) %>% 
  summarize(reads=mean(reads)) %>%
  rename(SubLocation=SubLocation,sscinames=sscinames,reads=reads)

seqtab_MawrowCreek_12SV5_tax_long_meta_clean_filtered<-seqtab_MawrowCreek_12SV5_tax_long_meta_clean_filtered %>% 
  group_by(SubLocation,sscinames,Control.,sample_id,fasta_id) %>% 
  summarize(reads=mean(reads)) %>%
  rename(SubLocation=SubLocation,sscinames=sscinames, reads=reads)

#RELATIVE AMPLICON FREQUENCY ABUNDANCE
ggplot(seqtab_MawrowCreek_12SV5_tax_long_meta_clean_filtered, aes(x = SubLocation, fill = sscinames, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "white")+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none",
  legend.text = element_text(face = "italic", size = 9))

ggplot(seqtab_MawrowCreek_12SV5_tax_long_meta_clean_filtered, aes(x = SubLocation, fill = sscinames, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "white")+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none",
        legend.text = element_text(face = "italic", size = 9))

#CONVERTING FILTERED LONG FORMAT INTO WIDE FORMAT

seqtab_MawrowCreek_12SV5_tax_wide_filtered <- seqtab_MawrowCreek_12SV5_tax_long_meta_clean_filtered %>%
  group_by(sample_id, fasta_id, sscinames, SubLocation, Control.) %>%
  summarize(reads = sum(reads), .groups = "drop") %>%
  spread(key = sample_id, value = reads, fill = 0)

seqtab_MawrowCreek_12SV5_tax_wide_filtered[] <- lapply(seqtab_MawrowCreek_12SV5_tax_wide_filtered, function(x) type.convert(as.character(x)))

# FILTERING BY NEGATIVE CONTROL
control_samples_MC <- seqtab_MawrowCreek_12SV5_tax_long_meta_clean_filtered %>% 
  filter(Control. == "Y") %>% 
  distinct(sample_id, Sublocation) %>%
  pull(sample_id)

real_samples_MC <- seqtab_MawrowCreek_12SV5_tax_long_meta_clean_filtered %>%
  filter(Control. == "N") %>%
  distinct(sample_id) %>%
  pull(sample_id)

# Decontam
numb.blanks_MC <- length(control_samples_MC)
numb.ind_MC <- seqtab_MawrowCreek_12SV5_tax_long_meta_clean_filtered %>%
  filter(Control. == "N") %>%
  distinct(sample_id, SubLocation) %>%
  count(SubLocation) %>%
  pull(n)

decontaminated_WD <- decon(
  data        = seqtab_MawrowCreek_12SV5_tax_wide_filtered_DECON_WD,
  numb.blanks = numb.blanks_WD,
  numb.ind    = numb.ind_WD,
  taxa        = FALSE
)
reads_removed_WD <- decontaminated_WD$reads.removed

seqtab_MawrowCreek_12SV5_tax_wide_filtered<-decontaminated$decon.table
seqtab_MawrowCreek_12SV5_tax_wide_filtered$Mean.blank<-NULL
seqtab_MawrowCreek_12SV5_tax_wide_filtered<-t(seqtab_MawrowCreek_12SV5_tax_wide_filtered)
seqtab_MawrowCreek_12SV5_tax_wide_filtered<-data.frame(seqtab_MawrowCreek_12SV5_tax_wide_filtered)

header.true <- function(seqtab_MawrowCreek_12SV5_tax_wide_filtered) {
  names(seqtab_MawrowCreek_12SV5_tax_wide_filtered) <- as.character(unlist(seqtab_MawrowCreek_12SV5_tax_wide_filtered[1,]))
  seqtab_MawrowCreek_12SV5_tax_wide_filtered[-1,]
}
seqtab_MawrowCreek_12SV5_tax_wide_filtered<-header.true(seqtab_12SV5_tax_wide_filtered)

i <- c(1:2484) 
seqtab_12SV5_tax_wide_filtered[ , i] <- apply(seqtab_12SV5_tax_wide_filtered[ , i], 2,  
                                              function(x) as.numeric(as.character(x)))

#REMOVE SAMPLES THAT NOW HAVE NO READS IN THEM
seqtab_MawrowCreek_12SV5_tax_wide_filtered<-seqtab_12SV5_tax_wide_filtered[rowSums(seqtab_MawrowCreek_12SV5_tax_wide_filtered[])>0,]

#CONVERTING INTO LONG FORMAT
seqtab_MawrowCreek_12SV5_tax_long_meta_clean_filtered<-seqtab_MawrowCreek_12SV5_tax_wide_filtered
seqtab_MawrowCreek_12SV5_tax_long_meta_clean_filtered <- tibble::rownames_to_column(seqtab_MawrowCreek_12SV5_tax_long_meta_clean_filtered, "rn")
seqtab_MawrowCreek_12SV5_tax_long_meta_clean_filtered <- gather(seqtab_MawrowCreek_12SV5_tax_long_meta_clean_filtered, species, reads, 2:2485, factor_key=TRUE)
seqtab_MawrowCreek_12SV5_tax_long_meta_clean_filtered$reads<-   as.numeric(seqtab_MawrowCreek_12SV5_tax_long_meta_clean_filtered$reads)
seqtab_MawrowCreek_12SV5_tax_long_meta_clean_filtered<-seqtab_MawrowCreek_12SV5_tax_long_meta_clean_filtered %>% rename(SampleSite_time = rn)

write.csv(seqtab_MawrowCreek_12SV5_tax_wide_filtered,"MawrowCreek_ASV_table_wide_filtered.csv")

##SPLITTING BY TAXONOMY
Arthropoda_count <- seqtab_12SV5_tax_long_meta_clean_filtered[grep("Arthropoda", seqtab_12SV5_tax_long_meta_clean_filtered$species), ]
Arthropoda_count<-sum(Arthropoda_count$reads)


#NORMALISING READS
seqtab_12SV5_tax_col_sum2<-rowSums(seqtab_12SV5_tax_wide_filtered[c(1:2484)])
seqtab_12SV5_tax_col_sum2<-data.frame(seqtab_12SV5_tax_col_sum2)
seqtab_12SV5_tax_col_sum2$SampleSite_time <- rownames(seqtab_12SV5_tax_col_sum2)

seqtab_12SV5_tax_normalised<- merge(seqtab_12SV5_tax_col_sum2[, c("SampleSite_time", "seqtab_12SV5_tax_col_sum2")],seqtab_12SV5_tax_long_meta_clean_filtered, by="SampleSite_time")

#GETTING THE SUM OF EACH SAMPLE FOR NORMALISING READS
seqtab_12SV5_tax_col_sum3<-colSums(seqtab_12SV5_tax_wide_filtered[c(1:2484)])
seqtab_12SV5_tax_col_sum3<-data.frame(seqtab_12SV5_tax_col_sum3)

seqtab_12SV5_tax_col_sum3$species <- rownames(seqtab_12SV5_tax_col_sum3)

seqtab_12SV5_tax_normalised<- merge(seqtab_12SV5_tax_col_sum2[, c("SampleSite_time", "seqtab_12SV5_tax_col_sum2")],seqtab_12SV5_tax_long_meta_clean_filtered, by="SampleSite_time")

seqtab_12SV5_tax_top<- merge(seqtab_12SV5_tax_col_sum3[, c("species", "seqtab_12SV5_tax_col_sum3")],seqtab_12SV5_tax_long_meta_clean_filtered, by="species")

seqtab_12SV5_tax_normalised <- transform(seqtab_12SV5_tax_normalised, normalised_reads = reads / seqtab_12SV5_tax_col_sum2)

ggplot(seqtab_12SV5_tax_normalised , aes(x = SampleSite_time, fill = species, y = normalised_reads)) + 
  geom_bar(stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")
##scale_fill_brewer(palette="Paired")

lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3_summary = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

seqtab_12SV5_tax_normalised <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "Days")],seqtab_12SV5_tax_normalised, by="SampleSite_time")
seqtab_12SV5_tax_normalised <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "SampleSite_code")],seqtab_12SV5_tax_normalised, by="SampleSite_time")

#PLOTTING TOP Taxon
seqtab_12SV5_tax_top %>%
  mutate(species = fct_reorder(species, desc(seqtab_12SV5_tax_col_sum3))) %>%
  ggplot(aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")+scale_fill_brewer(palette="Paired")

write.csv(seqtab_12SV5_tax_normalised,"NORM_WP2A_ASV_table_LONG_filtered.csv")

#PLOTTING PHYLA
seqtab_12SV5_tax_phyla <- data.frame(do.call('rbind', strsplit(as.character(seqtab_12SV5_tax_top$species),';',fixed=TRUE)))
seqtab_12SV5_tax_phyla<-cbind(seqtab_12SV5_tax_top,seqtab_12SV5_tax_phyla)

seqtab_12SV5_tax_phyla<-seqtab_12SV5_tax_phyla %>% 
  group_by(SampleSite_time,X2) %>% 
  summarize(reads=sum(reads),seqtab_12SV5_tax_col_sum3=sum(seqtab_12SV5_tax_col_sum3)) %>%
  rename(SampleSite_time=SampleSite_time,X2=X2, reads= reads,seqtab_12SV5_tax_col_sum3=seqtab_12SV5_tax_col_sum3)

seqtab_12SV5_tax_phyla <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "SampleSite_code")],seqtab_12SV5_tax_phyla, by="SampleSite_time")
unique(seqtab_12SV5_tax_phyla[c("X2")])
seqtab_12SV5_tax_phyla$X2<-factor(seqtab_12SV5_tax_phyla$X2,levels= c('Arthropoda','Rotifera','Mollusca','Chordata',
                                                                      'Annelida','Nematoda','Gastrotricha', 'Tardigrada','Cnidaria','Porifera','Echinodermata','Nemertea',
                                                                      'Xenacoelomorpha',
                                                                      'Nematomorpha', 'Acanthocephala <thorny-headed worm>','Chaetognatha','Hemichordata',
                                                                      'Ctenophora <comb jellyfish phylum>','Bryozoa','Platyhelminthes', 'Micrognathozoa'),ordered=TRUE)                                                                                         


ggplot(data=seqtab_12SV5_tax_phyla,aes(x = SampleSite_time, fill = X2, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  theme(axis.text.x = element_text(angle = 90))+facet_wrap(~ SampleSite_code,scales="free",ncol = 15)+scale_fill_manual(values = c("#A6CEE3",#Arthropods
                                                                                                                                   "#1F78B4", #Rotifera
                                                                                                                                   "#B2DF8A",#Molluscs
                                                                                                                                   "#33A02C",#Chordates
                                                                                                                                   "#FB9A99",#annelids
                                                                                                                                   "#E31A1C",#nematodes
                                                                                                                                   "#FDBF6F",#gastrotricha
                                                                                                                                   "#FF7F00",#tardigrada
                                                                                                                                   "#CAB2D6",#Cnideria
                                                                                                                                   "#6A3D9A",#porifera
                                                                                                                                   "#FFFF99",#echinodermata
                                                                                                                                   "#B15928",#nemertea,
                                                                                                                                   "white",
                                                                                                                                   "white","white","white","white","white","white","white","white"))
#PLOTTING TOP PHYLA
seqtab_12SV5_tax_top_phyla <- data.frame(do.call('rbind', strsplit(as.character(seqtab_12SV5_tax_top$species),';',fixed=TRUE)))
seqtab_12SV5_tax_top_phyla<-cbind(seqtab_12SV5_tax_top,seqtab_12SV5_tax_top_phyla)

seqtab_12SV5_tax_top_phyla<-seqtab_12SV5_tax_top_phyla %>% 
  group_by(SampleSite_time,X2) %>% 
  summarize(reads=sum(reads),seqtab_12SV5_tax_col_sum3=sum(seqtab_12SV5_tax_col_sum3)) %>%
  rename(SampleSite_time=SampleSite_time,X2=X2, reads= reads,seqtab_12SV5_tax_col_sum3=seqtab_12SV5_tax_col_sum3)

seqtab_12SV5_tax_top_phyla <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "SampleSite_code")],seqtab_12SV5_tax_top_phyla, by="SampleSite_time")

#seqtab_12SV5_tax_top_phyla %>%
#mutate(X2 = fct_reorder(X2, desc(seqtab_12SV5_tax_col_sum3))) %>%
#  ggplot(aes(x = SampleSite_time, fill = X2, y = reads)) + 
#  geom_bar(position="fill",stat = "identity", colour = "black")+
#  theme(axis.title.y=element_blank(),
#        axis.text.y=element_blank(),
#        axis.ticks.y=element_blank())+
# theme(axis.text.x = element_text(angle = 90))+scale_fill_brewer(palette="Paired")+facet_wrap(~ SampleSite_code,scales="free",ncol = 15)

seqtab_12SV5_tax_top_phyla<-seqtab_12SV5_tax_top_phyla %>% 
  group_by(X2) %>% 
  summarize(reads=sum(reads)) %>%
  rename(X2=X2, reads= reads)

#CONVERTING NORMALISED LONG FORMAT INTO WIDE FORMAT
seqtab_12SV5_tax_normalised_wide<-seqtab_12SV5_tax_normalised
seqtab_12SV5_tax_normalised_wide$read_filter <- NULL 
seqtab_12SV5_tax_normalised_wide$Days <- NULL 
seqtab_12SV5_tax_normalised_wide$Date_sampled <- NULL 
seqtab_12SV5_tax_normalised_wide$SampleSite_code <- NULL 
seqtab_12SV5_tax_normalised_wide$WP <- NULL 
seqtab_12SV5_tax_normalised_wide$Seq_length <- NULL 
seqtab_12SV5_tax_normalised_wide$reads <- NULL 
seqtab_12SV5_tax_normalised_wide$seqtabNoC_WP2_12SV5_tax_col_sum2 <- NULL 

seqtab_12SV5_tax_normalised_wide<-seqtab_12SV5_tax_normalised_wide %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, normalised_reads= normalised_reads,)

seqtab_12SV5_tax_normalised_wide <- spread(seqtab_12SV5_tax_normalised_wide,SampleSite_time,normalised_reads)
seqtab_12SV5_tax_normalised_wide[is.na(seqtab_12SV5_tax_normalised_wide)] <- 0
seqtab_12SV5_tax_normalised_wide <- data.frame(t(seqtab_12SV5_tax_normalised_wide))

names(seqtab_12SV5_tax_normalised_wide) <- as.matrix(seqtab_12SV5_tax_normalised_wide[1, ])
seqtab_12SV5_tax_normalised_wide <- seqtab_12SV5_tax_normalised_wide[-1, ]
seqtab_12SV5_tax_normalised_wide[] <- lapply(seqtab_12SV5_tax_normalised_wide, function(x) type.convert(as.character(x)))

write.csv(seqtab_12SV5_tax_normalised_wide,"NORM_WP2A_ASV_table_wide_filtered.csv")


# =============================================================================
# PERMANOVA (ALL TAXA, THEN PER-PHYLUM SUBSETS)
#
# vegan::adonis2() fits a permutation-based multivariate ANOVA on the
# Bray-Curtis dissimilarity matrix of the normalised community matrix.
# by='margin' fits each predictor after all others (marginal / type III SS),
# equivalent to asking 'does this variable explain additional variance
# beyond all other predictors?'. perm=999 gives p-value resolution of 0.001.
#
# Workflow for each taxon group (all taxa, then Arthropoda, Rotifera,
# Nematoda, Annelida, Mollusca):
#  1. Load the normalised wide table and remove outliers identified
#     interactively from the NMDS plot (IDs hard-coded after visual check).
#  2. For phylum subsets, columns are selected by grepl() on the taxonomy
#     string; all-zero columns and all-zero rows are then removed to avoid
#     rank-deficient dissimilarity matrices.
#  3. Join 16 environmental covariates from metadata (physical, chemical,
#     hydrological). Commented land-cover variables were tested but dropped.
#  4. Split the metadata (first 16 columns) from the community matrix before
#     passing both to adonis2().
#  5. Print the PERMANOVA result table (pseudo-F, R2, p-value per predictor).
# =============================================================================
###### PERMANOVA ON ALL TAXON ######
setwd("WP2_analysis/12SV5_SILVAngs")
NORM_WP2A_ASV_table_wide_filtered <- read_csv("NORM_WP2A_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2A_ASV_table_wide_filtered
#REMOVE OUTLIERS IDENTIFIED BY THE NMDS PLOT
adonis_data<-adonis_data[adonis_data$X1 != "E04_T07", ]  
adonis_data<-adonis_data[adonis_data$X1 != "E01_T19", ]  
adonis_data<-adonis_data[adonis_data$X1 != "E02_T01", ]  
adonis_data<-adonis_data[adonis_data$X1 != "E14_T01", ]  
adonis_data<-adonis_data[adonis_data$X1 != "E01_T09", ]  
adonis_data<-adonis_data[adonis_data$X1 != "E01_T18", ]  

adonis_data_samples<-adonis_data$X1
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(ID = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = X1)

lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","river_section","season","Days","Dist_from_lake","pH",
                                               "conductivity_uS_cm","gran_alk_ueq_L","daily_flow_by_catchment","monthly_flow_by_catchment","seasonal_flow_by_catchment",
                                               "week_before_flow_by_catchment","river_width_m","river_gradient_percent","rainfall_monthly_average",
                                               "temp_monthly_average"
                                               #"Arable_and_horticulture_ha","Broadleaf_woodland_ha",	"Coniferous_woodland_ha",
                                               #"Acid_grassland_ha",	"Neutral_grassland_ha",	"Improved_grassland_ha",	"Heather_grassland_ha",	"Heather_ha",
                                               #"Bog_ha",	"Freshwater_ha",	"Inland_rock_ha",	"Supralittoral_sediment_ha",	"Suburban_ha",	"Urban_ha"
)],adonis_data, by="SampleSite_time")

adonis_data<-na.omit(adonis_data)

adonis_data_meta<- adonis_data[c(1:16)]
adonis_data<-adonis_data[,-c(1:16)]

adonis_data_meta$Days<-as.numeric(adonis_data_meta$Days)

adon.results_WP2<-adonis2(adonis_data ~ adonis_data_meta$Dist_from_lake+ 
                            adonis_data_meta$season+
                            adonis_data_meta$daily_flow_by_catchment+
                            adonis_data_meta$monthly_flow_by_catchment+
                            adonis_data_meta$seasonal_flow_by_catchment+
                            adonis_data_meta$week_before_flow_by_catchment+
                            adonis_data_meta$river_width_m+
                            adonis_data_meta$river_gradient_percent+
                            adonis_data_meta$rainfall_monthly_average+
                            adonis_data_meta$temp_monthly_average+
                            adonis_data_meta$gran_alk_ueq_L+
                            adonis_data_meta$conductivity_uS_cm+
                            adonis_data_meta$pH,
                          method="bray",perm=999,by="margin")
print(adon.results_WP2)

#PLOTTING BETA DIVERSITY NMDS
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
data.scores_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","Blank","River","SampleSite_code","Days","Dist_from_lake","season")],data.scores, by="SampleSite_time")


p<-
  ggplot(data.scores_meta_data, aes(x = NMDS1, y = NMDS2,colour = SampleSite_time)) + 
  geom_point(size = 3.5, aes(colour = SampleSite_time))+
  theme_bw()+ theme(legend.position = "none") +
  theme(text = element_text(size=30))   
ggplotly(p)

gg <- merge(data.scores_meta_data,aggregate(cbind(NMDS1.x=NMDS1,NMDS2.y=NMDS2)~season,data.scores_meta_data,mean),by="season")
ggplot(gg, aes(x = NMDS1, y = NMDS2,colour = Dist_from_lake)) + scale_shape_manual(values=c(15,16,17,8,18))+
  geom_point(size = 4, aes(shape=season,colour = Dist_from_lake))+
  theme_bw()+  theme(axis.title.x=element_blank())+  theme(axis.title.y=element_blank())+ 
  geom_segment(aes(x=NMDS1.x, y=NMDS2.y, xend=NMDS1, yend=NMDS2))+
  geom_point(aes(x=NMDS1.x,y=NMDS2.y),size=3,colour="red")+
  theme(text = element_text(size=30))


#Canonical Correspondence Analysis (CCA)
ccamodel_adonis_data_meta<-adonis_data_meta
ccamodel_adonis_data_meta$SampleSite_time<-NULL

ccamodel<-cca(adonis_data ~.,ccamodel_adonis_data_meta)

ccamodel_adonis_data_meta$season<-as.factor(ccamodel_adonis_data_meta$season)
ccamodel_adonis_data_meta$SampleSite_code<-as.factor(ccamodel_adonis_data_meta$SampleSite_code)
ccamodel_adonis_data_meta$river_section<-as.factor(ccamodel_adonis_data_meta$river_section)

plot(ccamodel,scaling="sites")
points(ccamodel, pch=21, cex=1.2, scaling="sites")
with(ccamodel_adonis_data_meta, points(ccamodel, scaling="sites", pch=21, col=1, bg=river_section))
with(ccamodel_adonis_data_meta, legend("bottomleft", levels(river_section), pch=21, pt.bg=1:12, bty="n"))



# --- Per-phylum PERMANOVA / NMDS / CCA: Arthropoda ------------------
# Community matrix subsetted to columns containing 'Arthropoda' in the
# taxonomy string via grepl(). All-zero columns and rows are removed to
# avoid rank-deficient dissimilarity matrices. The PERMANOVA, NMDS, and
# CCA workflow from Sections 11-13 is then applied identically.
# This block repeats for Rotifera, Nematoda, Annelida, and Mollusca.
# For 12S: replace 'Arthropoda' with an order name e.g. 'Cypriniformes'.
# ---------------------------------------------------------------------
###### PERMANOVA ON ARTHROPODS ######
NORM_WP2A_ASV_table_wide_filtered <- read_csv("NORM_WP2A_ASV_table_wide_filtered.csv")
adonis_data_Arthropoda<-NORM_WP2A_ASV_table_wide_filtered

#REMOVE OUTLIERS IDENTIFIED BY THE NMDS PLOT
adonis_data_Arthropoda<-adonis_data_Arthropoda[adonis_data_Arthropoda$X1 != "E12_T03", ]  
adonis_data_Arthropoda<-adonis_data_Arthropoda[adonis_data_Arthropoda$X1 != "E05_T10", ]  
adonis_data_Arthropoda<-adonis_data_Arthropoda[adonis_data_Arthropoda$X1 != "E04_T05", ]  
adonis_data_Arthropoda<-adonis_data_Arthropoda[adonis_data_Arthropoda$X1 != "E06_T10", ]  

adonis_data_Arthropoda<-adonis_data_Arthropoda %>% rename(SampleSite_time_Arthropoda = X1)
adonis_data_Arthropoda <- adonis_data_Arthropoda[,grepl("Arthropoda", colnames(adonis_data_Arthropoda))]
adonis_data_Arthropoda<-adonis_data_Arthropoda[, colSums(adonis_data_Arthropoda != 0) > 0]
adonis_data_Arthropoda<-adonis_data_Arthropoda[apply(adonis_data_Arthropoda[,-1], 1, function(x) !all(x==0)),]

adonis_data_Arthropoda_samples<-adonis_data_Arthropoda$SampleSite_time_Arthropoda
adonis_data_Arthropoda_samples<-data.frame(adonis_data_Arthropoda_samples)
adonis_data_Arthropoda_samples<-adonis_data_Arthropoda_samples %>% rename(ID = adonis_data_Arthropoda_samples)
adonis_data_Arthropoda<-adonis_data_Arthropoda %>% rename(SampleSite_time = SampleSite_time_Arthropoda)

lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data_Arthropoda<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","river_section","season","Days","Dist_from_lake","pH",
                                                          "conductivity_uS_cm","gran_alk_ueq_L","daily_flow_by_catchment","monthly_flow_by_catchment","seasonal_flow_by_catchment",
                                                          "week_before_flow_by_catchment","river_width_m","river_gradient_percent","rainfall_monthly_average",
                                                          "temp_monthly_average")],adonis_data_Arthropoda, by="SampleSite_time")
adonis_data_Arthropoda<-na.omit(adonis_data_Arthropoda)

adonis_data_Arthropoda_meta<- adonis_data_Arthropoda[c(1:16)]
adonis_data_Arthropoda<-adonis_data_Arthropoda[,-c(1:16)]

adonis_data_Arthropoda_meta$Days<-as.numeric(adonis_data_Arthropoda_meta$Days)

adon.results_WP2_Arthropoda<-adonis2(adonis_data_Arthropoda ~ 
                                       adonis_data_Arthropoda_meta$season+
                                       adonis_data_Arthropoda_meta$Days+
                                       adonis_data_Arthropoda_meta$Dist_from_lake+
                                       adonis_data_Arthropoda_meta$daily_flow_by_catchment+
                                       adonis_data_Arthropoda_meta$monthly_flow_by_catchment+
                                       adonis_data_Arthropoda_meta$week_before_flow_by_catchment+
                                       adonis_data_Arthropoda_meta$river_gradient_percent+
                                       adonis_data_Arthropoda_meta$rainfall_monthly_average+
                                       adonis_data_Arthropoda_meta$temp_monthly_average+
                                       adonis_data_Arthropoda_meta$conductivity_uS_cm+
                                       adonis_data_Arthropoda_meta$pH,
                                     method="bray",perm=999,by="margin")
print(adon.results_WP2_Arthropoda)

#PLOTTING BETA DIVERSITY NMDS
com_Arthropoda <- adonis_data_Arthropoda[,col(adonis_data_Arthropoda)]
m_com_Arthropoda <- as.matrix(com_Arthropoda)
set.seed(123)
nmds_Arthropoda = metaMDS(m_com_Arthropoda, distance = "bray", try = 1000, trymax = 1000, k = 3)

#extract NMDS scores (x and y coordinates)
data.scores_Arthropoda = as.data.frame(scores(nmds_Arthropoda))
#add columns to data frame 
data.scores_Arthropoda$SampleSite_time = adonis_data_Arthropoda_meta[,1]

lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]
data.scores_Arthropoda_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","Blank","SampleSite_code","Days","Dist_from_lake")],data.scores_Arthropoda, by="SampleSite_time")

p<-
  ggplot(data.scores_Arthropoda_meta_data, aes(x = NMDS1, y = NMDS2,colour = SampleSite_time)) + 
  geom_point(size = 3.5, aes(colour = SampleSite_time))+
  theme_bw()+
  theme(text = element_text(size=30))   + theme(legend.position = "none") 
ggplotly(p)


#Canonical Correspondence Analysis (CCA) ARTHROPODA
ccamodel_adonis_data_Arthropoda_meta<-adonis_data_Arthropoda_meta
ccamodel_adonis_data_Arthropoda_meta$SampleSite_time<-NULL

ccamodel_Arthropoda<-cca(adonis_data_Arthropoda ~.,ccamodel_adonis_data_Arthropoda_meta)

ccamodel_adonis_data_Arthropoda_meta$season<-as.factor(ccamodel_adonis_data_Arthropoda_meta$season)
ccamodel_adonis_data_Arthropoda_meta$SampleSite_code<-as.factor(ccamodel_adonis_data_Arthropoda_meta$SampleSite_code)
ccamodel_adonis_data_Arthropoda_meta$river_section<-as.factor(ccamodel_adonis_data_Arthropoda_meta$river_section)

plot(ccamodel_Arthropoda,scaling="sites")
points(ccamodel_Arthropoda, pch=21, cex=1.2, scaling="sites")
with(ccamodel_adonis_data_Arthropoda_meta, points(ccamodel_Arthropoda, scaling="sites", pch=21, col=1, bg=river_section))
with(ccamodel_adonis_data_Arthropoda_meta, legend("bottomleft", levels(river_section), pch=21, pt.bg=1:12, bty="n"))

# =============================================================================
# SECTION 14: ALPHA DIVERSITY METRICS & GAM MODELLING
#
# Alpha (within-sample) diversity is calculated for all taxa combined,
# then separately for each major phylum, using the normalised wide matrix
# as input to vegan::diversity().
#
# Three metrics are computed:
#  - Shannon index (H'): H' = -sum(p_i * ln(p_i)). Sensitive to both
#    species richness and evenness. The most commonly reported metric.
#  - Simpson index (D): D = 1 - sum(p_i^2). Weighted toward dominant
#    taxa; ranges 0 (no diversity) to 1 (maximum diversity).
#  - Species richness: count of ASVs with reads > 0 per sample.
#
# For each phylum (Arthropoda, Rotifera, Nematoda, Mollusca, Annelida,
# Porifera, Chordata):
#  1. Subset the normalised wide matrix to phylum-specific columns via grepl().
#  2. Compute Shannon diversity on the submatrix.
#  3. Attach metadata (Dist_from_lake, Days, season, env. covariates).
#  4. Add a 'phyla' label column so rows can be identified after rbind().
#  5. Plot a boxplot: Shannon by site, coloured by distance from lake.
#     (outlier.shape=NA + geom_jitter avoids double-plotting outliers.)
#  6. Exploratory single-predictor GAM scatter plots (one per covariate)
#     to assess which variables warrant inclusion in the full model.
#  7. Fit a full multivariate GAM (mgcv::gam) with:
#     - season: parametric factor, shifts the intercept per season.
#     - s(Dist_from_lake): non-linear spatial diversity gradient.
#     - s(Dist_from_lake, by=season): season-specific spatial gradients
#       (varying-coefficient model -- tests whether the shape of the
#       spatial gradient differs between seasons).
#     - s(Days): temporal trend within the sampling window.
#     - Smooth terms for all remaining env. covariates.
#  8. Extract the smooth-term summary table (s.table), drop edf/Ref.df,
#     and format each term as 'F = X.XX, p = Y.YY' strings.
#     'p = 0' is replaced with 'p < 0.01' to avoid false precision.
#
# All five phylum GAM summaries are then combined into a single
# publication-ready table via reduce(left_join, by='Factor'), transposed
# so rows = phyla and columns = smooth terms, and written to CSV.
#
# Multi-panel figures are assembled with ggpubr::ggarrange() for both
# the boxplots (p_*) and the GAM Shannon vs. Dist_from_lake plots
# (gam_*_plot) across all five phyla.
#
# Per-phylum alpha data are stacked with rbind() and a combined-phyla
# exploratory GAM model is fit with phyla and season-specific spatial
# and temporal smooth interactions.
# =============================================================================
##########   ALPHA DIVERSITY METRICS ################

shannon_index<-diversity(seqtab_12SV5_tax_normalised_wide, index = "shannon")
shannon_index<-data.frame(shannon_index)
ID <- rownames(shannon_index)
shannon_index <- cbind(shannon_index,ID)

simpson_index<-diversity(seqtab_12SV5_tax_normalised_wide, index = "simpson")
simpson_index<-data.frame(simpson_index)
ID <- rownames(simpson_index)
simpson_index <- cbind(simpson_index,ID)

species_count <- rowSums(seqtab_12SV5_tax_normalised_wide >0)
species_count<-data.frame(species_count)
species_count$ID <- rownames(species_count)

alpha_12SV5_WP2A<-merge(shannon_index,simpson_index,by="ID")
alpha_12SV5_WP2A<-merge(alpha_12SV5_WP2A,species_count,by="ID")

write.csv(alpha_12SV5_WP2A,"alpha_12SV5_WP2A.csv")

lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

alpha_12SV5_WP2A_meta_data <- alpha_12SV5_WP2A %>% rename(SampleSite_time= ID)
alpha_12SV5_WP2A_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","SampleSite_code","Date_sampled","Dist_from_lake","Days","season")],alpha_12SV5_WP2A_meta_data, by="SampleSite_time")

ggplot(alpha_12SV5_WP2A_meta_data, aes(x=Dist_from_lake, y=shannon_index, color=season)) +
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+
  theme_bw()


#ALPHA DIVERSITY Arthropoda
seqtab_12SV5_tax_normalised_wide_Arthropoda <- seqtab_12SV5_tax_normalised_wide[ , grepl( "Arthropoda" , names( seqtab_12SV5_tax_normalised_wide ) ) ]
shannon_index_Arthropoda<-diversity(seqtab_12SV5_tax_normalised_wide_Arthropoda, index = "shannon")
shannon_index_Arthropoda<-data.frame(shannon_index_Arthropoda)
ID_Arthropoda <- rownames(shannon_index_Arthropoda)
shannon_index_Arthropoda <- cbind(shannon_index_Arthropoda,ID_Arthropoda)

lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

shannon_index_Arthropoda_meta_data <- shannon_index_Arthropoda %>% rename(SampleSite_time= ID_Arthropoda)
shannon_index_Arthropoda_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","SampleSite_code","Date_sampled","Dist_from_lake","Days")],shannon_index_Arthropoda_meta_data, by="SampleSite_time")
shannon_index_Arthropoda_meta_data$phyla <- paste0("Arthropoda", shannon_index_Arthropoda_meta_data$phyla)
shannon_index_Arthropoda_meta_data<-shannon_index_Arthropoda_meta_data %>% rename(shannon_index = shannon_index_Arthropoda)

p_Arthropoda<-ggplot(shannon_index_Arthropoda_meta_data, aes(x=SampleSite_code, y=shannon_index,colour=Dist_from_lake)) +
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter()+ 
  theme_bw()

shannon_index_Arthropoda_meta_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time",
                                                                      "season",
                                                                      "temp_monthly_average",
                                                                      "conductivity_uS_cm",
                                                                      "gran_alk_ueq_L","daily_flow_by_catchment",
                                                                      "monthly_flow_by_catchment",
                                                                      "seasonal_flow_by_catchment",
                                                                      "week_before_flow_by_catchment",
                                                                      "river_width_m",
                                                                      "river_gradient_percent",
                                                                      "rainfall_monthly_average",
                                                                      "pH")],shannon_index_Arthropoda_meta_data, by="SampleSite_time")

shannon_index_Arthropoda_meta_data$Days<-as.numeric(shannon_index_Arthropoda_meta_data$Days)

shannon_index_Arthropoda_meta_data$season<-as.factor(shannon_index_Arthropoda_meta_data$season)

gam_Arthropoda_plot<-ggplot(shannon_index_Arthropoda_meta_data, aes(Dist_from_lake, shannon_index,colour=season)) + 
  geom_point(size=3) + 
  geom_smooth(method = "gam", formula = y ~ s(x),size=2)+theme_bw()+
  theme(legend.position="none")+
  theme(text = element_text(size=30))+ theme(axis.title.y = element_blank())

ggplot(shannon_index_Arthropoda_meta_data, aes(daily_flow_by_catchment, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+theme_bw()+
  theme(legend.position="none")

ggplot(shannon_index_Arthropoda_meta_data, aes(week_before_flow_by_catchment, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+theme_bw()+
  theme(legend.position="none")

ggplot(shannon_index_Arthropoda_meta_data, aes(Days, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+theme_bw()+
  theme(legend.position="none")

ggplot(shannon_index_Arthropoda_meta_data, aes(pH, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+theme_bw()+
  theme(legend.position="none")

mod_gam_Arthropoda <- mgcv::gam(shannon_index ~ 
                                  season+
                                  s(Dist_from_lake)+
                                  s(Days)+
                                  s(Dist_from_lake, by=season)+
                                  s(pH)+
                                  s(daily_flow_by_catchment)+
                                  s(monthly_flow_by_catchment)+
                                  s(week_before_flow_by_catchment)+
                                  s(river_gradient_percent)+
                                  s(rainfall_monthly_average)+
                                  s(temp_monthly_average)+
                                  s(conductivity_uS_cm), 
                                data=shannon_index_Arthropoda_meta_data)
summary(mod_gam_Arthropoda)
anova(mod_gam_Arthropoda)

mod_gam_Arthropoda_summary<-summary(mod_gam_Arthropoda)
mod_gam_Arthropoda_summary<-mod_gam_Arthropoda_summary[["s.table"]]
mod_gam_Arthropoda_summary<-data.frame(mod_gam_Arthropoda_summary)
mod_gam_Arthropoda_summary$edf<-NULL
mod_gam_Arthropoda_summary$Ref.df<-NULL
mod_gam_Arthropoda_summary<-round(mod_gam_Arthropoda_summary, digits = 2)
mod_gam_Arthropoda_summary$F<-paste("F = ",mod_gam_Arthropoda_summary$F,  sep = "")
mod_gam_Arthropoda_summary$p.value<-paste("p = ",mod_gam_Arthropoda_summary$p.value,  sep = "")

mod_gam_Arthropoda_summary$p.value<-sub('^p = 0$', 'p < 0.01', mod_gam_Arthropoda_summary$p.value)

mod_gam_Arthropoda_summary$Arthropoda<-paste(mod_gam_Arthropoda_summary$F,mod_gam_Arthropoda_summary$p.value,
                                             sep = ", ")
mod_gam_Arthropoda_summary$F<-NULL
mod_gam_Arthropoda_summary$p.value<-NULL


#ALPHA DIVERSITY Porifera
seqtab_12SV5_tax_normalised_wide_Porifera <- seqtab_12SV5_tax_normalised_wide[ , grepl( "Porifera" , names( seqtab_12SV5_tax_normalised_wide ) ) ]
shannon_index_Porifera<-diversity(seqtab_12SV5_tax_normalised_wide_Porifera, index = "shannon")
shannon_index_Porifera<-data.frame(shannon_index_Porifera)
ID_Porifera <- rownames(shannon_index_Porifera)
shannon_index_Porifera <- cbind(shannon_index_Porifera,ID_Porifera)

lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

shannon_index_Porifera_meta_data <- shannon_index_Porifera %>% rename(SampleSite_time= ID_Porifera)
shannon_index_Porifera_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","SampleSite_code","Date_sampled","Dist_from_lake","Days")],shannon_index_Porifera_meta_data, by="SampleSite_time")
shannon_index_Porifera_meta_data$phyla <- paste0("Porifera", shannon_index_Porifera_meta_data$phyla)
shannon_index_Porifera_meta_data<-shannon_index_Porifera_meta_data %>% rename(shannon_index = shannon_index_Porifera)

p_Porifera<-ggplot(shannon_index_Porifera_meta_data, aes(x=SampleSite_code, y=shannon_index,colour=Dist_from_lake)) +
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter()+ 
  theme_bw()

#combining all of the GAM outputs from each of the phyla

mod_gam_Arthropoda_summary <- tibble::rownames_to_column(mod_gam_Arthropoda_summary, "Factor")
mod_gam_Rotifera_summary <- tibble::rownames_to_column(mod_gam_Rotifera_summary, "Factor")
mod_gam_Nematoda_summary <- tibble::rownames_to_column(mod_gam_Nematoda_summary, "Factor")
mod_gam_Annelida_summary <- tibble::rownames_to_column(mod_gam_Annelida_summary, "Factor")
mod_gam_Mollusca_summary <- tibble::rownames_to_column(mod_gam_Mollusca_summary, "Factor")

mod_gam_combined_summary<-list(mod_gam_Arthropoda_summary ,
                               mod_gam_Rotifera_summary ,
                               mod_gam_Nematoda_summary,
                               mod_gam_Mollusca_summary, 
                               mod_gam_Annelida_summary)%>% reduce(left_join, by = "Factor")

mod_gam_combined_summary[is.na(mod_gam_combined_summary)] <- "-"


mod_gam_combined_summary<-t(mod_gam_combined_summary)
mod_gam_combined_summary<-data.frame(mod_gam_combined_summary)

write.csv(mod_gam_combined_summary,"mod_gam_combined_summary.csv")

library(ggpubr)
ggarrange(p_Annelida,p_Arthropoda,p_Mollusca,p_Nematoda,p_Rotifera,ncol = 2,nrow=3)
ggarrange(gam_Arthropoda_plot,gam_Rotifera_plot,gam_Nematoda_plot,gam_Annelida_plot,
          gam_Mollusca_plot,ncol = 2,nrow=3)


#combine alpha diversity from different phyla
shannon_index_phlya_split_data<-rbind(shannon_index_Nematoda_meta_data,
                                      shannon_index_Annelida_meta_data,
                                      shannon_index_Rotifera_meta_data,
                                      shannon_index_Arthropoda_meta_data)

######## WP2A alpha diversity modeling ##########
#library(mgcv)

shannon_index_phlya_split_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time",
                                                                  "season",
                                                                  "temp_monthly_average",
                                                                  "conductivity_uS_cm",
                                                                  "gran_alk_ueq_L","daily_flow_by_catchment",
                                                                  "monthly_flow_by_catchment",
                                                                  "seasonal_flow_by_catchment",
                                                                  "week_before_flow_by_catchment",
                                                                  "river_width_m",
                                                                  "river_gradient_percent",
                                                                  "rainfall_monthly_average",
                                                                  "pH")],shannon_index_phlya_split_data, by="SampleSite_time")

shannon_index_phlya_split_data$Days<-as.numeric(shannon_index_phlya_split_data$Days)
shannon_index_phlya_split_data<-shannon_index_phlya_split_data[shannon_index_phlya_split_data$shannon_index != 0, ] 

#GAMs
ggplot(shannon_index_phlya_split_data, aes(Days, shannon_index,colour=phyla)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+
  facet_wrap(~SampleSite_code)

ggplot(shannon_index_phlya_split_data, aes(temp_monthly_average, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+
  facet_wrap(~SampleSite_code)

ggplot(shannon_index_phlya_split_data, aes(daily_flow_by_catchment, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+
  facet_wrap(~SampleSite_code)

ggplot(shannon_index_phlya_split_data, aes(conductivity_uS_cm, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+
  facet_wrap(~SampleSite_code)

ggplot(shannon_index_phlya_split_data, aes(pH, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+
  facet_wrap(~SampleSite_code)

ggplot(shannon_index_phlya_split_data, aes(gran_alk_ueq_L, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+
  facet_wrap(~SampleSite_code)

ggplot(shannon_index_phlya_split_data, aes(daily_flow_by_catchment, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+
  facet_wrap(~SampleSite_code)

ggplot(shannon_index_phlya_split_data, aes(rainfall_monthly_average, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+
  facet_wrap(~SampleSite_code)

ggplot(shannon_index_phlya_split_data, aes(monthly_flow_by_catchment, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+
  facet_wrap(~SampleSite_code)

ggplot(shannon_index_phlya_split_data, aes(river_width_m, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))

ggplot(shannon_index_phlya_split_data, aes(river_gradient_percent, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))

ggplot(shannon_index_phlya_split_data, aes(season, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+
  facet_wrap(~SampleSite_code)

ggplot(shannon_index_phlya_split_data, aes(Dist_from_lake, shannon_index,colour=season)) + 
  geom_point() + 
  #geom_smooth()+
  geom_smooth(method = "gam", formula = y ~ s(x))+
  facet_wrap(~phyla)+theme_bw()

shannon_index_phlya_split_data$phyla<-as.factor(shannon_index_phlya_split_data$phyla)
shannon_index_phlya_split_data$season<-as.factor(shannon_index_phlya_split_data$season)

mod_gam <- gam(shannon_index ~ phyla+
                 s(Days)+
                 season+
                 Dist_from_lake+
                 s(Dist_from_lake, by=season)+
                 s(Dist_from_lake, by=phyla)+
                 s(Days, by=phyla)+
                 s(Dist_from_lake, by=pH)+
                 daily_flow_by_catchment+
                 monthly_flow_by_catchment+
                 seasonal_flow_by_catchment+
                 week_before_flow_by_catchment+
                 river_width_m+
                 river_gradient_percent+
                 rainfall_monthly_average+
                 temp_monthly_average+
                 s(gran_alk_ueq_L)+
                 s(conductivity_uS_cm), 
               data=shannon_index_phlya_split_data)
summary(mod_gam)

plot(mod_gam, shade = TRUE, pages = 1, scale = 0)

plot_smooths(
  model = mod_gam2,
  series = Dist_from_lake,
  comparison = Days,
  facet_terms = phyla,
  split = list(Daysphyla = c("Days", "phyla"))
)


# =============================================================================
# BETA-DIVERSITY PARTITIONING (betapart / Sorensen)
#
# betapart::beta.pair() decomposes pairwise Sorensen beta-diversity (B_sor)
# into two biologically interpretable additive components:
#
#   [[1]] B_sim (Simpson / turnover): the fraction of beta-diversity driven
#         by species replacement between sites -- i.e. communities differ
#         because different species are present, not just more or fewer.
#
#   [[2]] B_sne (nestedness): the fraction of beta-diversity driven purely
#         by richness differences -- i.e. species-poor sites are subsets of
#         richer sites. High B_sne indicates a nested community structure.
#
#   [[3]] B_sor (total Sorensen): sum of turnover + nestedness.
#
# Workflow (repeated for each river: WP2A, Conwy, Glatt, Gwash, Towy, Skan):
#  1. Load the normalised wide table and join metadata including
#     Dist_from_lake and Days for each site.
#  2. Binarise the community matrix (presence/absence) using mutate_if;
#     Sorensen decomposition requires binary data.
#  3. Drop all-zero columns (ASVs absent from this river subset).
#  4. Run beta.pair() to get the three pairwise distance matrices.
#  5. Melt each matrix to long format (SiteID_1, SiteID_2, Beta_*).
#  6. Attach site 1 and site 2 metadata (Dist_from_lake, Days).
#  7. Compute Dist_from_lake_DIFF = site2 - site1 (positive = downstream).
#  8. Plot each beta component vs. Dist_from_lake_DIFF with geom_smooth();
#     a positive relationship indicates distance-decay of similarity.
#     scale_x_continuous(limits=c(0,NA)) restricts to the downstream
#     direction only.
#
# After all rivers are processed, results are combined with rbind(),
# the River label is joined from metadata, and only downstream-directed
# pairs (Dist_from_lake_DIFF >= 0) are retained to avoid reciprocal
# duplicates. Multi-river distance-decay plots colour points by River
# using a colourblind-friendly palette (Wong 2011).
#
# NOTE for 12S / Surama Creek: this analysis directly tests whether fish
# community turnover increases with downstream distance from Surama Lake.
# High B_sim + low B_sne = communities differ by species replacement
# (typical of upstream-downstream gradients in Amazonian rivers).
# High B_sne + low B_sim = downstream communities are nested subsets of
# upstream ones (typical of filtered / impoverished assemblages).
# =============================================================================
######## Turnover and nestedness ########
setwd("WP2_analysis/12SV5_SILVAngs")
NORM_WP2A_ASV_table_wide_filtered <- read_csv("NORM_WP2A_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2A_ASV_table_wide_filtered

SampleSite_time<-NORM_WP2A_ASV_table_wide_filtered$X1
SampleSite_time<-data.frame(SampleSite_time)
SampleSite_time <- tibble::rownames_to_column(SampleSite_time, "SiteID_1")

adonis_data_samples<-adonis_data$X1
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(ID = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = X1)

lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","river_section","season","Days","Dist_from_lake","pH",
                                               "conductivity_uS_cm","gran_alk_ueq_L","daily_flow_by_catchment","monthly_flow_by_catchment","seasonal_flow_by_catchment",
                                               "week_before_flow_by_catchment","river_width_m","river_gradient_percent","rainfall_monthly_average",
                                               "temp_monthly_average"
)],adonis_data, by="SampleSite_time")

SampleSite_time<-merge(lofresh_metadata_WP2_3[, c("SampleSite_time","Dist_from_lake","Days")],SampleSite_time, by="SampleSite_time")

adonis_data_meta<- adonis_data[c(1:16)]
adonis_data<-adonis_data[,-c(1:16)]

adonis_data<-na.omit(adonis_data)

adonis_data_meta$Days<-as.numeric(adonis_data_meta$Days)
#out <- nestedtemp(adonis_data)
#plot(out)
#plot(out, kind="incid")

library(betapart)
#library(qgraph)
adonis_data<-adonis_data %>% mutate_if(is.numeric, ~1 * (. != 0))
adonis_data<-adonis_data[, colSums(adonis_data != 0) > 0]
sorensen.beta<-beta.pair(adonis_data, index.family="sorensen")

nestedness <- as.matrix(sorensen.beta[[3]])
nestedness <- melt(nestedness)
names(nestedness) <- c("SiteID_1", "SiteID_2", "Beta_Distance")

###-----------------------
## get nestedness matrix (this is [[2]]) into a data frame
###-----------------------
m.beta.sorn <- as.matrix(sorensen.beta[[2]])
m.beta.sorn <- melt(m.beta.sorn)
names(m.beta.sorn) <- c("SiteID_1", "SiteID_2", "Beta_nest")
beta.measures.sor <- merge(nestedness, m.beta.sorn, by = c("SiteID_1","SiteID_2")) 

###-----------------------
## get turnover matrix (here simpson, this is [[1]]) into a data frame
###-----------------------
m.beta.sort <- as.matrix(sorensen.beta[[1]])
m.beta.sort  <- melt(m.beta.sort )
names(m.beta.sort ) <- c("SiteID_1", "SiteID_2", "Beta_turn")
beta.measures.sor <- merge(beta.measures.sor, m.beta.sort , by = c("SiteID_1","SiteID_2")) 

beta.measures.sor<-merge(SampleSite_time,beta.measures.sor,by="SiteID_1")
beta.measures.sor<-beta.measures.sor %>% rename(SampleSite_time_1 = SampleSite_time)
beta.measures.sor<-beta.measures.sor %>% rename(Dist_from_lake_1 = Dist_from_lake)
beta.measures.sor<-beta.measures.sor %>% rename(Days_1 = Days)

SampleSite_time<-SampleSite_time %>% rename(SiteID_2 = SiteID_1)
beta.measures.sor<-merge(SampleSite_time,beta.measures.sor,by="SiteID_2")
beta.measures.sor<-beta.measures.sor %>% rename(SampleSite_time_2 = SampleSite_time)
beta.measures.sor<-beta.measures.sor %>% rename(Dist_from_lake_2 = Dist_from_lake)
beta.measures.sor<-beta.measures.sor %>% rename(Days_2 = Days)

beta.measures.sor$Dist_from_lake_DIFF<-(beta.measures.sor$Dist_from_lake_2-beta.measures.sor$Dist_from_lake_1)
beta.measures.sor$Days_2<-as.numeric(beta.measures.sor$Days_2)
beta.measures.sor$Days_1<-as.numeric(beta.measures.sor$Days_1)
beta.measures.sor$Days_DIFF<-(beta.measures.sor$Days_2-beta.measures.sor$Days_1)
beta.measures.sor <- beta.measures.sor[beta.measures.sor$Days_DIFF >= 0, ]

ggplot(beta.measures.sor, aes(x=(Dist_from_lake_DIFF), y=Beta_turn, colour=Days_DIFF)) +
  geom_point(size=3.5, shape=23)+theme_bw()+
  geom_smooth()+scale_color_gradientn(colours = rainbow(5))+scale_x_continuous(limits = c(0, NA))+
  theme(text = element_text(size=30))

ggplot(beta.measures.sor, aes(x=(Dist_from_lake_DIFF), y=Beta_nest, colour=Days_DIFF)) +
  geom_point(size=3.5, shape=23)+theme_bw()+
  geom_smooth()+scale_color_gradientn(colours = rainbow(5))+scale_x_continuous(limits = c(0, NA))+
  theme(text = element_text(size=30))

ggplot(beta.measures.sor, aes(x=(Dist_from_lake_DIFF), y=Beta_Distance, colour=Days_DIFF)) +
  geom_point(size=3.5, shape=23)+theme_bw()+
  geom_smooth()+scale_color_gradientn(colours = rainbow(5))+scale_x_continuous(limits = c(0, NA))+
  theme(text = element_text(size=30))


#hist(beta.measures.sor$Days_DIFF_sqrd)

######## Sparse partial least squares analysis ##########
library(mixOmics)
library(ggdendro)
library(plotly)

setwd("WP2_analysis/12SV5_SILVAngs")
NORM_WP2A_ASV_table_wide_filtered <- read_csv("NORM_WP2A_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2A_ASV_table_wide_filtered
#REMOVE OUTLIERS IDENTIFIED BY THE NMDS PLOT
adonis_data<-adonis_data[adonis_data$X1 != "E04_T07", ]  
adonis_data<-adonis_data[adonis_data$X1 != "E01_T19", ]  
adonis_data<-adonis_data[adonis_data$X1 != "E02_T01", ]  
adonis_data<-adonis_data[adonis_data$X1 != "E14_T01", ]  
adonis_data<-adonis_data[adonis_data$X1 != "E01_T09", ]  
adonis_data<-adonis_data[adonis_data$X1 != "E01_T18", ]  

adonis_data_samples<-adonis_data$X1
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(ID = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = X1)

lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time",
                                               "temp_monthly_average",
                                               "conductivity_uS_cm",
                                               "daily_flow_by_catchment",
                                               "monthly_flow_by_catchment",
                                               "week_before_flow_by_catchment",
                                               "river_gradient_percent",
                                               "rainfall_monthly_average",
                                               "pH" 
)],adonis_data, by="SampleSite_time")

adonis_data<-na.omit(adonis_data)

adonis_data_meta<- adonis_data[c(1:11)]
adonis_data<-adonis_data[,-c(1:11)]

adonis_data_meta$Days<-as.numeric(adonis_data_meta$Days)
adonis_data_meta$SampleSite_time<-NULL

ncomp = 1
result.spls <- spls(adonis_data, adonis_data_meta, ncomp = ncomp,  mode = 'canonical')

design <- data.frame(samp = adonis_data_meta$sample)

cim_plot<-cim(result.spls,
              comp = 1:1,
              margins = c(15, 5))

cl<-colnames(adonis_data)
cl<-data.frame(cl)
cl<- data.frame(do.call('rbind', strsplit(as.character(cl$cl),';',fixed=TRUE)))
cl<-subset(cl,select=c(X2))

unique(cl[c("X2")])

stim.col<-c("darkblue","purple", "green4","red3","orange", "cyan","darkviolet","yellow",
            "darkorange","red","green","black","darkcyan","darkred","darkorchid","darksalmon","gray","violet","springgreen4")#19
cl$X2<-as.factor(cl$X2)
stim.col <- stim.col[as.numeric(cl$X2)]

cim_plot<-cim(result.spls,row.sideColors = stim.col,
              comp = 1:1,
              margins = c(15, 5))

spls_correlation<-cim_plot[["mat"]]
spls_correlation<-data.frame(spls_correlation)
spls_correlation <- tibble::rownames_to_column(spls_correlation, "taxa")


#Calculating 95% confidence interval on correlations
conf_95_spls_correlation<-spls_correlation
conf_95_spls_correlation$row_number<-NULL
conf_95_spls_correlation$taxa<-NULL
conf_95_spls_correlation<-data.frame(x=unlist(conf_95_spls_correlation[]))
hist(conf_95_spls_correlation$x)
boxplot(conf_95_spls_correlation$x)

quantile(conf_95_spls_correlation$x, probs = seq(0, 1, 0.01))
mean(conf_95_spls_correlation$x)

#investigating the 95% percentile of correlations 
spls_correlation_upper <- gather(spls_correlation, env_variable, correlation, Inland_rock_ha:rainfall_monthly_average, factor_key=TRUE)
spls_correlation_upper<-spls_correlation_upper %>%
  filter(correlation > 0.2695705083)
spls_correlation_upper$env_variable<-NULL
spls_correlation_upper$correlation<-NULL
spls_correlation_upper<-unique(spls_correlation_upper$taxa)
spls_correlation_upper<-data.frame(spls_correlation_upper)
spls_correlation_upper$keep<-"keep"
spls_correlation_upper<-spls_correlation_upper %>% rename(species = spls_correlation_upper)


setwd("WP2_analysis/12SV5_SILVAngs")
seqtab_12SV5_tax_long_meta_clean_filtered <- read_csv("WP2A_ASV_table_long_filtered_family.csv")

#aggregating replicates
seqtab_12SV5_tax_long_meta_clean_filtered<-seqtab_12SV5_tax_long_meta_clean_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=mean(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads)

spls_correlation_upper_abundance_plot<-seqtab_12SV5_tax_long_meta_clean_filtered
spls_correlation_upper_abundance_plot<-merge(seqtab_12SV5_tax_long_meta_clean_filtered,spls_correlation_upper,by="species")

spls_correlation_upper_abundance_plot$species<-gsub("NA;", "", spls_correlation_upper_abundance_plot$species)
spls_correlation_upper_abundance_plot$species<-gsub(";NA", "", spls_correlation_upper_abundance_plot$species)
spls_correlation_upper_abundance_plot$species<-gsub(" <Metazoa>", "", spls_correlation_upper_abundance_plot$species)
spls_correlation_upper_abundance_plot$species<-gsub(" <vertebrate>", "", spls_correlation_upper_abundance_plot$species)
spls_correlation_upper_abundance_plot$species<-gsub(" <chordata>", "", spls_correlation_upper_abundance_plot$species)

ggplot(spls_correlation_upper_abundance_plot, aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))+guides(fill=guide_legend(ncol=1))

# --- Cross-river compilation -----------------------------------------
# Combine per-river beta partitioning results, join River labels, filter
# to downstream-directed pairs (Dist_from_lake_DIFF >= 0), and produce
# multi-river distance-decay plots coloured by River.
# ---------------------------------------------------------------------
#Compiling nestedness for different rivers
beta.measures.sor_all<-rbind(beta.measures.sor_Conwy,beta.measures.sor_Glatt,beta.measures.sor_Gwash,
                             beta.measures.sor_Towy,beta.measures.sor_Skan)

lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]
beta.measures.sor_all<-beta.measures.sor_all %>% rename(SampleSite_time = SampleSite_time_1)

beta.measures.sor_all<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","River")],beta.measures.sor_all, by="SampleSite_time")

beta.measures.sor_all<-beta.measures.sor_all %>%
  group_by(SampleSite_time,SampleSite_time_2) %>%
  filter(all(Dist_from_lake_DIFF>=0))
#beta.measures.sor_all<-beta.measures.sor_all[beta.measures.sor_all$Dist_from_lake_DIFF != 0, ]


ggplot(beta.measures.sor_all, aes(x=(Dist_from_lake_DIFF), y=Beta_turn,colour=River)) +
  geom_point(size=3.5, shape=23)+theme_bw()+
  geom_smooth()+scale_x_continuous(limits = c(0, NA))+
  scale_color_manual(values=c("#56B4E9","#F0E442","#E69F00","#D55E00","#0072B2"))+
  theme(text = element_text(size=30))

ggplot(beta.measures.sor_all, aes(x=(Dist_from_lake_DIFF), y=Beta_nest,colour=River)) +
  geom_point(size=3.5, shape=23)+theme_bw()+
  geom_smooth()+scale_x_continuous(limits = c(0, NA))+
  scale_color_manual(values=c("#56B4E9","#F0E442","#E69F00","#D55E00","#0072B2"))+
  theme(text = element_text(size=30))

ggplot(beta.measures.sor_all, aes(x=(Dist_from_lake_DIFF), y=Beta_Distance,colour=River)) +
  geom_point(size=3.5, shape=23)+theme_bw()+
  geom_smooth()+scale_x_continuous(limits = c(0, NA))+
  scale_color_manual(values=c("#56B4E9","#F0E442","#E69F00","#D55E00","#0072B2"))+
  theme(text = element_text(size=30))
