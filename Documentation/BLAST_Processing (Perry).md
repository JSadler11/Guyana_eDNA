
```
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
#________________________________________________________
```
```

#Removing 12S sequences in the taxtable that are the wrong size,then BLASTing with a broad 12S database to remove
#sequences that are assigned to other things that are not fish

#________________________________________________________

taxTab_WP2_12S <- read_table2("taxTab_WP2_12S.txt")
taxTab_WP2_12S$Kingdom <-NULL
taxTab_WP2_12S$Phylum <-NULL
taxTab_WP2_12S$Class <-NULL
taxTab_WP2_12S$Order <-NULL
taxTab_WP2_12S$Family <-NULL
taxTab_WP2_12S$Genus <-NULL
taxTab_WP2_12S$Species <-NULL

#sequence length of the ASVs
taxTab_WP2_12S_length<-aggregate(Seq~Seq_length, transform(taxTab_WP2_12S, Seq_length=Seq),
                                        FUN=function(x) nchar(unique(x)))
hist(taxTab_WP2_12S_length$Seq)

#removing ASVs that are the wrong size, over 20% larger than the average amplicon size the MiFish amplicon should be (172bp)
taxTab_WP2_12S_length<-taxTab_WP2_12S_length[!(taxTab_WP2_12S_length$Seq>206.2),]

#removing ASVs that are the wrong size, under 20% smaller than the average amplicon size the MiFish amplicon should be (172bp)
taxTab_WP2_12S_length<-taxTab_WP2_12S_length[!(taxTab_WP2_12S_length$Seq<137.6),]
hist(taxTab_WP2_12S_length$Seq)

taxTab_WP2_12S_length$Seq_number <- 1:nrow(taxTab_WP2_12S_length) 
taxTab_WP2_12S_length$Seq<-"SEQ"
taxTab_WP2_12S_length$div<-">"

taxTab_WP2_12S_length$Seq_number <-paste(taxTab_WP2_12S_length$div,taxTab_WP2_12S_length$Seq,taxTab_WP2_12S_length$Seq_number)
taxTab_WP2_12S_length$Seq_number = as.character(gsub(" ", "", taxTab_WP2_12S_length$Seq_number))

taxTab_WP2_12S_length$Seq <-NULL
taxTab_WP2_12S_length$div <-NULL

write.table(taxTab_WP2_12S_length,"length_filtered_SEQ_to_sequence_lookup_table.txt",quote = FALSE)

taxTab_WP2_12S_length$Seq <-paste(taxTab_WP2_12S_length$Seq_number,taxTab_WP2_12S_length$Seq_length)

taxTab_WP2_12S_length$Seq <- as.character(gsub(" ", "\n", taxTab_WP2_12S_length$Seq))

length_adjust_Seq_number<-taxTab_WP2_12S_length$Seq_number
length_adjust_Seq_number<-data.frame(length_adjust_Seq_number)
length_adjust_Seq_number$Seq_length<-taxTab_WP2_12S_length$Seq_length

taxTab_WP2_12S_length$Seq_number<-NULL
taxTab_WP2_12S_length$Seq_length<-NULL

#produce the sequence length adjusted BLAST input file
write.csv(taxTab_WP2_12S_length,"length_filtered_taxTab_WP2_12S_BLAST_INPUT.fasta", row.names=FALSE,quote = FALSE)

#now BLAST this on the cluster using the "12S rRNA gene" database using -perc_identity 70 -qcov_hsp_perc 80 -evalue 1

#________________________________________________________

#Read in the length adjusted BLAST output file from the first BLAST of NCBI 
#you need objects created in the code before this (^) in order for it to work
#________________________________________________________

length_filtered_BLAST_OUTPUT_12S <- read_table2("length_filtered_BLAST_OUTPUT_12S.txt", 
                                                col_names = FALSE)

TAXONOMY_filtered_BLAST_OUTPUT_12S <- read_csv("TAXONOMY_filtered_BLAST_OUTPUT_12S.csv")
length_filtered_BLAST_OUTPUT_12S<-cbind(length_filtered_BLAST_OUTPUT_12S,TAXONOMY_filtered_BLAST_OUTPUT_12S)

length_adjust_Seq_number$length_adjust_Seq_number<-as.character(gsub(">", "", length_adjust_Seq_number$length_adjust_Seq_number))
length_adjust_Seq_number<-length_adjust_Seq_number %>% rename(Seq = "length_adjust_Seq_number")
length_filtered_BLAST_OUTPUT_12S<-length_filtered_BLAST_OUTPUT_12S %>% rename(Seq = "X1")

length_adjust_Seq_number <- length_adjust_Seq_number %>% full_join(length_filtered_BLAST_OUTPUT_12S, by = c("Seq"))

#remove ASVs that could not be assigned using broad BLAST
length_adjust_Seq_number<-na.omit(length_adjust_Seq_number)

#remove ASvs that were classified as non-fish
length_adjust_Seq_number_fish<-length_adjust_Seq_number %>% filter(class == "Actinopteri")
length_adjust_Seq_number_mammals<-length_adjust_Seq_number %>% filter(class == "Mammalia")

#Format for the second BLAST

BLAST_2_mito<-length_adjust_Seq_number_fish$Seq
BLAST_2_mito<-data.frame(BLAST_2_mito)
BLAST_2_mito$Seq<-length_adjust_Seq_number_fish$Seq_length

BLAST_2_mito$Seq_BLAST <-paste(">",BLAST_2_mito$BLAST_2_mito,BLAST_2_mito$Seq)
BLAST_2_mito$Seq_BLAST<- gsub("> ", ">", BLAST_2_mito$Seq_BLAST)
BLAST_2_mito$Seq_BLAST<- as.character(gsub(" ", "\n", BLAST_2_mito$Seq_BLAST))

BLAST_2_mito$Seq<-NULL
BLAST_2_mito$BLAST_2_mito<-NULL

#produce the sequence length adjusted BLAST input file
write.csv(BLAST_2_mito,"length_filtered_taxTab_WP2_12S_BLAST_INPUT_FISH_ONLY.fasta", row.names=FALSE,quote = FALSE)
remove(BLAST_2_mito, length_adjust_Seq_number)

#now BLAST this on the cluster using the mitogenome database using -perc_identity 90 -qcov_hsp_perc 90 -evalue 0.001

#________________________________________________________

#Read in fish only BLAST output using mitogenome database

#________________________________________________________


BLAST_OUTPUT_12S <- read_table2("BLAST_OUTPUT_FISH_ONLY_12S.txt", col_names = FALSE)
BLAST_OUTPUT_12S<-BLAST_OUTPUT_12S %>% rename(rn = X2 )

#Adding the taxonomy linked with the mitogenome Ids
mitogenome_annotations <- read_csv("mitogenome_annotations.csv")
BLAST_OUTPUT_12S$rn <- gsub("ref|", "", BLAST_OUTPUT_12S$rn)
BLAST_OUTPUT_12S$rn <- gsub("\\|", "", BLAST_OUTPUT_12S$rn)

BLAST_OUTPUT_12S <- BLAST_OUTPUT_12S %>% full_join(mitogenome_annotations, by = c("rn"))
BLAST_OUTPUT_12S<-na.omit(BLAST_OUTPUT_12S)

SEQ_lookup<-length_adjust_Seq_number_fish$Seq
SEQ_lookup<-data.frame(SEQ_lookup)
SEQ_lookup$Seq<-length_adjust_Seq_number_fish$Seq_length

SEQ_lookup<-SEQ_lookup %>% rename(X1 = SEQ_lookup )

BLAST_OUTPUT_12S <- BLAST_OUTPUT_12S %>% full_join(SEQ_lookup, by = c("X1"))

BLAST_OUTPUT_12S<-na.omit(BLAST_OUTPUT_12S)

write.table(BLAST_OUTPUT_12S,"BLAST_OUTPUT_12S_IDs_tax.txt",row.names = FALSE,col.names=FALSE,quote = FALSE)
remove(mitogenome_annotations,length_adjust_Seq_number_fish,length_filtered_BLAST_OUTPUT_12S,SEQ_lookup,
       TAXONOMY_filtered_BLAST_OUTPUT_12S,taxTab_WP2_12S,taxTab_WP2_12S_length)
gc()

seqtabNoC_WP2_12S <- read_table2("seqtabNoC_WP2_12S.txt")
seqtabNoC_WP2_12S_transpose<-t(seqtabNoC_WP2_12S)
seqtabNoC_WP2_12S_transpose<-as.data.frame(seqtabNoC_WP2_12S_transpose)
seqtabNoC_WP2_12S_transpose<-setDT(seqtabNoC_WP2_12S_transpose, keep.rownames = TRUE)
names(seqtabNoC_WP2_12S_transpose)[names(seqtabNoC_WP2_12S_transpose) == "rn"] <- "Seq"
seqtabNoC_WP2_12S_transpose[] <- lapply(seqtabNoC_WP2_12S_transpose, gsub, pattern='"', replacement='')
seqtabNoC_WP2_12S_transpose <- na.omit(transform(seqtabNoC_WP2_12S_transpose, Seq = c("Seq", Seq[-nrow(seqtabNoC_WP2_12S_transpose)])))
header.true <- function(seqtabNoC_WP2_12S_transpose) {
  names(seqtabNoC_WP2_12S_transpose) <- as.character(unlist(seqtabNoC_WP2_12S_transpose[1,]))
  seqtabNoC_WP2_12S_transpose[-1,]
}
seqtabNoC_WP2_12S_transpose<-header.true(seqtabNoC_WP2_12S_transpose)
remove(seqtabNoC_WP2_12S)

#REMOVING SAMPLES WITH LOW READ DEPTH
seqtabNoC_WP2_12S_transpose<-data.frame(seqtabNoC_WP2_12S_transpose)
i <- c(2:1224) 
seqtabNoC_WP2_12S_transpose[ , i] <- apply(seqtabNoC_WP2_12S_transpose[ , i], 2,  
                                           function(x) as.numeric(as.character(x)))

rarefy_col_sum<-colSums(seqtabNoC_WP2_12S_transpose[c(2:1224)])
rarefy_col_sum<-data.frame(rarefy_col_sum)
rarefy_col_sum<-setDT(rarefy_col_sum, keep.rownames = TRUE)[]
rarefy_col_sum<-data.frame(rarefy_col_sum)

rarefy_col_sum<-rarefy_col_sum[!(rarefy_col_sum$rarefy_col_sum>=1000),]
drop<-rarefy_col_sum$rn

seqtabNoC_WP2_12S_transpose<-seqtabNoC_WP2_12S_transpose[,!(names(seqtabNoC_WP2_12S_transpose) %in% drop)]

#Replace the literal sequences from the seqtabNoc table and replace them with the corresponding SEQ Ids

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

#CALCULATING 0.05% READ THREASHOLD FOR LATER FILTEIRNG
seqtabNoC_WP2_12S_tax_col_sum<-colSums(seqtabNoC_WP2_12S_tax[c(2:992)])
seqtabNoC_WP2_12S_tax_col_sum<-data.frame(seqtabNoC_WP2_12S_tax_col_sum)
seqtabNoC_WP2_12S_tax_col_sum$read_filter<-seqtabNoC_WP2_12S_tax_col_sum$seqtabNoC_WP2_12S_tax_col_sum *0.0005
seqtabNoC_WP2_12S_tax_col_sum<- seqtabNoC_WP2_12S_tax_col_sum %>% rownames_to_column("ID")

#REMOVE ASVS THAT CANNOT BE IDENTIFIED WITH BLAST
seqtabNoC_WP2_12S_tax<-seqtabNoC_WP2_12S_tax %>% drop_na(species)
seqtabNoC_WP2_12S_tax$X1<-NULL

#CONVERITNG ASV TABLE INTO LONG FORMAT
seqtabNoC_WP2_12S_tax_long <- gather(seqtabNoC_WP2_12S_tax, Sample, reads, 2:992, factor_key=TRUE)
seqtabNoC_WP2_12S_tax_long$reads<-   as.numeric(seqtabNoC_WP2_12S_tax_long$reads)

#READ IN METADATA
lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3_summary<-table(lofresh_metadata_WP2_3$Days,lofresh_metadata_WP2_3$SampleSite_code)
lofresh_metadata_WP2_3_summary<-data.frame(lofresh_metadata_WP2_3_summary)
write.csv(lofresh_metadata_WP2_3_summary,"lofresh_metadata_WP2_3_summary.csv")
seqtabNoC_WP2_12S_tax_long<-seqtabNoC_WP2_12S_tax_long %>% rename(ID=Sample)

#ADD METADATA TO ASV TABLE
seqtabNoC_WP2_12S_tax_long_meta <- merge(lofresh_metadata_WP2_3[, c("ID", "WP")],seqtabNoC_WP2_12S_tax_long, by="ID")
sum_reads<- sum(seqtabNoC_WP2_12S_tax_long_meta$reads)

seqtabNoC_WP2_12S_tax_long_meta_clean_filtered <- merge(seqtabNoC_WP2_12S_tax_col_sum[, c("ID", "read_filter")],seqtabNoC_WP2_12S_tax_long_meta, by="ID")

#REMOVE READS LESS THAN THE 0.005% CUT OFF
seqtabNoC_WP2_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2_12S_tax_long_meta_clean_filtered$reads<=seqtabNoC_WP2_12S_tax_long_meta_clean_filtered$read_filter),]
#REMOVE READS LESS THAN 20
seqtabNoC_WP2_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2_12S_tax_long_meta_clean_filtered$reads<20),]

write.csv(seqtabNoC_WP2_12S_tax_long_meta_clean_filtered,"WP2_ASV_table_long_filtered_family.csv")
seqtabNoC_WP2_12S_tax_long_meta_clean_filtered <- read_csv("WP2_ASV_table_long_filtered_family.csv")

#CONVERTING FILTERED LONG FORMAT INTO WIDE FORMAT
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

#NORMALISING READS

#YOU MAY NEED TO CHANGE NUMBERS HERE

seqtabNoC_WP2_12S_tax_col_sum2<-rowSums(seqtabNoC_WP2_12S_tax_wide_filtered[c(1:68)])
seqtabNoC_WP2_12S_tax_col_sum2<-data.frame(seqtabNoC_WP2_12S_tax_col_sum2)
seqtabNoC_WP2_12S_tax_col_sum2$ID <- rownames(seqtabNoC_WP2_12S_tax_col_sum2)

seqtabNoC_WP2_12S_tax_normalised<- merge(seqtabNoC_WP2_12S_tax_col_sum2[, c("ID", "seqtabNoC_WP2_12S_tax_col_sum2")],seqtabNoC_WP2_12S_tax_long_meta_clean_filtered, by="ID")

seqtabNoC_WP2_12S_tax_normalised <- transform(seqtabNoC_WP2_12S_tax_normalised, normalised_reads = reads / seqtabNoC_WP2_12S_tax_col_sum2)

ggplot(seqtabNoC_WP2_12S_tax_normalised , aes(x = ID, fill = species, y = normalised_reads)) + 
  geom_bar(stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))+ theme(legend.position = "none")

#CONVERTING NORMALISED LONG FORMAT INTO WIDE FORMAT
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

seqtabNoC_WP2C_12S_tax_long_meta_clean_Gwash<-seqtabNoC_WP2C_12S_tax_long_meta_clean %>% filter(Region == "Leicestershire")
seqtabNoC_WP2C_12S_tax_long_meta_clean_Towy<-seqtabNoC_WP2C_12S_tax_long_meta_clean %>% filter(Region == "Carmarthenshire")
seqtabNoC_WP2C_12S_tax_long_meta_clean_Skan<-seqtabNoC_WP2C_12S_tax_long_meta_clean %>% filter(Region == "New_York")
seqtabNoC_WP2C_12S_tax_long_meta_clean_Glatt<-seqtabNoC_WP2C_12S_tax_long_meta_clean %>% filter(Region == "Switzerland")

write.csv(seqtabNoC_WP2C_12S_tax_long_meta_clean_Gwash,"WP2C_Gwash_ASV_table_long_filtered_family.csv")

write.csv(seqtabNoC_WP2C_12S_tax_long_meta_clean_Towy,"WP2C_Towy_ASV_table_long_filtered_family.csv")

write.csv(seqtabNoC_WP2C_12S_tax_long_meta_clean_Skan,"WP2C_Skan_ASV_table_long_filtered_family.csv")

write.csv(seqtabNoC_WP2C_12S_tax_long_meta_clean_Glatt,"WP2C_Glatt_ASV_table_long_filtered_family.csv")

#________________________________________________________

########WP2A BLAST  TAXONOMY ANALYSIS######

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

#aggregating replicates
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% 
 group_by(SampleSite_time,species) %>% 
  summarize(reads=mean(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads)

seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% 
group_by(SampleSite_time,species) %>% 
summarize(reads=mean(reads)) %>%
rename(SampleSite_time=SampleSite_time,species=species, reads= reads)

#REMOVE CLEAR FISH THAT ARE LIEKLY CONSUMED BY HUAMNS AND NOT NATIVE
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
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species<-sub('Salmo_obtusirostris','Salmo_trutta',seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species)
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species<-sub('Anguilla_rostrata','Anguilla_anguilla',seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species)
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species<-sub('Salvelinus_fontinalis_x_Salvelinus_malma','Salvelinus_malma',seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species)

#Remove negative controls
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$SampleSite_time=="ENC_T10"),]
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$SampleSite_time=="ENC_T17"),]
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$SampleSite_time=="ENC_T18"),]
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$SampleSite_time=="ENC_T19"),]


#CONVERTING FILTERED LONG FORMAT INTO WIDE FORMAT
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

#rarefaction_curve<-seqtabNoC_WP2A_12S_tax_wide_filtered
#rownames(rarefaction_curve) <- NULL

#i <- c(1:16) 
#rarefaction_curve[ , i] <- apply(rarefaction_curve[ , i], 2,  
                                                  # function(x) as.integer(as.character(x)))
#(raremax <- min(rowSums(rarefaction_curve)))
#rarecurve(rarefaction_curve, col = c("red", "blue", "orange", "yellow", "green"), step = 20, sample = raremax, label = FALSE)

#write.csv(seqtabNoC_WP2A_12S_tax_wide_filtered,"WP2A_ASV_table_wide_filtered_replicates.csv")
write.csv(seqtabNoC_WP2A_12S_tax_wide_filtered,"WP2A_ASV_table_wide_filtered.csv")

#RELATIVE AMPLICON FREQUENCY ABUNDANCE
ggplot(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered, aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))


#NORMALISING READS
seqtabNoC_WP2A_12S_tax_col_sum2<-rowSums(seqtabNoC_WP2A_12S_tax_wide_filtered[c(1:16)])
WP2A_12S_sps_read_summary<-colSums(seqtabNoC_WP2A_12S_tax_wide_filtered[c(1:16)])
WP2A_12S_sps_read_summary<-data.frame(WP2A_12S_sps_read_summary)

seqtabNoC_WP2A_12S_tax_col_sum2<-data.frame(seqtabNoC_WP2A_12S_tax_col_sum2)
seqtabNoC_WP2A_12S_tax_col_sum2$SampleSite_time <- rownames(seqtabNoC_WP2A_12S_tax_col_sum2)

seqtabNoC_WP2A_12S_tax_normalised<- merge(seqtabNoC_WP2A_12S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2A_12S_tax_col_sum2")],seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered, by="SampleSite_time")

#GETTING THE SUM OF EACH SAMPLE FOR NORMALISING READS
seqtabNoC_WP2A_12S_tax_col_sum3<-colSums(seqtabNoC_WP2A_12S_tax_wide_filtered[c(1:16)])
seqtabNoC_WP2A_12S_tax_col_sum3<-data.frame(seqtabNoC_WP2A_12S_tax_col_sum3)

seqtabNoC_WP2A_12S_tax_col_sum3$species <- rownames(seqtabNoC_WP2A_12S_tax_col_sum3)

seqtabNoC_WP2A_12S_tax_normalised<- merge(seqtabNoC_WP2A_12S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2A_12S_tax_col_sum2")],seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered, by="SampleSite_time")
seqtabNoC_WP2A_12S_tax_top<- merge(seqtabNoC_WP2A_12S_tax_col_sum3[, c("species", "seqtabNoC_WP2A_12S_tax_col_sum3")],seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered, by="species")

seqtabNoC_WP2A_12S_tax_normalised <- transform(seqtabNoC_WP2A_12S_tax_normalised, normalised_reads = reads / seqtabNoC_WP2A_12S_tax_col_sum2)

ggplot(seqtabNoC_WP2A_12S_tax_normalised , aes(x = SampleSite_time, fill = species, y = normalised_reads)) + 
  geom_bar(stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")
##scale_fill_brewer(palette="Paired")

lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3_summary = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

seqtabNoC_WP2A_12S_tax_normalised <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "Days")],seqtabNoC_WP2A_12S_tax_normalised, by="SampleSite_time")
seqtabNoC_WP2A_12S_tax_normalised <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "SampleSite_code")],seqtabNoC_WP2A_12S_tax_normalised, by="SampleSite_time")

########SEASONAL IMPACT ON SPECIES DETECTION########
sps_seasonal_detection<-unique(seqtabNoC_WP2A_12S_tax_normalised[c("SampleSite_time", "species")])
sps_seasonal_detection <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "season")],sps_seasonal_detection, by="SampleSite_time")

sps_seasonal_detection_summary<-count(sps_seasonal_detection, season, species)
sps_seasonal_detection_summary$season <- factor(sps_seasonal_detection_summary$season, levels = c("Winter", "Spring", "Summer", "Autumn"))

#sps_seasonal_detection_summary<-subset(sps_seasonal_detection_summary,duplicated(species) | duplicated(species, fromLast=TRUE))

ggplot(data=sps_seasonal_detection_summary, aes(x=season, y=n)) +
  geom_line()+
  geom_point()+facet_grid((.~species))


#PLOTTING TOP Taxon
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

#CONVERTING NORMALISED LONG FORMAT INTO WIDE FORMAT
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

########PERMANOVA ON TAXON########
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


#PLOTTING BETA DIVERSITY NMDS
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

gg <- merge(data.scores_meta_data,aggregate(cbind(NMDS1.x=NMDS1,NMDS2.y=NMDS2)~season,data.scores_meta_data,mean),by="season")
ggplot(gg, aes(x = NMDS1, y = NMDS2,colour = Dist_from_lake)) + scale_shape_manual(values=c(15,16,17,8,18))+
  geom_point(size = 4, aes(shape=season,colour = Dist_from_lake))+
  theme_bw()+  theme(axis.title.x=element_blank())+  theme(axis.title.y=element_blank())+ 
  geom_segment(aes(x=NMDS1.x, y=NMDS2.y, xend=NMDS1, yend=NMDS2))+
  geom_point(aes(x=NMDS1.x,y=NMDS2.y),size=3,colour="red")


########ALPHA DIVERSITY METRICS########
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

#write.csv(alpha_12S_WP2A,"alpha_12S_WP2A_replicates.csv")
write.csv(alpha_12S_WP2A,"alpha_12S_WP2A.csv")

######## WP2A alpha diversity modeling ##########
library(mgcv)

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

#Linear models
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


######## Sparse partial least squares analysis ##########
library(mixOmics)

NORM_WP2A_ASV_table_wide_filtered <- read_csv("NORM_WP2A_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2A_ASV_table_wide_filtered
#REMOVE OUTLIERS IDENTIFIED BY THE NMDS PLOT
#adonis_data<-adonis_data[adonis_data$X1 != "E04_T07", ]  

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
#result.spls <- spls(adonis_data, adonis_data_meta, ncomp = ncomp,  mode = 'regression')
result.spls <- spls(adonis_data, adonis_data_meta, ncomp = ncomp,  mode = 'canonical')

design <- data.frame(samp = adonis_data_meta$sample)

tune.spls <- perf(result.spls, validation = 'Mfold', folds = 5,
                  criterion = 'all',nrepeat = 99, progressBar = TRUE)

tune.spls$R2
plot(tune.spls$Q2.total)
abline(h = 0.0975)

plot(tune.spls, criterion="R2", type = 'l')
plot(tune.spls, criterion="Q2", type = 'l')
plot(tune.spls)

cim_plot<-cim(result.spls,
              comp = 1:1,
              margins = c(15, 15))

######## Mantel test ##########
adonis_data_meta$Days<-NULL
adonis_data_meta$Dist_from_lake<-NULL

dist.abund = vegdist(adonis_data, method = "bray")
dist.meta = dist(adonis_data_meta, method = "euclidean")

abund_meta = mantel(dist.abund, dist.meta, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_meta

aa = as.vector(dist.abund)
mm = as.vector(dist.meta)
mat = data.frame(aa,mm)

ggplot(mat, aes(y = aa, x = mm)) + 
  geom_point(size = 1) + 
  theme_bw()+ geom_smooth(method = lm)+ scale_y_continuous(limits = c(0, 1))

#________________________________________________________

########WP2C_Glatt BLAST  TAXONOMY ANALYSIS########

#________________________________________________________

setwd("WP2_analysis/12S_BLAST")
seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean_filtered <- read_csv("WP2C_Glatt_ASV_table_long_filtered_family.csv")

#aggregating replicates
seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=mean(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads)

#REMOVE CLEAR FISH THAT ARE LIEKLY CONSUMED BY HUAMNS AND NOT NATIVE
seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean_filtered$species=="Sardina_pilchardus"),]

#CONVERTING FILTERED LONG FORMAT INTO WSampleSite_timeE FORMAT
seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered<-seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean_filtered
seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered$read_filter <- NULL 
seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered$Days <- NULL 
seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered$Date_sampled <- NULL 
seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered$SampleSite_code <- NULL 
seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered$WP <- NULL 
seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered$Seq_length <- NULL 

seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered<-seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=sum(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads,)

seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered <- spread(seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered,SampleSite_time,reads)
seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered[is.na(seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered)] <- 0
seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered <- data.frame(t(seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered))

names(seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered) <- as.matrix(seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered[1, ])
seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered <- seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered[-1, ]
seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered[] <- lapply(seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered, function(x) type.convert(as.character(x)))

#write.csv(seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered,"WP2C_Glatt_ASV_table_wide_filtered_replicates.csv")
write.csv(seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered,"WP2C_Glatt_ASV_table_wide_filtered.csv")

#RELATIVE AMPLICON FREQUENCY ABUNDANCE
ggplot(seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean_filtered, aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))

#ABSOLUTE AMPLICON ABUNDANCE
#ggplot(seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean , aes(x = SampleSite_time, fill = species, y = reads)) + 
# geom_bar(stat = "identity", colour = "black")+
#theme(axis.text.x = element_text(angle = 90),legend.position = "none")
##scale_fill_brewer(palette="Paired")

##SPLITTING BY SPECIES
Sardina_pilchardus<-seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean_filtered %>% filter(species == "Sardina_pilchardus")

Anguilla_anguilla_sum_reads<-Anguilla_anguilla %>% 
  group_by(SampleSite_time) %>% 
  summarize(reads=sum(reads)) %>%
  rename(SampleSite_time_sum =SampleSite_time,reads_sum= reads)

ggplot(data=Anguilla_anguilla_sum_reads, aes(x=SampleSite_time_sum, y=reads_sum))+
  geom_col()

#Venn diagram
#seqtabNoC_WP2C_Glatt_12S_Venn<-seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean[!(seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean$reads==0),]
#hist(seqtabNoC_WP2C_Glatt_12S_Venn$reads, breaks=1000,xlim=c(0,50000))

#E01_Venn<-seqtabNoC_WP2C_Glatt_12S_Venn %>% filter(SampleSite_code == "E01")
#E01_Venn<-E01_Venn$species

#E02_Venn<-seqtabNoC_WP2C_Glatt_12S_Venn %>% filter(SampleSite_code == "E02")
#E02_Venn<-E02_Venn$species

#E03_Venn<-seqtabNoC_WP2C_Glatt_12S_Venn %>% filter(SampleSite_code == "E03")
#E03_Venn<-E03_Venn$species

#E04_Venn<-seqtabNoC_WP2C_Glatt_12S_Venn %>% filter(SampleSite_code == "E04")
#E04_Venn<-E04_Venn$species

#E05_Venn<-seqtabNoC_WP2C_Glatt_12S_Venn %>% filter(SampleSite_code == "E05")
#E05_Venn<-E05_Venn$species


#venn.diagram(
# x = list(E01_Venn, E02_Venn,E03_Venn,E04_Venn,E05_Venn),
#category.names = c("E01" , "E02","E03","E04","E05"),
#filename = 'WP2C_Glatt_venn_diagramm.png',
#output=TRUE)


#NORMALISING READS

seqtabNoC_WP2C_Glatt_12S_tax_col_sum2<-rowSums(seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered[c(1:25)])
seqtabNoC_WP2C_Glatt_12S_tax_col_sum2<-data.frame(seqtabNoC_WP2C_Glatt_12S_tax_col_sum2)
seqtabNoC_WP2C_Glatt_12S_tax_col_sum2$SampleSite_time <- rownames(seqtabNoC_WP2C_Glatt_12S_tax_col_sum2)

seqtabNoC_WP2C_Glatt_12S_tax_normalised<- merge(seqtabNoC_WP2C_Glatt_12S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2C_Glatt_12S_tax_col_sum2")],seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean_filtered, by="SampleSite_time")

#GETTING THE SUM OF EACH SAMPLE FOR NORMALISING READS
seqtabNoC_WP2C_Glatt_12S_tax_col_sum3<-colSums(seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered[c(1:25)])
seqtabNoC_WP2C_Glatt_12S_tax_col_sum3<-data.frame(seqtabNoC_WP2C_Glatt_12S_tax_col_sum3)

seqtabNoC_WP2C_Glatt_12S_tax_col_sum3$species <- rownames(seqtabNoC_WP2C_Glatt_12S_tax_col_sum3)

seqtabNoC_WP2C_Glatt_12S_tax_normalised<- merge(seqtabNoC_WP2C_Glatt_12S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2C_Glatt_12S_tax_col_sum2")],seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean_filtered, by="SampleSite_time")
seqtabNoC_WP2C_Glatt_12S_tax_top<- merge(seqtabNoC_WP2C_Glatt_12S_tax_col_sum3[, c("species", "seqtabNoC_WP2C_Glatt_12S_tax_col_sum3")],seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean_filtered, by="species")

seqtabNoC_WP2C_Glatt_12S_tax_normalised <- transform(seqtabNoC_WP2C_Glatt_12S_tax_normalised, normalised_reads = reads / seqtabNoC_WP2C_Glatt_12S_tax_col_sum2)

ggplot(seqtabNoC_WP2C_Glatt_12S_tax_normalised , aes(x = SampleSite_time, fill = species, y = normalised_reads)) + 
  geom_bar(stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")
##scale_fill_brewer(palette="Paired")


#PLOTTING TOP Taxon
seqtabNoC_WP2C_Glatt_12S_tax_top$SampleSite_time<-factor(seqtabNoC_WP2C_Glatt_12S_tax_top$SampleSite_time,levels= c('GL01',
                         'GL13','GL02','GL03','GL04','GL05','GL06','GL07','GL08','GL09','GL10','GL11','GL12'),ordered=TRUE)

seqtabNoC_WP2C_Glatt_12S_tax_top %>%
  mutate(species = fct_reorder(species, desc(seqtabNoC_WP2C_Glatt_12S_tax_col_sum3))) %>%
  ggplot(aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))+scale_fill_brewer(palette="Paired")+ 
  theme(legend.text=element_text(size=10))+
  theme(legend.key.size = unit(0.3, 'cm'))

#write.csv(seqtabNoC_WP2C_Glatt_12S_tax_normalised,"NORM_WP2C_Glatt_ASV_table_LONG_filtered_replicates.csv")
write.csv(seqtabNoC_WP2C_Glatt_12S_tax_normalised,"NORM_WP2C_Glatt_ASV_table_LONG_filtered.csv")

#CONVERTING NORMALISED LONG FORMAT INTO SampleSite_time FORMAT
seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide<-seqtabNoC_WP2C_Glatt_12S_tax_normalised
seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide$read_filter <- NULL 
seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide$Days <- NULL 
seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide$Date_sampled <- NULL 
seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide$SampleSite_code <- NULL 
seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide$WP <- NULL 
seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide$Seq_length <- NULL 
seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide$reads <- NULL 
seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide$seqtabNoC_WP2_12S_tax_col_sum2 <- NULL 


seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide<-seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, normalised_reads= normalised_reads,)

seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide <- spread(seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide,SampleSite_time,normalised_reads)
seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide[is.na(seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide)] <- 0
seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide <- data.frame(t(seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide))

names(seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide) <- as.matrix(seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide[1, ])
seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide <- seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide[-1, ]
seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide[] <- lapply(seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide, function(x) type.convert(as.character(x)))

#write.csv(seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide,"NORM_WP2C_Glatt_ASV_table_wide_filtered_replicates.csv")
write.csv(seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide,"NORM_WP2C_Glatt_ASV_table_wide_filtered.csv")


#PERMANOVA ON TAXON
NORM_WP2C_Glatt_ASV_table_wide_filtered <- read_csv("NORM_WP2C_Glatt_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2C_Glatt_ASV_table_wide_filtered
adonis_data_samples<-adonis_data$X1
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(SampleSite_time = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = X1)

lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","WP","Days","Catchment_flow_m3_s","Average_Soil_Temp_5cm_oC","Total_Rainfall","Dist_from_lake","pH","conductivity_uS_cm","gran_alk_ueq_L")],adonis_data, by="SampleSite_time")

adonis_data<-na.omit(adonis_data)
adonis_data_meta<- adonis_data[c(1:10)]
adonis_data<-adonis_data[,-c(1:10)]

adonis_data_meta$Days<-as.numeric(adonis_data_meta$Days)

adon.results_WP2<-adonis(adonis_data ~ adonis_data_meta$Days+
                           adonis_data_meta$Dist_from_lake+
                           adonis_data_meta$Average_Soil_Temp_5cm_oC+
                           adonis_data_meta$Catchment_flow_m3_s+
                           adonis_data_meta$gran_alk_ueq_L+
                           adonis_data_meta$Total_Rainfall+
                           adonis_data_meta$conductivity_uS_cm+
                           adonis_data_meta$pH,
                         method="bray",perm=999)

print(adon.results_WP2)

#PLOTTING BETA DIVERSITY NMDS
abund_table_12S_SampleSite_time <- read.csv("NORM_WP2C_Glatt_ASV_table_wide_filtered.csv", header = TRUE)
abund_table_12S <- abund_table_12S_SampleSite_time[,-1]
com <- abund_table_12S[,col(abund_table_12S)]
m_com <- as.matrix(com)
set.seed(123)
nmds = metaMDS(m_com, distance = "bray", try = 1000, trymax = 1000, k = 3)

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds))
#add columns to data frame 
data.scores$SampleSite_time = abund_table_12S_SampleSite_time[,1]

lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]
data.scores_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","SampleSite_code","Days","Dist_from_lake")],data.scores, by="SampleSite_time")


p<-ggplot(data.scores_meta_data, aes(x = NMDS1, y = NMDS2,colour = Dist_from_lake)) + 
  geom_point(size = 2, aes(colour = Dist_from_lake))+
  theme_bw()
ggplotly(p)


#ALPHA DIVERSITY METRICS
shannon_index<-diversity(seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide, index = "shannon")
shannon_index<-data.frame(shannon_index)
ID <- rownames(shannon_index)
shannon_index <- cbind(shannon_index,ID)

simpson_index<-diversity(seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide, index = "simpson")
simpson_index<-data.frame(simpson_index)
ID <- rownames(simpson_index)
simpson_index <- cbind(simpson_index,ID)

species_count <- rowSums(seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide >0)
species_count<-data.frame(species_count)
species_count$ID <- rownames(species_count)

alpha_12S_WP2C_Glatt<-merge(shannon_index,simpson_index,by="ID")
alpha_12S_WP2C_Glatt<-merge(alpha_12S_WP2C_Glatt,species_count,by="ID")

#write.csv(alpha_12S_WP2C_Glatt,"alpha_12S_WP2C_Glatt_replicates.csv")
write.csv(alpha_12S_WP2C_Glatt,"alpha_12S_WP2C_Glatt.csv")

#________________________________________________________

########WP2C_Towy BLAST  TAXONOMY ANALYSIS########

#________________________________________________________

setwd("WP2_analysis/12S_BLAST")
seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean_filtered <- read_csv("WP2C_Towy_ASV_table_long_filtered_family.csv")

#aggregating replicates
seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=mean(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads)

#REMOVE CLEAR FISH THAT ARE LIEKLY CONSUMED BY HUAMNS AND NOT NATIVE
seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean_filtered$species=="Platichthys_stellatus"),]

#CONVERTING FILTERED LONG FORMAT INTO WSampleSite_timeE FORMAT
seqtabNoC_WP2C_Towy_12S_tax_wide_filtered<-seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean_filtered
seqtabNoC_WP2C_Towy_12S_tax_wide_filtered$read_filter <- NULL 
seqtabNoC_WP2C_Towy_12S_tax_wide_filtered$Days <- NULL 
seqtabNoC_WP2C_Towy_12S_tax_wide_filtered$Date_sampled <- NULL 
seqtabNoC_WP2C_Towy_12S_tax_wide_filtered$SampleSite_code <- NULL 
seqtabNoC_WP2C_Towy_12S_tax_wide_filtered$WP <- NULL 
seqtabNoC_WP2C_Towy_12S_tax_wide_filtered$Seq_length <- NULL 

seqtabNoC_WP2C_Towy_12S_tax_wide_filtered<-seqtabNoC_WP2C_Towy_12S_tax_wide_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=sum(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads,)

seqtabNoC_WP2C_Towy_12S_tax_wide_filtered <- spread(seqtabNoC_WP2C_Towy_12S_tax_wide_filtered,SampleSite_time,reads)
seqtabNoC_WP2C_Towy_12S_tax_wide_filtered[is.na(seqtabNoC_WP2C_Towy_12S_tax_wide_filtered)] <- 0
seqtabNoC_WP2C_Towy_12S_tax_wide_filtered <- data.frame(t(seqtabNoC_WP2C_Towy_12S_tax_wide_filtered))

names(seqtabNoC_WP2C_Towy_12S_tax_wide_filtered) <- as.matrix(seqtabNoC_WP2C_Towy_12S_tax_wide_filtered[1, ])
seqtabNoC_WP2C_Towy_12S_tax_wide_filtered <- seqtabNoC_WP2C_Towy_12S_tax_wide_filtered[-1, ]
seqtabNoC_WP2C_Towy_12S_tax_wide_filtered[] <- lapply(seqtabNoC_WP2C_Towy_12S_tax_wide_filtered, function(x) type.convert(as.character(x)))

#write.csv(seqtabNoC_WP2C_Towy_12S_tax_wide_filtered,"WP2C_Towy_ASV_table_wide_filtered_replicates.csv")
write.csv(seqtabNoC_WP2C_Towy_12S_tax_wide_filtered,"WP2C_Towy_ASV_table_wide_filtered.csv")

#RELATIVE AMPLICON FREQUENCY ABUNDANCE
ggplot(seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean_filtered, aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))

#ABSOLUTE AMPLICON ABUNDANCE
#ggplot(seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean , aes(x = SampleSite_time, fill = species, y = reads)) + 
# geom_bar(stat = "identity", colour = "black")+
#theme(axis.text.x = element_text(angle = 90),legend.position = "none")
##scale_fill_brewer(palette="Paired")

##SPLITTING BY SPECIES
Sardina_pilchardus<-seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean_filtered %>% filter(species == "Sardina_pilchardus")

Anguilla_anguilla_sum_reads<-Anguilla_anguilla %>% 
  group_by(SampleSite_time) %>% 
  summarize(reads=sum(reads)) %>%
  rename(SampleSite_time_sum =SampleSite_time,reads_sum= reads)

ggplot(data=Anguilla_anguilla_sum_reads, aes(x=SampleSite_time_sum, y=reads_sum))+
  geom_col()

#Venn diagram
#seqtabNoC_WP2C_Towy_12S_Venn<-seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean[!(seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean$reads==0),]
#hist(seqtabNoC_WP2C_Towy_12S_Venn$reads, breaks=1000,xlim=c(0,50000))

#E01_Venn<-seqtabNoC_WP2C_Towy_12S_Venn %>% filter(SampleSite_code == "E01")
#E01_Venn<-E01_Venn$species

#E02_Venn<-seqtabNoC_WP2C_Towy_12S_Venn %>% filter(SampleSite_code == "E02")
#E02_Venn<-E02_Venn$species

#E03_Venn<-seqtabNoC_WP2C_Towy_12S_Venn %>% filter(SampleSite_code == "E03")
#E03_Venn<-E03_Venn$species

#E04_Venn<-seqtabNoC_WP2C_Towy_12S_Venn %>% filter(SampleSite_code == "E04")
#E04_Venn<-E04_Venn$species

#E05_Venn<-seqtabNoC_WP2C_Towy_12S_Venn %>% filter(SampleSite_code == "E05")
#E05_Venn<-E05_Venn$species


#venn.diagram(
# x = list(E01_Venn, E02_Venn,E03_Venn,E04_Venn,E05_Venn),
#category.names = c("E01" , "E02","E03","E04","E05"),
#filename = 'WP2C_Towy_venn_diagramm.png',
#output=TRUE)


#NORMALISING READS
seqtabNoC_WP2C_Towy_12S_tax_col_sum2<-rowSums(seqtabNoC_WP2C_Towy_12S_tax_wide_filtered[c(1:10)])
seqtabNoC_WP2C_Towy_12S_tax_col_sum2<-data.frame(seqtabNoC_WP2C_Towy_12S_tax_col_sum2)
seqtabNoC_WP2C_Towy_12S_tax_col_sum2$SampleSite_time <- rownames(seqtabNoC_WP2C_Towy_12S_tax_col_sum2)

seqtabNoC_WP2C_Towy_12S_tax_normalised<- merge(seqtabNoC_WP2C_Towy_12S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2C_Towy_12S_tax_col_sum2")],seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean_filtered, by="SampleSite_time")

#GETTING THE SUM OF EACH SAMPLE FOR NORMALISING READS
seqtabNoC_WP2C_Towy_12S_tax_col_sum3<-colSums(seqtabNoC_WP2C_Towy_12S_tax_wide_filtered[c(1:10)])
seqtabNoC_WP2C_Towy_12S_tax_col_sum3<-data.frame(seqtabNoC_WP2C_Towy_12S_tax_col_sum3)

seqtabNoC_WP2C_Towy_12S_tax_col_sum3$species <- rownames(seqtabNoC_WP2C_Towy_12S_tax_col_sum3)

seqtabNoC_WP2C_Towy_12S_tax_normalised<- merge(seqtabNoC_WP2C_Towy_12S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2C_Towy_12S_tax_col_sum2")],seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean_filtered, by="SampleSite_time")
seqtabNoC_WP2C_Towy_12S_tax_top<- merge(seqtabNoC_WP2C_Towy_12S_tax_col_sum3[, c("species", "seqtabNoC_WP2C_Towy_12S_tax_col_sum3")],seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean_filtered, by="species")

seqtabNoC_WP2C_Towy_12S_tax_normalised <- transform(seqtabNoC_WP2C_Towy_12S_tax_normalised, normalised_reads = reads / seqtabNoC_WP2C_Towy_12S_tax_col_sum2)

ggplot(seqtabNoC_WP2C_Towy_12S_tax_normalised , aes(x = SampleSite_time, fill = species, y = normalised_reads)) + 
  geom_bar(stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))
##scale_fill_brewer(palette="Paired")


#PLOTTING TOP Taxon
seqtabNoC_WP2C_Towy_12S_tax_top %>%
  mutate(species = fct_reorder(species, desc(seqtabNoC_WP2C_Towy_12S_tax_col_sum3))) %>%
  ggplot(aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))+scale_fill_brewer(palette="Paired")+ 
  theme(legend.text=element_text(size=10))+
  theme(legend.key.size = unit(0.3, 'cm'))

#write.csv(seqtabNoC_WP2C_Towy_12S_tax_normalised,"NORM_WP2C_Towy_ASV_table_LONG_filtered_replicates.csv")
write.csv(seqtabNoC_WP2C_Towy_12S_tax_normalised,"NORM_WP2C_Towy_ASV_table_LONG_filtered.csv")


#CONVERTING NORMALISED LONG FORMAT INTO SampleSite_time FORMAT
seqtabNoC_WP2C_Towy_12S_tax_normalised_wide<-seqtabNoC_WP2C_Towy_12S_tax_normalised
seqtabNoC_WP2C_Towy_12S_tax_normalised_wide$read_filter <- NULL 
seqtabNoC_WP2C_Towy_12S_tax_normalised_wide$Days <- NULL 
seqtabNoC_WP2C_Towy_12S_tax_normalised_wide$Date_sampled <- NULL 
seqtabNoC_WP2C_Towy_12S_tax_normalised_wide$SampleSite_code <- NULL 
seqtabNoC_WP2C_Towy_12S_tax_normalised_wide$WP <- NULL 
seqtabNoC_WP2C_Towy_12S_tax_normalised_wide$Seq_length <- NULL 
seqtabNoC_WP2C_Towy_12S_tax_normalised_wide$reads <- NULL 
seqtabNoC_WP2C_Towy_12S_tax_normalised_wide$seqtabNoC_WP2_12S_tax_col_sum2 <- NULL 


seqtabNoC_WP2C_Towy_12S_tax_normalised_wide<-seqtabNoC_WP2C_Towy_12S_tax_normalised_wide %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, normalised_reads= normalised_reads,)

seqtabNoC_WP2C_Towy_12S_tax_normalised_wide <- spread(seqtabNoC_WP2C_Towy_12S_tax_normalised_wide,SampleSite_time,normalised_reads)
seqtabNoC_WP2C_Towy_12S_tax_normalised_wide[is.na(seqtabNoC_WP2C_Towy_12S_tax_normalised_wide)] <- 0
seqtabNoC_WP2C_Towy_12S_tax_normalised_wide <- data.frame(t(seqtabNoC_WP2C_Towy_12S_tax_normalised_wide))

names(seqtabNoC_WP2C_Towy_12S_tax_normalised_wide) <- as.matrix(seqtabNoC_WP2C_Towy_12S_tax_normalised_wide[1, ])
seqtabNoC_WP2C_Towy_12S_tax_normalised_wide <- seqtabNoC_WP2C_Towy_12S_tax_normalised_wide[-1, ]
seqtabNoC_WP2C_Towy_12S_tax_normalised_wide[] <- lapply(seqtabNoC_WP2C_Towy_12S_tax_normalised_wide, function(x) type.convert(as.character(x)))

#write.csv(seqtabNoC_WP2C_Towy_12S_tax_normalised_wide,"NORM_WP2C_Towy_ASV_table_wide_filtered_replicates.csv")
write.csv(seqtabNoC_WP2C_Towy_12S_tax_normalised_wide,"NORM_WP2C_Towy_ASV_table_wide_filtered.csv")


#PERMANOVA ON TAXON
NORM_WP2C_Towy_ASV_table_wide_filtered <- read_csv("NORM_WP2C_Towy_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2C_Towy_ASV_table_wide_filtered
adonis_data_samples<-adonis_data$X1
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(ID = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = X1)

lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","WP","Days","Catchment_flow_m3_s","Average_Soil_Temp_5cm_oC","Total_Rainfall","Dist_from_lake","pH","conductivity_uS_cm","gran_alk_ueq_L")],adonis_data, by="SampleSite_time")

adonis_data<-na.omit(adonis_data)
adonis_data_meta<- adonis_data[c(1:10)]
adonis_data<-adonis_data[,-c(1:10)]

adonis_data_meta$Days<-as.numeric(adonis_data_meta$Days)

adon.results_WP2<-adonis(adonis_data ~ adonis_data_meta$Days+
                           adonis_data_meta$Dist_from_lake+
                           adonis_data_meta$Catchment_flow_m3_s+
                           adonis_data_meta$gran_alk_ueq_L+
                           adonis_data_meta$Total_Rainfall+
                           adonis_data_meta$conductivity_uS_cm+
                           adonis_data_meta$pH,
                         method="bray",perm=999)

#PLOTTING BETA DIVERSITY NMDS
abund_table_12S_ID <- read.csv("NORM_WP2C_Towy_ASV_table_wide_filtered.csv", header = TRUE)
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
data.scores_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","SampleSite_code","Days","Dist_from_lake")],data.scores, by="SampleSite_time")


p<-ggplot(data.scores_meta_data, aes(x = NMDS1, y = NMDS2,colour = Dist_from_lake)) + 
  geom_point(size = 2, aes(colour = Dist_from_lake))+
  theme_bw()
ggplotly(p)


#ALPHA DIVERSITY METRICS
shannon_index<-diversity(seqtabNoC_WP2C_Towy_12S_tax_normalised_wide, index = "shannon")
shannon_index<-data.frame(shannon_index)
ID <- rownames(shannon_index)
shannon_index <- cbind(shannon_index,ID)

simpson_index<-diversity(seqtabNoC_WP2C_Towy_12S_tax_normalised_wide, index = "simpson")
simpson_index<-data.frame(simpson_index)
ID <- rownames(simpson_index)
simpson_index <- cbind(simpson_index,ID)

species_count <- rowSums(seqtabNoC_WP2C_Towy_12S_tax_normalised_wide >0)
species_count<-data.frame(species_count)
species_count$ID <- rownames(species_count)

alpha_12S_WP2C_Towy<-merge(shannon_index,simpson_index,by="ID")
alpha_12S_WP2C_Towy<-merge(alpha_12S_WP2C_Towy,species_count,by="ID")

#write.csv(alpha_12S_WP2C_Towy,"alpha_12S_WP2C_Towy_replicates.csv")
write.csv(alpha_12S_WP2C_Towy,"alpha_12S_WP2C_Towy.csv")


#________________________________________________________

########WP2C_Gwash BLAST  TAXONOMY ANALYSIS#########

#________________________________________________________

setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/12S_BLAST")
seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean_filtered <- read_csv("WP2C_Gwash_ASV_table_long_filtered_family.csv")

#aggregating replicates
seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean_filtered %>% 
  group_by(SampleSite_time,species) %>% 
 summarize(reads=mean(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads)

#REMOVE CLEAR FISH THAT ARE LIEKLY CONSUMED BY HUAMNS AND NOT NATIVE
seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean_filtered$species=="Sardina_pilchardus"),]

#CONVERTING FILTERED LONG FORMAT INTO WIDE FORMAT
seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered<-seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean_filtered
seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered$read_filter <- NULL 
seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered$Days <- NULL 
seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered$Date_sampled <- NULL 
seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered$SampleSite_code <- NULL 
seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered$WP <- NULL 
seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered$Seq_length <- NULL 

seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered<-seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=sum(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads,)

seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered <- spread(seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered,SampleSite_time,reads)
seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered[is.na(seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered)] <- 0
seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered <- data.frame(t(seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered))

names(seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered) <- as.matrix(seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered[1, ])
seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered <- seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered[-1, ]
seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered[] <- lapply(seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered, function(x) type.convert(as.character(x)))

#write.csv(seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered,"WP2C_Gwash_ASV_table_wide_filtered_replicates.csv")
write.csv(seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered,"WP2C_Gwash_ASV_table_wide_filtered.csv")


#RELATIVE AMPLICON FREQUENCY ABUNDANCE
ggplot(seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean_filtered, aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))

#ABSOLUTE AMPLICON ABUNDANCE
#ggplot(seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean , aes(x = SampleSite_time, fill = species, y = reads)) + 
# geom_bar(stat = "identity", colour = "black")+
#theme(axis.text.x = element_text(angle = 90),legend.position = "none")
##scale_fill_brewer(palette="Paired")

seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered_transpose<-t(seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered)

##SPLITTING BY SPECIES
Anguilla_anguilla<-seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean_filtered %>% filter(species == "Anguilla_anguilla")

Sardina_pilchardus<-seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean_filtered %>% filter(species == "Sardina_pilchardus")

Anguilla_anguilla_sum_reads<-Anguilla_anguilla %>% 
  group_by(SampleSite_time) %>% 
  summarize(reads=sum(reads)) %>%
  rename(SampleSite_time_sum =SampleSite_time,reads_sum= reads)

ggplot(data=Anguilla_anguilla_sum_reads, aes(x=SampleSite_time_sum, y=reads_sum))+
  geom_col()

#Venn diagram
#seqtabNoC_WP2C_Gwash_12S_Venn<-seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean[!(seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean$reads==0),]
#hist(seqtabNoC_WP2C_Gwash_12S_Venn$reads, breaks=1000,xlim=c(0,50000))

#E01_Venn<-seqtabNoC_WP2C_Gwash_12S_Venn %>% filter(SampleSite_code == "E01")
#E01_Venn<-E01_Venn$species

#E02_Venn<-seqtabNoC_WP2C_Gwash_12S_Venn %>% filter(SampleSite_code == "E02")
#E02_Venn<-E02_Venn$species

#E03_Venn<-seqtabNoC_WP2C_Gwash_12S_Venn %>% filter(SampleSite_code == "E03")
#E03_Venn<-E03_Venn$species

#E04_Venn<-seqtabNoC_WP2C_Gwash_12S_Venn %>% filter(SampleSite_code == "E04")
#E04_Venn<-E04_Venn$species

#E05_Venn<-seqtabNoC_WP2C_Gwash_12S_Venn %>% filter(SampleSite_code == "E05")
#E05_Venn<-E05_Venn$species


#venn.diagram(
# x = list(E01_Venn, E02_Venn,E03_Venn,E04_Venn,E05_Venn),
#category.names = c("E01" , "E02","E03","E04","E05"),
#filename = 'WP2C_Gwash_venn_diagramm.png',
#output=TRUE)


#NORMALISING READS
seqtabNoC_WP2C_Gwash_12S_tax_col_sum2<-rowSums(seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered[c(1:23)])
seqtabNoC_WP2C_Gwash_12S_tax_col_sum2<-data.frame(seqtabNoC_WP2C_Gwash_12S_tax_col_sum2)
seqtabNoC_WP2C_Gwash_12S_tax_col_sum2$SampleSite_time <- rownames(seqtabNoC_WP2C_Gwash_12S_tax_col_sum2)

seqtabNoC_WP2C_Gwash_12S_tax_normalised<- merge(seqtabNoC_WP2C_Gwash_12S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2C_Gwash_12S_tax_col_sum2")],seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean_filtered, by="SampleSite_time")

#GETTING THE SUM OF EACH SAMPLE FOR NORMALISING READS
seqtabNoC_WP2C_Gwash_12S_tax_col_sum3<-colSums(seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered[c(1:23)])
seqtabNoC_WP2C_Gwash_12S_tax_col_sum3<-data.frame(seqtabNoC_WP2C_Gwash_12S_tax_col_sum3)

seqtabNoC_WP2C_Gwash_12S_tax_col_sum3$species <- rownames(seqtabNoC_WP2C_Gwash_12S_tax_col_sum3)

seqtabNoC_WP2C_Gwash_12S_tax_normalised<- merge(seqtabNoC_WP2C_Gwash_12S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2C_Gwash_12S_tax_col_sum2")],seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean_filtered, by="SampleSite_time")
seqtabNoC_WP2C_Gwash_12S_tax_top<- merge(seqtabNoC_WP2C_Gwash_12S_tax_col_sum3[, c("species", "seqtabNoC_WP2C_Gwash_12S_tax_col_sum3")],seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean_filtered, by="species")

seqtabNoC_WP2C_Gwash_12S_tax_normalised <- transform(seqtabNoC_WP2C_Gwash_12S_tax_normalised, normalised_reads = reads / seqtabNoC_WP2C_Gwash_12S_tax_col_sum2)

ggplot(seqtabNoC_WP2C_Gwash_12S_tax_normalised , aes(x = SampleSite_time, fill = species, y = normalised_reads)) + 
  geom_bar(stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))
##scale_fill_brewer(palette="Paired")


#PLOTTING TOP Taxon
seqtabNoC_WP2C_Gwash_12S_tax_top %>%
  mutate(species = fct_reorder(species, desc(seqtabNoC_WP2C_Gwash_12S_tax_col_sum3))) %>%
  ggplot(aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))+scale_fill_brewer(palette="Paired")+ 
  theme(legend.text=element_text(size=10))+
  theme(legend.key.size = unit(0.3, 'cm'))


#write.csv(seqtabNoC_WP2C_Gwash_12S_tax_normalised,"NORM_WP2C_Gwash_ASV_table_LONG_filtered_replicates.csv")
write.csv(seqtabNoC_WP2C_Gwash_12S_tax_normalised,"NORM_WP2C_Gwash_ASV_table_LONG_filtered.csv")


#CONVERTING NORMALISED LONG FORMAT INTO SampleSite_time FORMAT
seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide<-seqtabNoC_WP2C_Gwash_12S_tax_normalised
seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide$read_filter <- NULL 
seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide$Days <- NULL 
seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide$Date_sampled <- NULL 
seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide$SampleSite_code <- NULL 
seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide$WP <- NULL 
seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide$Seq_length <- NULL 
seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide$reads <- NULL 
seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide$seqtabNoC_WP2_12S_tax_col_sum2 <- NULL 


seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide<-seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, normalised_reads= normalised_reads,)

seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide <- spread(seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide,SampleSite_time,normalised_reads)
seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide[is.na(seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide)] <- 0
seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide <- data.frame(t(seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide))

names(seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide) <- as.matrix(seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide[1, ])
seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide <- seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide[-1, ]
seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide[] <- lapply(seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide, function(x) type.convert(as.character(x)))

#write.csv(seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide,"NORM_WP2C_Gwash_ASV_table_wide_filtered_replicates.csv")
write.csv(seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide,"NORM_WP2C_Gwash_ASV_table_wide_filtered.csv")


#PERMANOVA ON TAXON
NORM_WP2C_Gwash_ASV_table_wide_filtered <- read_csv("NORM_WP2C_Gwash_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2C_Gwash_ASV_table_wide_filtered
adonis_data_samples<-adonis_data$X1
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(SampleSite_time = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = X1)

lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","WP","Days","Catchment_flow_m3_s","Average_Soil_Temp_5cm_oC","Total_Rainfall","Dist_from_lake","pH","conductivity_uS_cm","gran_alk_ueq_L")],adonis_data, by="SampleSite_time")

adonis_data<-na.omit(adonis_data)
adonis_data_meta<- adonis_data[c(1:10)]
adonis_data<-adonis_data[,-c(1:10)]

adonis_data_meta$Days<-as.numeric(adonis_data_meta$Days)

adon.results_WP2<-adonis(adonis_data ~ adonis_data_meta$Days+
                           adonis_data_meta$Dist_from_lake+
                           adonis_data_meta$Average_Soil_Temp_5cm_oC+
                           adonis_data_meta$Catchment_flow_m3_s+
                           adonis_data_meta$gran_alk_ueq_L+
                           adonis_data_meta$Total_Rainfall+
                           adonis_data_meta$conductivity_uS_cm+
                           adonis_data_meta$pH,
                         method="bray",perm=999)

print(adon.results_WP2)

#PLOTTING BETA DIVERSITY NMDS
abund_table_12S_SampleSite_time <- read.csv("NORM_WP2C_Gwash_ASV_table_wide_filtered.csv", header = TRUE)
abund_table_12S <- abund_table_12S_SampleSite_time[,-1]
com <- abund_table_12S[,col(abund_table_12S)]
m_com <- as.matrix(com)
set.seed(123)
nmds = metaMDS(m_com, distance = "bray", try = 1000, trymax = 1000, k = 3)

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds))
#add columns to data frame 
data.scores$SampleSite_time = abund_table_12S_SampleSite_time[,1]

lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]
data.scores_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","SampleSite_code","Days","Dist_from_lake")],data.scores, by="SampleSite_time")


p<-ggplot(data.scores_meta_data, aes(x = NMDS1, y = NMDS2,colour = Dist_from_lake)) + 
  geom_point(size = 2, aes(colour = Dist_from_lake))+
  theme_bw()
ggplotly(p)


#ALPHA DIVERSITY METRICS
shannon_index<-diversity(seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide, index = "shannon")
shannon_index<-data.frame(shannon_index)
ID <- rownames(shannon_index)
shannon_index <- cbind(shannon_index,ID)

simpson_index<-diversity(seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide, index = "simpson")
simpson_index<-data.frame(simpson_index)
ID <- rownames(simpson_index)
simpson_index <- cbind(simpson_index,ID)

species_count <- rowSums(seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide >0)
species_count<-data.frame(species_count)
species_count$ID <- rownames(species_count)

alpha_12S_WP2C_Gwash<-merge(shannon_index,simpson_index,by="ID")
alpha_12S_WP2C_Gwash<-merge(alpha_12S_WP2C_Gwash,species_count,by="ID")

#write.csv(alpha_12S_WP2C_Gwash,"alpha_12S_WP2C_Gwash_replicates.csv")
write.csv(alpha_12S_WP2C_Gwash,"alpha_12S_WP2C_Gwash.csv")


#________________________________________________________

#########WP2C_Skan BLAST  TAXONOMY ANALYSIS##########

#________________________________________________________

setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/12S_BLAST")
seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean_filtered <- read_csv("WP2C_Skan_ASV_table_long_filtered_family.csv")

#aggregating replicates
seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean_filtered %>% 
  group_by(SampleSite_time,species) %>% 
 summarize(reads=mean(reads)) %>%
 rename(SampleSite_time=SampleSite_time,species=species, reads= reads)

#REMOVE CLEAR FISH THAT ARE LIEKLY CONSUMED BY HUAMNS AND NOT NATIVE
#seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean_filtered$species=="Platichthys_stellatus"),]

#CONVERTING FILTERED LONG FORMAT INTO WSampleSite_timeE FORMAT
seqtabNoC_WP2C_Skan_12S_tax_wide_filtered<-seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean_filtered
seqtabNoC_WP2C_Skan_12S_tax_wide_filtered$read_filter <- NULL 
seqtabNoC_WP2C_Skan_12S_tax_wide_filtered$Days <- NULL 
seqtabNoC_WP2C_Skan_12S_tax_wide_filtered$Date_sampled <- NULL 
seqtabNoC_WP2C_Skan_12S_tax_wide_filtered$SampleSite_code <- NULL 
seqtabNoC_WP2C_Skan_12S_tax_wide_filtered$WP <- NULL 
seqtabNoC_WP2C_Skan_12S_tax_wide_filtered$Seq_length <- NULL 

seqtabNoC_WP2C_Skan_12S_tax_wide_filtered<-seqtabNoC_WP2C_Skan_12S_tax_wide_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=sum(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads,)

seqtabNoC_WP2C_Skan_12S_tax_wide_filtered <- spread(seqtabNoC_WP2C_Skan_12S_tax_wide_filtered,SampleSite_time,reads)
seqtabNoC_WP2C_Skan_12S_tax_wide_filtered[is.na(seqtabNoC_WP2C_Skan_12S_tax_wide_filtered)] <- 0
seqtabNoC_WP2C_Skan_12S_tax_wide_filtered <- data.frame(t(seqtabNoC_WP2C_Skan_12S_tax_wide_filtered))

names(seqtabNoC_WP2C_Skan_12S_tax_wide_filtered) <- as.matrix(seqtabNoC_WP2C_Skan_12S_tax_wide_filtered[1, ])
seqtabNoC_WP2C_Skan_12S_tax_wide_filtered <- seqtabNoC_WP2C_Skan_12S_tax_wide_filtered[-1, ]
seqtabNoC_WP2C_Skan_12S_tax_wide_filtered[] <- lapply(seqtabNoC_WP2C_Skan_12S_tax_wide_filtered, function(x) type.convert(as.character(x)))

#write.csv(seqtabNoC_WP2C_Skan_12S_tax_wide_filtered,"WP2C_Skan_ASV_table_wide_filtered_replicates.csv")
write.csv(seqtabNoC_WP2C_Skan_12S_tax_wide_filtered,"WP2C_Skan_ASV_table_wide_filtered.csv")

#RELATIVE AMPLICON FREQUENCY ABUNDANCE
ggplot(seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean_filtered, aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))

#ABSOLUTE AMPLICON ABUNDANCE
#ggplot(seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean , aes(x = SampleSite_time, fill = species, y = reads)) + 
# geom_bar(stat = "identity", colour = "black")+
#theme(axis.text.x = element_text(angle = 90),legend.position = "none")
##scale_fill_brewer(palette="Paired")

##SPLITTING BY SPECIES
#Alosa_sapidissima<-seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean_filtered %>% filter(species == "Alosa_sapidissima")

Anguilla_anguilla_sum_reads<-Anguilla_anguilla %>% 
  group_by(SampleSite_time) %>% 
  summarize(reads=sum(reads)) %>%
  rename(SampleSite_time_sum =SampleSite_time,reads_sum= reads)

ggplot(data=Anguilla_anguilla_sum_reads, aes(x=SampleSite_time_sum, y=reads_sum))+
  geom_col()

#Venn diagram
#seqtabNoC_WP2C_Skan_12S_Venn<-seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean[!(seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean$reads==0),]
#hist(seqtabNoC_WP2C_Skan_12S_Venn$reads, breaks=1000,xlim=c(0,50000))

#E01_Venn<-seqtabNoC_WP2C_Skan_12S_Venn %>% filter(SampleSite_code == "E01")
#E01_Venn<-E01_Venn$species

#E02_Venn<-seqtabNoC_WP2C_Skan_12S_Venn %>% filter(SampleSite_code == "E02")
#E02_Venn<-E02_Venn$species

#E03_Venn<-seqtabNoC_WP2C_Skan_12S_Venn %>% filter(SampleSite_code == "E03")
#E03_Venn<-E03_Venn$species

#E04_Venn<-seqtabNoC_WP2C_Skan_12S_Venn %>% filter(SampleSite_code == "E04")
#E04_Venn<-E04_Venn$species

#E05_Venn<-seqtabNoC_WP2C_Skan_12S_Venn %>% filter(SampleSite_code == "E05")
#E05_Venn<-E05_Venn$species


#venn.diagram(
# x = list(E01_Venn, E02_Venn,E03_Venn,E04_Venn,E05_Venn),
#category.names = c("E01" , "E02","E03","E04","E05"),
#filename = 'WP2C_Skan_venn_diagramm.png',
#output=TRUE)


#NORMALISING READS

seqtabNoC_WP2C_Skan_12S_tax_col_sum2<-rowSums(seqtabNoC_WP2C_Skan_12S_tax_wide_filtered[c(1:22)])
seqtabNoC_WP2C_Skan_12S_tax_col_sum2<-data.frame(seqtabNoC_WP2C_Skan_12S_tax_col_sum2)
seqtabNoC_WP2C_Skan_12S_tax_col_sum2$SampleSite_time <- rownames(seqtabNoC_WP2C_Skan_12S_tax_col_sum2)

seqtabNoC_WP2C_Skan_12S_tax_normalised<- merge(seqtabNoC_WP2C_Skan_12S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2C_Skan_12S_tax_col_sum2")],seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean_filtered, by="SampleSite_time")

#GETTING THE SUM OF EACH SAMPLE FOR NORMALISING READS
seqtabNoC_WP2C_Skan_12S_tax_col_sum3<-colSums(seqtabNoC_WP2C_Skan_12S_tax_wide_filtered[c(1:22)])
seqtabNoC_WP2C_Skan_12S_tax_col_sum3<-data.frame(seqtabNoC_WP2C_Skan_12S_tax_col_sum3)

seqtabNoC_WP2C_Skan_12S_tax_col_sum3$species <- rownames(seqtabNoC_WP2C_Skan_12S_tax_col_sum3)

seqtabNoC_WP2C_Skan_12S_tax_normalised<- merge(seqtabNoC_WP2C_Skan_12S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2C_Skan_12S_tax_col_sum2")],seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean_filtered, by="SampleSite_time")
seqtabNoC_WP2C_Skan_12S_tax_top<- merge(seqtabNoC_WP2C_Skan_12S_tax_col_sum3[, c("species", "seqtabNoC_WP2C_Skan_12S_tax_col_sum3")],seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean_filtered, by="species")

seqtabNoC_WP2C_Skan_12S_tax_normalised <- transform(seqtabNoC_WP2C_Skan_12S_tax_normalised, normalised_reads = reads / seqtabNoC_WP2C_Skan_12S_tax_col_sum2)

ggplot(seqtabNoC_WP2C_Skan_12S_tax_normalised , aes(x = SampleSite_time, fill = species, y = normalised_reads)) + 
  geom_bar(stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))

#PLOTTING TOP Taxon
seqtabNoC_WP2C_Skan_12S_tax_top %>%
  mutate(species = fct_reorder(species, desc(seqtabNoC_WP2C_Skan_12S_tax_col_sum3))) %>%
  ggplot(aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))+scale_fill_brewer(palette="Paired")+ 
  theme(legend.text=element_text(size=10))+
  theme(legend.key.size = unit(0.3, 'cm'))

#write.csv(seqtabNoC_WP2C_Skan_12S_tax_normalised,"NORM_WP2C_Skan_ASV_table_LONG_filtered_replicates.csv")
write.csv(seqtabNoC_WP2C_Skan_12S_tax_normalised,"NORM_WP2C_Skan_ASV_table_LONG_filtered.csv")


#CONVERTING NORMALISED LONG FORMAT INTO WSampleSite_timeE FORMAT
seqtabNoC_WP2C_Skan_12S_tax_normalised_wide<-seqtabNoC_WP2C_Skan_12S_tax_normalised
seqtabNoC_WP2C_Skan_12S_tax_normalised_wide$read_filter <- NULL 
seqtabNoC_WP2C_Skan_12S_tax_normalised_wide$Days <- NULL 
seqtabNoC_WP2C_Skan_12S_tax_normalised_wide$Date_sampled <- NULL 
seqtabNoC_WP2C_Skan_12S_tax_normalised_wide$SampleSite_code <- NULL 
seqtabNoC_WP2C_Skan_12S_tax_normalised_wide$WP <- NULL 
seqtabNoC_WP2C_Skan_12S_tax_normalised_wide$Seq_length <- NULL 
seqtabNoC_WP2C_Skan_12S_tax_normalised_wide$reads <- NULL 
seqtabNoC_WP2C_Skan_12S_tax_normalised_wide$seqtabNoC_WP2_12S_tax_col_sum2 <- NULL 


seqtabNoC_WP2C_Skan_12S_tax_normalised_wide<-seqtabNoC_WP2C_Skan_12S_tax_normalised_wide %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, normalised_reads= normalised_reads,)

seqtabNoC_WP2C_Skan_12S_tax_normalised_wide <- spread(seqtabNoC_WP2C_Skan_12S_tax_normalised_wide,SampleSite_time,normalised_reads)
seqtabNoC_WP2C_Skan_12S_tax_normalised_wide[is.na(seqtabNoC_WP2C_Skan_12S_tax_normalised_wide)] <- 0
seqtabNoC_WP2C_Skan_12S_tax_normalised_wide <- data.frame(t(seqtabNoC_WP2C_Skan_12S_tax_normalised_wide))

names(seqtabNoC_WP2C_Skan_12S_tax_normalised_wide) <- as.matrix(seqtabNoC_WP2C_Skan_12S_tax_normalised_wide[1, ])
seqtabNoC_WP2C_Skan_12S_tax_normalised_wide <- seqtabNoC_WP2C_Skan_12S_tax_normalised_wide[-1, ]
seqtabNoC_WP2C_Skan_12S_tax_normalised_wide[] <- lapply(seqtabNoC_WP2C_Skan_12S_tax_normalised_wide, function(x) type.convert(as.character(x)))

#write.csv(seqtabNoC_WP2C_Skan_12S_tax_normalised_wide,"NORM_WP2C_Skan_ASV_table_wide_filtered_replicates.csv")
write.csv(seqtabNoC_WP2C_Skan_12S_tax_normalised_wide,"NORM_WP2C_Skan_ASV_table_wide_filtered.csv")

#PERMANOVA ON TAXON
NORM_WP2C_Skan_ASV_table_wide_filtered <- read_csv("NORM_WP2C_Skan_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2C_Skan_ASV_table_wide_filtered
adonis_data_samples<-adonis_data$X1
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(SampleSite_time = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = X1)

lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","WP","Days","Catchment_flow_m3_s","Average_Soil_Temp_5cm_oC","Total_Rainfall","Dist_from_lake","pH","conductivity_uS_cm","gran_alk_ueq_L")],adonis_data, by="SampleSite_time")

adonis_data<-na.omit(adonis_data)
adonis_data_meta<- adonis_data[c(1:10)]
adonis_data<-adonis_data[,-c(1:10)]

adonis_data_meta$Days<-as.numeric(adonis_data_meta$Days)

adon.results_WP2<-adonis(adonis_data ~ adonis_data_meta$Days+
                           adonis_data_meta$Dist_from_lake+
                           adonis_data_meta$Average_Soil_Temp_5cm_oC+
                           adonis_data_meta$Catchment_flow_m3_s+
                           adonis_data_meta$gran_alk_ueq_L+
                           adonis_data_meta$Total_Rainfall+
                           adonis_data_meta$conductivity_uS_cm+
                           adonis_data_meta$pH,
                         method="bray",perm=999)

print(adon.results_WP2)

#PLOTTING BETA DIVERSITY NMDS
abund_table_12S_SampleSite_time <- read.csv("NORM_WP2C_Skan_ASV_table_wide_filtered.csv", header = TRUE)
abund_table_12S <- abund_table_12S_SampleSite_time[,-1]
com <- abund_table_12S[,col(abund_table_12S)]
m_com <- as.matrix(com)
set.seed(123)
nmds = metaMDS(m_com, distance = "bray", try = 1000, trymax = 1000, k = 3)

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds))
#add columns to data frame 
data.scores$SampleSite_time = abund_table_12S_SampleSite_time[,1]

lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]
data.scores_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","SampleSite_code","Days","Dist_from_lake")],data.scores, by="SampleSite_time")


p<-ggplot(data.scores_meta_data, aes(x = NMDS1, y = NMDS2,colour = Dist_from_lake)) + 
  geom_point(size = 2, aes(colour = Dist_from_lake))+
  theme_bw()
ggplotly(p)


#ALPHA DIVERSITY METRICS
shannon_index<-diversity(seqtabNoC_WP2C_Skan_12S_tax_normalised_wide, index = "shannon")
shannon_index<-data.frame(shannon_index)
ID <- rownames(shannon_index)
shannon_index <- cbind(shannon_index,ID)

simpson_index<-diversity(seqtabNoC_WP2C_Skan_12S_tax_normalised_wide, index = "simpson")
simpson_index<-data.frame(simpson_index)
ID <- rownames(simpson_index)
simpson_index <- cbind(simpson_index,ID)

species_count <- rowSums(seqtabNoC_WP2C_Skan_12S_tax_normalised_wide >0)
species_count<-data.frame(species_count)
species_count$ID <- rownames(species_count)

alpha_12S_WP2C_Skan<-merge(shannon_index,simpson_index,by="ID")
alpha_12S_WP2C_Skan<-merge(alpha_12S_WP2C_Skan,species_count,by="ID")

#write.csv(alpha_12S_WP2C_Skan,"alpha_12S_WP2C_Skan_replicates.csv")
write.csv(alpha_12S_WP2C_Skan,"alpha_12S_WP2C_Skan.csv")


#________________________________________________________

########COMPILING WP2 CURATED NORMALISED OTU TABLES (APART FROM WP2B) FOR ALPHA AND BETA DIVERSITY #####

#________________________________________________________

alpha_12S_WP2A <- read_csv("alpha_12S_WP2A.csv")
alpha_12S_WP2C_Glatt <- read_csv("alpha_12S_WP2C_Glatt.csv")
alpha_12S_WP2C_Gwash <- read_csv("alpha_12S_WP2C_Gwash.csv")
alpha_12S_WP2C_Towy <- read_csv("alpha_12S_WP2C_Towy.csv")
#alpha_12S_WP2C_Skan <- read_csv("alpha_12S_WP2C_Skan.csv")

alpha_12S_WP2<-rbind(alpha_12S_WP2A,alpha_12S_WP2C_Glatt,alpha_12S_WP2C_Gwash,alpha_12S_WP2C_Towy)#,alpha_12S_WP2C_Skan)
alpha_12S_WP2<-alpha_12S_WP2[!(alpha_12S_WP2$shannon_index==0),]

lofresh_metadata_WP2_3 <- read_csv("metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3_summary = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]  

alpha_12S_WP2<-alpha_12S_WP2 %>% rename(SampleSite_time = ID)

alpha_12S_WP2 <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "Dist_from_lake")],alpha_12S_WP2, by="SampleSite_time")
alpha_12S_WP2 <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "River")],alpha_12S_WP2, by="SampleSite_time")
alpha_12S_WP2 <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "Date_sampled")],alpha_12S_WP2, by="SampleSite_time")
alpha_12S_WP2 <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "pH")],alpha_12S_WP2, by="SampleSite_time")
alpha_12S_WP2 <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "conductivity_uS_cm")],alpha_12S_WP2, by="SampleSite_time")

#Only keep the conwy summer samples that were coordinated with the Gwash, Glatt and Towy
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

ggplot(alpha_12S_WP2, aes(x=Dist_from_lake, y=species_count, color=River)) +
  geom_point() + 
  geom_smooth(method = "lm",fullrange = T)+theme_bw()+scale_color_manual(values=c("#56B4E9","#F0E442","#E69F00","#0072B2","#D55E00","#999999"))


hist(alpha_12S_WP2$shannon_index)

#ggplot(alpha_12S_WP2, aes(x=Dist_from_lake, y=shannon_index, color=River)) +
 #geom_point() + 
  #geom_smooth(method = "lm", formula = y ~ poly(x, 2))+theme_bw()+scale_color_manual(values=c("#56B4E9","#F0E442","#E69F00","#0072B2","#D55E00","#999999"))

#lm_12S_WP2 = lm(shannon_index~I(Dist_from_lake^2)*River, data = alpha_12S_WP2)
#lm_12S_WP2.1 = lm(shannon_index~Dist_from_lake*River, data = alpha_12S_WP2)
#AIC(lm_12S_WP2)
#AIC(lm_12S_WP2.1)
#anova(lm_12S_WP2,lm_12S_WP2.1)

lm_12S_WP2.1 = lm(species_count~Dist_from_lake*River, data = alpha_12S_WP2)
anova(lm_12S_WP2.1)


#BETA DIVERSITY - PERMANOVA ON TAXON
setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/12S_BLAST")
NORM_WP2A_ASV_table_LONG_filtered <- read_csv("NORM_WP2A_ASV_table_LONG_filtered.csv")
NORM_WP2C_Glatt_ASV_table_LONG_filtered <- read_csv("NORM_WP2C_Glatt_ASV_table_LONG_filtered.csv")
NORM_WP2C_Gwash_ASV_table_LONG_filtered <- read_csv("NORM_WP2C_Gwash_ASV_table_LONG_filtered.csv")
NORM_WP2C_Towy_ASV_table_LONG_filtered <- read_csv("NORM_WP2C_Towy_ASV_table_LONG_filtered.csv")
NORM_WP2C_Skan_ASV_table_LONG_filtered <- read_csv("NORM_WP2C_Skan_ASV_table_LONG_filtered.csv")

NORM_WP2A_ASV_table_LONG_filtered$Days<-NULL
NORM_WP2A_ASV_table_LONG_filtered$SampleSite_code<-NULL

NORM_WP2A_ASV_table_LONG_filtered[,3] <- NULL
NORM_WP2C_Glatt_ASV_table_LONG_filtered[,3] <- NULL
NORM_WP2C_Gwash_ASV_table_LONG_filtered[,3] <- NULL
NORM_WP2C_Towy_ASV_table_LONG_filtered[,3] <- NULL
NORM_WP2C_Skan_ASV_table_LONG_filtered[,3] <- NULL

NORM_WP2_ASV_table_LONG_filtered<-rbind(NORM_WP2A_ASV_table_LONG_filtered,NORM_WP2C_Glatt_ASV_table_LONG_filtered,
                                        NORM_WP2C_Gwash_ASV_table_LONG_filtered ,NORM_WP2C_Towy_ASV_table_LONG_filtered,NORM_WP2C_Skan_ASV_table_LONG_filtered)

#CONVERTING NORMALISED LONG FORMAT INTO WIDE FORMAT
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

adon.results_WP2<-adonis2(adonis_data ~ adonis_data_meta$River*adonis_data_meta$Dist_from_lake,
                         method="bray",perm=999)
print(adon.results_WP2)

#PLOTTING BETA DIVERSITY NMDS
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

gg <- merge(data.scores_meta_data,aggregate(cbind(NMDS1.x=NMDS1,NMDS2.y=NMDS2)~River,data.scores_meta_data,mean),by="River")
ggplot(gg, aes(x = NMDS1, y = NMDS2,colour = Dist_from_lake)) + scale_shape_manual(values=c(15,16,17,8,18))+
  geom_point(size = 4, aes(shape=River,colour = Dist_from_lake))+
  theme_bw()+  theme(axis.title.x=element_blank())+  theme(axis.title.y=element_blank())+ 
  geom_point(aes(x=NMDS1.x,y=NMDS2.y),size=3,colour="red")+
  geom_segment(aes(x=NMDS1.x, y=NMDS2.y, xend=NMDS1, yend=NMDS2))
```
