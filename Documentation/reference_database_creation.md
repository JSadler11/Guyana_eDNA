# Reference Database Curation & Classifier Building

These details were adapted from the Ecuador Bird Poop database created and shared by Dr. Sophie Picq and [Mike Allen](https://www.mikeallen-eco.com/news/2024/7/25/making-a-metabarcoding-reference-database-for-arthropods-part-1-downloading-sequences-1). Do not distribute outside of the lab without consulting Dr. Picq directly. Dr. Allen's blog is publicly available and I strongly recommend reading in full if this is new.

The general order of operations is as follows, pulled from [Devon O'Rourke's](https://forum.qiime2.org/t/building-a-coi-database-from-bold-references/16129) tutorial.

1. Collect database sequences and associated taxonomy information.
2. Reformat these data and import into QIIME
3. Filter these data in QIIME using options available via the RESCRIPt plugin
4. Evaluate the database using options in RESCRIPt
5. Create a naive Bayes classifier using RESCRIPt which can be used in downstream classification of your particular experiment.

Create a reference database that has sequences from BOLD, NCBI, and MIDORI on metazoa and plants of the area (qza).

Use QIIME to then make a tree of sample ASVs once they are retrieved from DADA2.
```
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs_single.qza \
  --output-dir mafft-fasttree-output

# turn into qzv to view
qiime tools export
  --input-path mafft-fastree-output/rooted_tree.qza \
  --output-path exported_tree

qiime tools export \
  --input-path mafft-fasttree-output/alignment.qza \
  --output-path exported_alignment

qiime tools export \
  --input-path mafft-fasttree-output/masked_alignment.qza \
  --output-path exported_masked_alignment
```
# 1. Download Reference Datasets
Download all reference database data from BOLD, EMBL, etc.

```
#Activate CRABS environment
#Downloading taxa from BOLD
export PATH="~/reference_database_creator:$PATH"

for i in "${TAXAnames[@]}";
do
  crabs
  --download-bold
  --taxon "$i"
  --output "bold_${i}.fasta"
done

cat *.fasta > bold_all_TAXA_MMDDYYYY.fasta

OR
# From Mike Allen

artnames=("XXXXXX" "XXXXXXXX" "XXXXXXXX")

for i in "${TAXAnames[@]}";
  do crabs db_download
  --source bold
  --database "$i"
  --output "bold_${i}.fasta"
  --keep_original no
done

#Download birds from BOLD (as an example)

birdnames=("Hellmayrea gularis" "Cinnycerthia unirufa" "Myioborus melanocephalus" "Myiothlypis nigrocristata" "Margarornis squamiger" "Ochthoeca diadema" "Myiothlypis coronata" "Henicorhina leucophrys" "Cinnycerthia olivascens" "Basileuterus tristriatus" "Myioborus miniatus" "Premnoplex brunnescens" "Pseudotriccus pelzelni" "Syndactyla subalaris" "Myiotriccus ornatus" "Myiobius villosus" "Xiphorhynchus erythropygius" "Myrmotherula schisticolor" "Glyphorynchus spirurus" "Myiothlypis chysogaster" "Dendrocincla fuliginosa" "Dendrocincla tyrannina")

for i in "${birdnames[@]}"; 
do crabs --download-bold --taxon "$i" --output "bold_${i}.fasta"
done

cat *.fasta > bold_all_birds_08252025.fasta
```

Downloading all arthropods from GenBank

```
crabs
  --download-ncbi
  --query '(Arthropoda [Organism] AND (cytochrome c oxidase subunit 1[Title] OR cytochrome c oxidase subunit       I[Title] OR cytochrome oxidase subunit 1[Title] OR cytochrome oxidase subunit I[Title] OR COX1[Title] OR CO1[Title] OR COI[Title]))'
  --output ~/crabs-downloads/NCBI/ncbi_PHYLUM_MMDDYYYY.fasta
  --email jsadler@wesleyan.edu
  --database nucleotide

# Ran on MSU HPCC as such:
# SBATCH --time=08:00:00
# SBATCH --nodes=1
# SBATCH --ntasks-per-node=1
# SBATCH --cpus-per-task=15
# SBATCH --mem-per-cpu=32G
# SBATCH --job-name ncbi_art

# write job information to SLURM output file
scontrol show job $SLURM_JOB_ID

#write resource usage to SLURM output file (uses a powertools command)
module load powertools
js -j $SLURM_JOB_ID

module purge
module load Miniforge3
conda activate py39
export PATH="~/reference_database_creator:$PATH"

crabs
  --download-ncbi
  --query '(txid##PHYLUM## [ORGN] AND (cytochrome c oxidase subunit 1[Title] OR cytochrome c oxidase subunit I[Title] OR cytochrome oxidase subunit 1[Title] OR cytochrome oxidase subunit I[Title] OR COX1[Title] OR CO1[Title] OR COI[Title]) NOT environmental sample[Title] NOT environmental samples[Title] NOT environmental[Title] NOT uncultured[Title] NOT unclassified[Title] NOT unidentified[Title] NOT unverified[Title] AND ("200"[SLEN] : "50000"[SLEN]))'
  --output ~/crabs-downloads/NCBI/ncbi_PHYLUM_MMDDYYYY.fasta
  --email jsadler@wesleyan.edu
  --database nucleotide
```
Example code from Sophie for birds:

```
crabs
  --download-ncbi
  --query '("Hellmayrea gularis" [ORGANISM] OR "Cinnycerthia unirufa"[ORGANISM] OR "Myioborus melanocephalus" [ORGANISM] OR "Myiothlypis nigrocristata" [ORGANISM] OR "Margarornis squamiger" [ORGANISM] OR "Ochthoeca diadema" [ORGANISM] OR "Myiothlypis coronata"[ORGANISM] OR "Henicorhina leucophrys" [ORGANISM] OR "Cinnycerthia olivascens" [ORGANISM] OR "Basileuterus tristriatus" [ORGANISM] OR "Myioborus miniatus" [ORGANISM] OR "Premnoplex brunnescens"[ORGANISM] OR "Pseudotriccus pelzelni"[ORGANISM] OR "Syndactyla subalaris" [ORGANISM] OR "Myiotriccus ornatus"[ORGANISM] OR "Myiobius villosus"[ORGANISM] OR "Xiphorhynchus erythropygius"[ORGANISM] OR "Myrmotherula schisticolor"[ORGANISM] OR "Glyphorynchus spirurus"[ORGANISM] OR "Myiothlypis chysogaster" [ORGANISM] OR "Dendrocincla fuliginosa"[ORGANISM] OR "Dendrocincla tyrannina" AND (cytochrome c oxidase subunit 1[Title] OR cytochrome c oxidase subunit I[Title] OR cytochrome oxidase subunit 1[Title] OR cytochrome oxidase subunit I[Title] OR COX1[Title] OR CO1[Title] OR COI[Title]) NOT environmental sample[Title] NOT environmental samples[Title] NOT environmental[Title] NOT uncultured[Title] NOT unclassified[Title] NOT unidentified[Title] NOT unverified[Title] AND ("200"[SLEN] : "50000"[SLEN]))'
  --output ~/crabs-downloads/ncbi_bird/ncbi_birds_08252025.fasta
  --email spicq@fieldmuseum.org
  --database nucleotide
```
Downloading from EMBL

```
# First, make a list of species suffixes

suf=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19")

# loop through the suffixes to download EMBL in manageable batches
for j in "${suf[@]}"; 
  do crabs
    --download-embl
    --taxon "TAX_${j}*"
    --output "embl_TAX_${j}.fasta"
done

cat *.fasta > embl_all_inv_08222025.fasta
```
To download CO1 MIDORI

```
# to download the MIDORI CO1 fasta file
wget https://www.reference-midori.info/download/Databases/GenBank267_2025_06-19/RAW/uniq/MIDORI_UNIQ_NUC_GB267_CO1_RAW.fasta.gz

# unzip it
gunzip *.fasta.gz

# convert periods into spaces so that CRABS can read the accession numbers
sed 's/\./ /g' MIDORI2_UNIQ_NUC_GB267_CO1_RAW.fasta > midori2.fix.fasta"

# import sequences into crabs
crabs db_import
  --input midori2.fix.fasta \
  --output midori2.fix.crb.fasta \
  --seq-header accession \
  --delim ' '

crabs
  --download-midori
  --output coi_total_267.fasta
  --gb-number 267_2025-06-19
  --gene CO1
  --gb-type total
```
Download taxonomy
```
crabs
  --download-taxonomy
  --output crabtax
```
# 2. Import databases into CRABS

### NCBI

```
crabs
  --import
  --import-format ncbi
  --input ~/crabs-downloads/NCBI_inv/ncbi_PHYLUM_MMDDYYYY.fasta
  --names names.dmp
  --nodes nodes.dmp
  --acc2tax nucl_gb.accession2taxid
  --output ~/EXAMPLEecuador_poop_novaseq/ref_db/NCBI_art_CRABS_08212025.txt
  --ranks 'superkingdom;phylum;class;order;family;genus;species'
```
### MIDORI

```
crabs
  --import
  --import-format midori
  --input ~/crabs-downloads/MIDORI/coi_total_267.fasta
  --names names.dmp
  --nodes nodes.dmp
  --acc2tax nucl_gb.accession2taxid
  --output ~/crabs-downloads/MIDORI_art_CRABS_267.txt
  --ranks 'superkingdom;phylum;class;order;family;genus;species'
```

### BOLD

```
crabs --import --import-format bold --input ~/crabs-downloads/BOLD/bold_all_art_08222025.fasta --names names.dmp --nodes nodes.dmp --acc2tax nucl_gb.accession2taxid --output ~/crabs-downloads/BOLD_art_CRABS_08212025.txt --ranks 'superkingdom;phylum;class;order;family;genus;species'
```

### EMBL

```
crabs --import --import-format embl --input ~/crabs-downloads/EMBL/embl_all_inv_08222025.fasta --names names.dmp --nodes nodes.dmp --acc2tax nucl_gb.accession2taxid --output ~/crabs-downloads/EMBL_art_CRABS_08212025.txt --ranks 'superkingdom;phylum;class;order;family;genus;species'
```
## Merge Sequences

```
crabs --merge --input 'NCBI_art_CRABS_08212025.txt;NCBI_birds_CRABS_08252025.txt;MIDORI_artbirds_CRABS_267.txt;BOLD_art_CRABS_08212025.txt;BOLD_birds_CRABS_08252025.txt;EMBL_art_CRABS_08212025.txt' --uniq --output ArtBirds_ALL.txt
```

### Extract regions through _in silico_ PCR

```
ulimit -n 65536

crabs --in-silico-pcr --input ArtBirds_ALL.txt --output ArtBirds_ALL_V5.txt --forward GGTCAACAAATCATAAAGATATTGG --reverse GGWACTAATCAATTTCCAAATCC --mismatch 4.5
```

### Retrieve amplicons without primer-binding regions
This can take some time. 
```
crabs --pairwise-global-alignment --input ArtBirds_ALL_c.txt --amplicons ArtBirds_ALL_V5.txt --output extra_V5_ArtBirds_ALL.txt --forward GGTCAACAAATCATAAAGATATTGG --reverse GGWACTAATCAATTTCCAAATCC  --percent-identity 0.60 --coverage 95
```

### Dereplicate the database
```
crabs --dereplicate --input extra_V5_ArtBirds_ALL.txt --output extra_V5_ArtBirds_ALL_derep.txt --dereplication-method 'unique_species'
```

### Filter (N>10)
```
crabs --filter --input extra_V5_ArtBirds_ALL_derep.txt --output extra_V5_ArtBirds_ALL_derep_filtered.txt --maximum-n 10
```
### Export into QIIME format
```
crabs
  --export
  --input extra_V5_ArtBirds_ALL_derep.txt
  --output QIIME2_Ecuador_tax.txt
  --export-format 'qiime-text'

crabs
  --export
  --input extra_V5_ArtBirds_ALL_derep.txt
  --output  QIIME2_Ecuador_seq.txt
  --export-format 'qiime-fasta'
```

### Import into QIIME
```
conda activate qiime2-amplicon-2024.10
## taxonomy file import
 qiime tools import \
   --type 'FeatureData[Taxonomy]' \
   --input-format HeaderlessTSVTaxonomyFormat \
   --input-path QIIME2_Ecuador_tax.txt \
   --output-path QIIME2_Ecuador_tax.qza

## fasta file import
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path QIIME2_Ecuador_seq.txt \
  --output-path QIIME2_Ecuador_seq.qza
```

# Building a classifier using RESCRIPT

This example, from [Devon O'Rourke](https://forum.qiime2.org/t/building-a-coi-database-from-bold-references/16129) uses the CO1 gene marker.

Trimming BOLD sequences involves a few steps before dereplication.
- Ambiguous nucleotide content (5 or more N's)
- Long homopolymer runs (12 or more)
- Very short (<250 bp) sequences
- Very long (>1600 bp) sequences





