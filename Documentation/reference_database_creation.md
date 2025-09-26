# Reference Database Curation

These details were adapted from the Ecuador Bird Poop database created and shared by Dr. Sophie Picq and [Mike Allen](https://www.mikeallen-eco.com/news/2024/7/25/making-a-metabarcoding-reference-database-for-arthropods-part-1-downloading-sequences-1). Do not distribute outside of the lab without consulting Dr. Picq directly. Dr. Allen's blog is publicly available and I strongly recommend reading in full if this is new.

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
  do crabs --download-embl --taxon "INV_${j}*" --output "embl_inv_${j}.fasta"
done

cat *.fasta > embl_all_inv_08222025.fasta
```
To download CO1 MIDORI

```
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


