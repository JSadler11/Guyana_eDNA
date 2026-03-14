
## Import into Qiime

```
conda activate qiime2-amplicon-2024.10

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path /Volumes/ROSALIND/eDNA_Pilots/Pilot/_manifest.tsv \
--output-path demux.qza \
--input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
--i-data demux.qza \
--p-n 200000 \
--o-visualization demux.qzv
```

## Trimming

Worth running FastQC before and after trimming and denoising at least for one set of samples to monitor behavior.

Trimming primers according to Sophie's pipeline. The primers I have listed are:

Forward: 12SV5F=ACTGGGATTAGATACCCC

Reverse: 12SV5R=TAGAACAGGCTCCTCTAG
Detecting the reverse primers (in RC) will get rid of the adapters too as they are before them 
(in the 5’ to 3’ orientation). Sophie put a fairly high error rate because there were errors 
while inspecting the fastq pre processing files:


Trim amplicon primers in two steps only the reverse for now

```
qiime cutadapt trim-paired \
--i-demultiplexed-sequences demux.qza \
--p-cores 6 \
--p-adapter-f CTAGAGGAGCCTGTTCTA \
--p-adapter-r GGGGTATCTAATCCCAGT \
--p-match-read-wildcards TRUE \
--p-error-rate 0.25 \
--o-trimmed-sequences eDNA_SC_12SV5.qza \
--verbose \
&

qiime demux summarize \
--i-data eDNA_SC_12SV5.qza \
--o-visualization eDNA_SC_12SV5.qzv

qiime cutadapt trim-paired \
--i-demultiplexed-sequences eDNA_SC_12SV5.qza \
--p-cores 6 \
--p-front-f ACTGGGATTAGATACCCC \
--p-front-r TAGAACAGGCTCCTCTAG \
--p-match-read-wildcards TRUE \
--p-discard-untrimmed TRUE \
--p-match-adapter-wildcards TRUE \
--o-trimmed-sequences eDNA_SC_12SV5.qza \
--verbose \
&

qiime demux summarize \
--i-data eDNA_SC_cutadapt_12SV5.qza \
--o-visualization eDNA_SC_cutadapt_12SV5.qzv
```
## Denoising
```
qiime dada2 denoise-paired \
--p-n-threads 20 \
--i-demultiplexed-seqs eDNA_SC_12SV5.qza \
--p-trim-left-f 0 \
--p-trim-left-r 0 \
--p-trunc-len-f 95 \
--p-trunc-len-r 95 \
--output-dir DADA2_eDNA_SC_12SV5_len95 \
--o-representative-sequences eDNA-SC_12SV5-repseqs.qza \
--verbose \
&

# Denoising stats
qiime metadata tabulate \
--m-input-file denoising_stats.qza \
--o-visualization denoising_stats.qzv

# Representative sequences
qiime feature-table tabulate-seqs \
--i-data rep_seqs_SC_nano.qza \
--o-visualization rep_seqs_SC_12SV5.qzv

# Feature table
qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv
```
The feature-table and tabulate-seqs commands will prove a mapping of feature IDs to sequences, and provide links for BLAST'ing each sequence against the NCBI nt database. 

You can also include the feature freq. information by passing it as metadata, but in this case we are looking at _feature_ metadata, and not _sample_ metadata. 

Next, you can filter your feature data. First, remove any ASVs that only occur once. Then, use the resulting feature table to filter the sample sequences to only the ones present in the new table. 

```
qiime feature-table filter-features \
```

## QZA Conversions in QIIME

We will need to convert the qza repseq and table output files to fasta (a type of text file) and tsv files for work in R. We will run the first steps in QIIME2 using the export feature, and also with biomformat in R. You can technically do this all in R, but I like using QIIME because there are fewer points to screw everything up. I'm looking at you row alignment.

```
qiime tools export \
--input-path repseqs.qza
--output-path repseqs
```
If you're running the repseqs file, you'll be able to import into R and convert (coming shortly below). Repeat with your table file. You should get a table.biom file. Move that down to the next step for biom conversion. 

```
biom convert \
-i table.biom \
-o table.tsv \
--to-tsv
```
For running biomformat in R, use the following to convert from biom to tsv. This is optional to the above conversion. 

```
library(biomformat)
biom_obj<-read_biom("table.biom")
asv_table<-as(biom_data(biom_obj), "matrix")
```

## FastQC

Before we go any further, here is how I ran my FastQC analyses. In the command line, run:
```
brew install fastqc
```
Next, cd'd in the directory with all of your .fastq.gz seq data run:
```
fastqc -t 4 -o /Volumes/ROSALIND/eDNA_Pilots/WWorth/QualityControlReports *.fastq.gz
```
## BLASTing ASVs

```
blastn -query /Volumes/ROSALIND/sequences.fasta
-task megablast
-db nt_euk
-out /Volumes/ROSALIND/blast_suramacreek.txt
-max_target_seqs 200
-perc_identity 90
-qcov_hsp_perc 95
-num_threads 8
-outfmt "6 delim=, std qlen slen staxids sscinames scomnames sfamilies sorders sskingdoms"
```
Now to filter our output ASV files we want to run the following in R:

```
library(tidyr)
library(dplyr)
```

For blast filtering we want to start with 98% identity and e-values less than 0.001
Limit results to the top five best-scoring hits per query (max_target_seqs 5)

Our outputs don't have the above though, so we will do that in R'

Read in BLAST results
```
blast_12S <- blast_suramacreek_11_14
blast_18S <- blast_suramacreek_18S_11_14
```

Filter by e-value and percent identity
```
filtered_12S <- blast_12S %>%
  filter(blast_12S$evalue <= 1e-3) %>%      # Keep only e-values <= 1e-3
  filter(blast_12S$pident >= 98)                # Keep only matches with >= 95% identity

filtered_18S <- blast_18S %>%
  filter(blast_18S$evalue <= 1e-3) %>%      # Keep only e-values <= 1e-3
  filter(blast_18S$pident >= 98)   

```
View the filtered results
```
head(filtered_12S)
head(filtered_18S)

filtered_12S <- filtered_12S %>%
  filter(sscinames != "Homo sapiens")

sum(filtered_12S$sscinames == "Homo sapiens")

unique(filtered_12S$sscinames)

species_counts_12S <- filtered_12S %>%
  count(sscinames, sort = TRUE) %>%
  rename(Species = sscinames, Count = n)

species_counts_18S <- filtered_18S %>%
  count(sscinames, sort = TRUE) %>%
  rename(Species = sscinames, Count = n)

```

Optional: Get the best hit per query
```
best_hits_12S <- filtered_12S %>%
  group_by(qseqid) %>%
  slice_min(evalue, n = 5, with_ties = TRUE)
  ungroup()

best_hits_18S <- filtered_18S %>%
  group_by(qseqid) %>%
  slice_min(evalue, n = 5, with_ties = TRUE) %>%       # Get lowest e-value per query
  ungroup()

#Other
## Create manifest

Running in Docker, first we want to mount the external drive: 

```
docker run -it --rm \
  -v "/Volumes/ROSALIND:/Volumes/ROSALIND" \
  quay.io/qiime2/amplicon:2024.10
```

```
LIB="project_12SV5"

pwd > "$LIB"_pwd.tmptxt
find . -name "*.gz" | sort -u | cut -d '/' -f 2 > "$LIB"_filenames.tmptxt
cut -f 1 -d "_" "$LIB"_filenames.tmptxt | sort -u > "$LIB"_col1.tmptxt
wc -l "$LIB"_col1.tmptxt | awk '{print $1}' > "$LIB"_lines.tmptxt
yes $(cat "$LIB"_pwd.tmptxt) | head -n $(cat "$LIB"_lines.tmptxt) > "$LIB"_dirpath.tmptxt

for i in $(cat "$LIB"_col1.tmptxt); do
  find . -name "${i}_*.1.fastq.gz" | head -n 1 | sed 's|^\./||'
done > "$LIB"_forward_files.tmptxt
paste -d "/" "$LIB"_dirpath.tmptxt "$LIB"_forward_files.tmptxt > "$LIB"_col2.tmptxt

for i in $(cat "$LIB"_col1.tmptxt); do
  find . -name "${i}_*.2.fastq.gz" | head -n 1 | sed 's|^\./||'
done > "$LIB"_reverse_files.tmptxt
paste -d "/" "$LIB"_dirpath.tmptxt "$LIB"_reverse_files.tmptxt > "$LIB"_col3.tmptxt

paste -d "," "$LIB"_col1.tmptxt "$LIB"_col2.tmptxt "$LIB"_col3.tmptxt > "$LIB"_manifest.tmptxt
echo 'sample-id,forward-absolute-filepath,reverse-absolute-filepath' | cat - "$LIB"_manifest.tmptxt > ../$(echo "$LIB").manifest.file
rm *.tmptxt

```

To view and convert the manifest file into a .tsv, add ".csv" extension to the file name and run the following to open in Excel. 

```
# Open with Excel
open -a "Microsoft Excel" /Volumes/ROSALIND/eDNA_Pilots/SCreek/SC_trimmed/12SV5.manifest.file.csv
```

  ```
