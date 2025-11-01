## Importing into QIIME

```
conda activate qiime2-amplicon-2024.10

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path /Volumes/ROSALIND/eDNA_Pilots/Pilot/_manifest.tsv \
--output-path demux.qza \
--input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
--i-data demux.qza \
--o-visualization demux.qzv
```
## Trimming
Trimming primers according to Sophie's pipeline. The primers I have listed are: 

  Forward:
  mlCO1intF=GGWACWGGWTGAACWGTWTAYCCYCC
  18SV7&8F=GYGGTGCATGGCCGTTSKTRGTT
  MiFishUF=GTCGGTAAAACTCGTGCCAGC
  ITS-S2F=ATGCGATACTTGGTGTGAAT
  16SV3&4F=GTGYCAGCMGCCGCGGTAA
  12SV5F=ACTGGGATTAGATACCCC

  Reverse:
  jgHCO2198=TAIACYTCIGGRTGICCRAARAAYCA
  18SV7&8R=GTGTGYACAAAGGBCAGGGAC
  MiFishUR=CATAGTGGGGTATCTAATCCCAGTTTG
  ITS-S3R=GACGCTTCTCCAGACTACAAT
  16SV3&4R=GGACTACNVGGGTWTCTAAT
  12SV5R=TAGAACAGGCTCCTCTAG

Detecting the reverse primers (in RC) will get rid of the adapters too 
as they are before them (in the 5’ to 3’ orientation). Sophie put a fairly 
high error rate because there were errors while inspecting the fastq 
pre processing files:
```

#Trim amplicon primers in two steps only the reverse for now
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
```

```

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
--verbose \
&

# Denoising stats
qiime metadata tabulate \
--m-input-file denoising_stats.qza \
--o-visualization denoising_stats.qzv

# Representative sequences
qiime feature-table tabulate-seqs \
--i-data rep_seqs_SC_nano.qza \
--o-visualization rep_seqs_SC_nano.qzv

# Feature table
qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv
```
