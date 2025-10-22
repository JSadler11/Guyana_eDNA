# QIIME2 & CRABS Workflow

These procedures are an adapted combination of the work done by Devon O'Rourke, in both [NH bat guano diversity analyses](https://github.com/devonorourke/nhguano/blob/master/docs/sequence_processing.md), and his tutorial on building out a [BOLD Arthropod Bayesion Classifier](https://forum.qiime2.org/t/building-a-coi-database-from-bold-references/16129). Additional details provided from the [QIIME2 Moving Pictures](https://amplicon-docs.qiime2.org/en/latest/tutorials/moving-pictures.html), and [Gut-to-soil Axis](https://amplicon-docs.qiime2.org/en/latest/tutorials/gut-to-soil.html) tutorials.

## FastQC

Before we go any further, here is how I ran my FastQC analyses. In the command line, run:
```
brew install fastqc
```
Next, cd'd in the directory with all of your .fastq.gz seq data run:
```
fastqc -t 4 -o /Volumes/ROSALIND/eDNA_Pilots/WWorth/QualityControlReports *.fastq.gz
```
## Back to CutAdapt

First we want to demultiplex our reads. These examples are using paired-end sequence data. Each read will be sorted according to primers. 

For command line history: cat ~/.zsh_history

$RAWDIR points to the libraries of stored sequence data.

```
RAWDIR=/Volumes/ROSALIND/eDNA_Pilots/WWorth/raw/all_reads/    
```

$OUTDIR points to an output directory where trimmed reads are stored.
```
OUTDIR=/Volumes/ROSALIND/eDNA_Pilots/WWorth/WW_trimmed
```
Demultiplexing all of the sequence data from primers using CutAdapt.
Full primer library is listed here, copy and paste only what you need for runs.
```
FORWARD_PRIMERS=("mlCO1intF=GGWACWGGWTGAACWGTWTAYCCYCC"
                "18SV7&8F=GYGGTGCATGGCCGTTSKTRGTT"
                "MiFishUF=GTCGGTAAAACTCGTGCCAGC"
                "ITS-S2F=ATGCGATACTTGGTGTGAAT"
                "16SV3&4F=GTGYCAGCMGCCGCGGTAA"
                "12SV5F=ACTGGGATTAGATACCCC")



REVERSE_PRIMERS=("jgHCO2198=TAIACYTCIGGRTGICCRAARAAYCA"
                "18SV7&8R=GTGTGYACAAAGGBCAGGGAC"
                "MiFishUR=CATAGTGGGGTATCTAATCCCAGTTTG"
                "ITS-S3R=GACGCTTCTCCAGACTACAAT"
                "16SV3&4R=GGACTACNVGGGTWTCTAAT"
                "12SV5R=TAGAACAGGCTCCTCTAG")

MARKER_NAMES=("CO1" "18S" "MiFishU" "ITS2" "16S" "12SV5")
```
To check to see what markers are logged in your terminal, run this echo command:
```
echo "Markers: ${MARKER_NAMES[@]}"
echo "Forward: ${FORWARD_PRIMERS[@]}"
echo "Reverse: ${REVERSE_PRIMERS[@]}"
```
Full script is below. We are using a zsh index, so don't use "[$i-1]" when indexing, use "[$i]". Additionally, zsh indexes will allow us to only use i when assigning array indices. However, to run the same code on bash we will need to use idx.
```
FORWARD_PRIMERS=("mlCO1intF=GGWACWGGWTGAACWGTWTAYCCYCC"
                "18SV7&8F=GYGGTGCATGGCCGTTSKTRGTT"
                "MiFishUF=GTCGGTAAAACTCGTGCCAGC"
                "ITS-S2F=ATGCGATACTTGGTGTGAAT"
                "16SV3&4F=GTGYCAGCMGCCGCGGTAA"
                "12SV5F=ACTGGGATTAGATACCCC")

REVERSE_PRIMERS=("jgHCO2198=TAIACYTCIGGRTGICCRAARAAYCA"
                "18SV7&8R=GTGTGYACAAAGGBCAGGGAC"
                "MiFishUR=CATAGTGGGGTATCTAATCCCAGTTTG"
                "ITS-S3R=GACGCTTCTCCAGACTACAAT"
                "16SV3&4R=GGACTACNVGGGTWTCTAAT"
                "12SV5R=TAGAACAGGCTCCTCTAG")

MARKER_NAMES=("CO1" "18S" "MiFishU" "ITS2" "16S" "12SV5")


## Make individual folders for Primers

for MARKER in "${MARKER_NAMES[@]}"; do
    mkdir -p "$OUTDIR/$MARKER"
done

## Now separate R1 from R2 and run.
## Replace N with number of primers you are running, (EVEN if you are only running one primer, you'll put "for i in {1..1}

for SAMPLE in $(ls $RAWDIR | grep "_R1.fastq.gz" | sed 's/_R1.fastq.gz//' | sort -u); do
  echo "Processing sample: $SAMPLE"

for i in {1..N}; do
  MARKER=${MARKER_NAMES[$i]}
  FWD=$(echo "${FORWARD_PRIMERS[$i]}" | cut -d'=' -f2)
  REV=$(echo "${REVERSE_PRIMERS[$i]}" | cut -d'=' -f2)

  echo "$MARKER : $FWD / $REV"

  cutadapt \
    -g "$FWD" \
    -G "$REV" \
    --cores=8 \
    --discard-untrimmed \
    -o "$OUTDIR/$MARKER/${SAMPLE}_${MARKER}.1.fastq.gz" \
    -p "$OUTDIR/$MARKER/${SAMPLE}_${MARKER}.2.fastq.gz" \
    "$RAWDIR/${SAMPLE}_R1.fastq.gz" \
    "$RAWDIR/${SAMPLE}_R2.fastq.gz" \
    > "$OUTDIR/$MARKER/${SAMPLE}_${MARKER}_cutadapt_report.txt" 2>&1
  done
done
```
Note that while we have trimmed the primers and sorted the reads according to their primers, we have not trimmed for quality. That will be done after importing to QIIME2 for Quality Control checks. 

To import into QIIME2, we need to look into making a manifest file [O'Rourke](https://github.com/devonorourke/nhguano/blob/master/docs/sequence_processing.md).

To activate QIIME: 

```
conda activate qiime2-amplicon-2024.10
```

Next, you want to use a QIIME import function to make a .qza artifact containing all the unjoined paired end data. $QIIMEDIR refers to the path to the output path for the resulting QIIME artifact outputs.

```
## create manifest file according to O'Rourke (filename:absolute-path:direction)

pwd > "$LIB"_pwd.tmptxt
## This collects the current working directory and saves it to a temporary file named {LIB}_pwd.tmptxt, where $LIB is a variable.

find . -name "*.gz" | sort -u | cut -d '/' -f 2 > "$LIB"_filenames.tmptxt
##Finds all .gz files in current directory and subdirectories
##sorts the list and removes all duplicates
##Splits paths by / and takes the second field (aka, extracts the filenames from the ./filename)
##Saves results

cut -f 1 -d "_" "$LIB"_filenames.tmptxt > "$LIB"_col1.tmptxt
## Extracts sample ID by taking everything before the first underscore in the filename

wc -l "$LIB"_col1.tmptxt | cut -f 1 -d ' ' > "$LIB"_lines.tmptxt
## Counts the number of lines in each file, extracts just the number, and saves count to LIB

seq $(echo $(cat "$LIB"_lines.tmptxt)) | xargs -Iz echo $(cat "$LIB"_pwd.tmptxt) > "$LIB"_dirpath.tmptxt
##Generates a sequence from 1 to N number of files. For each number, it repeats the directory path N times, and creates a file with the same directory repeated for each file.

paste "$LIB"_dirpath.tmptxt "$LIB"_filenames.tmptxt -d "/" > "$LIB"_col2.tmptxt
## Combines the directory path and the filename with / as a separator, creating full absolute file paths.

for i in $(cat "$LIB"_filenames.tmptxt); do
if [[ $i == *.1.fastq.gz ]];
then
  echo "forward"
else
  echo "reverse"
fi;
done > "$LIB"_col3.tmptxt
paste "$LIB"_col1.tmptxt "$LIB"_col2.tmptxt "$LIB"_col3.tmptxt -d "," > "$LIB"_manifest.tmptxt
echo 'sample-id,absolute-filepath,direction' | cat - "$LIB"_manifest.tmptxt > ../$(echo "$LIB").manifest.file
rm *.tmptxt
## Filter forward (.1.) and reverse (.2.) reads and combine all three columns into a CSV format while removing all temporary files

```
All together:

```
##For MiFishU

LIB="MiFishU"

pwd > "$LIB"_pwd.tmptxt
find . -name "*.gz" | sort -u | cut -d '/' -f 2 > "$LIB"_filenames.tmptxt
cut -f 1 -d "_" "$LIB"_filenames.tmptxt > "$LIB"_col1.tmptxt
wc -l "$LIB"_col1.tmptxt | awk '{print $1}' > "$LIB"_lines.tmptxt
yes $(cat "$LIB"_pwd.tmptxt) | head -n $(cat "$LIB"_lines.tmptxt) > "$LIB"_dirpath.tmptxt
paste -d "/" "$LIB"_dirpath.tmptxt "$LIB"_filenames.tmptxt > "$LIB"_col2.tmptxt
for i in $(cat "$LIB"_filenames.tmptxt); do
if [[ $i == *.1.fastq.gz ]];
then
  echo "forward"
else
  echo "reverse"
fi;
done > "$LIB"_col3.tmptxt
paste -d "," "$LIB"_col1.tmptxt "$LIB"_col2.tmptxt "$LIB"_col3.tmptxt > "$LIB"_manifest.tmptxt
echo 'sample-id,absolute-filepath,direction' | cat - "$LIB"_manifest.tmptxt > ../$(echo "$LIB").manifest.file
rm *.tmptxt

```

However, the O'Rourke manifest may not work for paired end datasets (Not sure why, as his IS paired end). So, we reformatted and wrote this up. 

```
LIB="MiFishU"

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

Import sequence data [O'Rourke](https://github.com/devonorourke/nhguano/blob/master/docs/sequence_processing.md)

```
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path /Volumes/ROSALIND/eDNA_Pilots/WWorth/ww_manifest.tsv \
--output-path demux.qza \
--input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
--i-data demux.qza \
--p-n 50000 \\
--o-visualization demux-sumry.qzv
```
qiime demux paired end reads with barcodes and metadata file. Two steps are listed below, first for primer filtering, then for trimming.
```
qiime cutadapt demux-paired \
  --i-seqs demux.qza \
  --m-forward-barcodes-file primer-metadata.tsv \
  --m-forward-barcodes-column forward-primer-sequence \
  --m-reverse-barcodes-file primer-metadata.tsv \
  --m-reverse-barcodes-column reverse-primer-sequence \
  --o-per-sample-sequences demuxed-by-primer.qza \
  --o-untrimmed-sequences untrimmed.qza

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux.qza \
  --p-front-f GGWACWGGWTGAACWGTWTAYCCYCC \ # CO1
  --p-front-f GTCGGTAAAACTCGTGCCAGC \ # MiFishU
  --p-front-f GTGYCAGCMGCCGCGGTAA \ # 16S
  --p-front-f ATGCGATACTTGGTGTGAAT \ # ITS2
  --p-front-f GYGGTGCATGGCCGTTSKTRGTT \ #18S
  --p-front-r TAIACYTCIGGRTGICCRAARAAYCA \ # CO1
  --p-front-r CATAGTGGGGTATCTAATCCCAGTTTG \ # MiFishU
  --p-front-r GGACTACNVGGGTWTCTAAT \ # 16S
  --p-front-r GACGCTTCTCCAGACTACAAT \ # ITS2
  --p-front-r GTGTGYACAAAGGBCAGGGAC \ # 18S
  --p-discard-untrimmed \
  --p-cores 8 \
  --o-trimmed-sequences demux-trimmed.qza
```
## DADA2 Denoise Paired-End

```
## generate repseqs, the repseqs and MiFish-table are also your ASVs
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs MiFish-demux.qza \
  --p-n-threads 8 \
  --p-trim-left-f 0 \
  --p-trunc-len-f 275 \
  --p-trim-left-r 0 \
  --p-trunc-len-r 275 \
  --o-denoising-stats MiFish-denoising-stats.qza \
  --o-table MiFish-table.qza \ 
  --o-representative-sequences MiFish-repSeqs.qza

## generate summary visualization for metadata
qiime metadata tabulate \
  --m-input-file "$LIB".denoisingStats.qza
  --o-visualization "$LIB".denoisingStats.qzv  

## Summarize feature table and feature data

qiime feature-table summarize-plus \
  --i-table MiFish-table.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-summary asv-table.qzv \
  --o-sample-frequencies sample-frequencies.qza \
  --o-feature-frequencies asv-frequencies.qza

qiime feature-table tabulate-seqs \
  --i-data asv-seqs.qza \
  --m-metadata-file asv-frequencies.qza \
  --o-visualization asv-seqs.qzv

```

The feature-table and tabulate-seqs commands will prove a mapping of feature IDs to sequences, and provide links for BLAST'ing each sequence against the NCBI nt database. 

You can also include the feature freq. information by passing it as metadata, but in this case we are looking at _feature_ metadata, and not _sample_ metadata. 

Next, you can filter your feature data. First, remove any ASVs that only occur once. Then, use the resulting feature table to filter the sample sequences to only the ones present in the new table. 

```
qiime feature-table filter-features \
  --i-table asv-table.qza \
  --p-min-samples 2 \
  --o-filtered-table asv-table-ms2.qza

qiime feature-table filter-seqs \
  --i-data asv-seqs.qza \
  --i-table asv-table-ms2.qza \
  --o-filtered-data asv-seqs-ms2.qza
```
Now, summarize your filtered data.

```
qiime feature-table summarize-plus \
  --i-table asv-table-ms2.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-summary asv-table-ms2.qzv \
  --o-sample-frequencies sample-frequencies-ms2.qza \
  --o-feature-frequencies asv-frequencies-ms2.qza
```

## Combining DADA2 Datasets

Each library is processed separately in DADA2 and therefore all ASV tables and representative seqeunce .qza files are combined into a single pair of artifacts. Note that $PFX refers to the parent directory path to the individual libraries with DADA2-processed .qza tables and rep seq files. 

```
# tables
qiime feature-table merge
  --i-tables "$PFX"/lib12/p12.raw_table.qza \
  --i-tables "$PFX"/lib31/p31.raw_table.qza
  --i-tables "$PFX"/lib32/p32.raw_table.qza \
  --i-tables "$PFX"/lib41/p41.raw_table.qza
  --i-tables "$PFX"/lib42/p42.raw_table.qza \
  --i-tables "$PFX"/lib51/p51.raw_table.qza
  --i-tables "$PFX"/lib52/p52.raw_table.qza \
  --i-tables "$PFX"/lib71/p71.raw_table.qza
  --i-tables "$PFX"/lib72/p72.raw_table.qza \
  --o-merged-table tmp.raw_table.qza

# sequences
qiime feature-table merge-seqs
  --i-data "$PFX"/lib12/p12.raw_repSeqs.qza \
  --i-data "$PFX"/lib31/p31.raw_repSeqs.qza
  --i-data "$PFX"/lib32/p32.raw_repSeqs.qza \
  --i-data "$PFX"/lib41/p41.raw_repSeqs.qza
  --i-data "$PFX"/lib42/p42.raw_repSeqs.qza \
  --i-data "$PFX"/lib51/p51.raw_repSeqs.qza
  --i-data "$PFX"/lib52/p52.raw_repSeqs.qza \
  --i-data "$PFX"/lib71/p71.raw_repSeqs.qza
  --i-data "$PFX"/lib72/p72.raw_repSeqs.qza \
  --o-merged-data tmp.raw_repSeqs.qza
```

## Next Steps (O'Rourke)

1) The tmp.raw_table.qza file will serve as input into the contamination overview outlined in the contamination workflow doc from O'Rourke. This includes evaluating seqeunce variants for potential wet-bench cross-contamination and sequencing platform contamination.
2) If there is little evidence of pervasive contamination, then you may proceed to cluster the exact sequence variants (ASVs) into groups of representative variants sharing AT LEAST **98.5%** similarity.
3) Subsequent cluseters are then classified using both alignment and kmer-based approaches, allowing for identification of target taxa.
4) Subsequent OTUs are then analyzed according to diversity metrics. 



### extra stuff

```

qiime demux filter-samples \                                   
--i-demux demux-subsample.qza \     
--m-metadata-file ./demux/subsample/per-sample-fastq-counts.tsv \
--p-where 'CAST([forward sequence count] AS INT) > 100' \
--o-filtered-demux demux.qza

qiime dada2 denoise-paired \                         
--i-demultiplexed-seqs demux-paired.qza \
--p-trim-left-f 10 \                      
--p-trunc-len-f 285 \
--p-trim-left-r 10 \
--p-trunc-len-r 285 \
--p-n-threads 8 \
--o-representative-sequences asv-seqs.qza \
--o-table asv-table.qza \
--o-denoising-stats stats.qza

qiime metadata tabulate \                                      
--m-input-file denoising-stats.qza \
--o-visualization denoising-stats.qzv


```



Filtering sequence data by barcode:
```
qiime tools import \
  --type EMPPairedEndSequences \
  --input-path /path/to/your/fastq_folder/ \
  --output-path emp-paired-end.qza

  qiime demux emp-paired \
  --i-seqs emp-paired-end.qza \
  --m-barcodes-file sample-metadata.tsv \
  --m-barcodes-column BarcodeSequence \
  --o-per-sample-sequences demux-paired.qza \
  --o-error-correction-details demux-paired-details.qza

  qiime demux summarize \
  --i-data demux-paired.qza \
  --o-visualization demux-paired.qzv

  qiime tools export \
  --input-path demux-paired.qza \
  --output-path demultiplexed-sequences
  ```
For running PCA on sequence data that has been split:

```
for file in split_chunk_*; do                                                          
output = "pga_CO1_${file}.txt"
echo "Running CRABS on $file -> $output"
crabs --pairwise-global-alignment --input "$file" --amplicons insilico_CO1_test.txt --output "$output" --forward GGWACWGGWTGAACWGTWTAYCCYCC --reverse TAIACYTCIGGRTGICCRAARAAYCA --size-select 10000 --percent-identity 0.95 --coverage 95 --threads 4
done

declare -A REV=( ["18S"]="GTGTGYACAAAGGBCAGGGAC" ["MiFish"]="CATAGTGGGGTATCTAATCCCAGTTTG")
declare -A AMP=( ["18S"]="pga_18S_test.txt" ["MiFish"]="pga_MiFish_test.txt")
    
THREADS=4
    
for primer in "${"; do
  echo "Running primer set: $primer"
  for file in split_chunk_*.txt; do
    output="pga_${primer}_${file}"
    echo "Running CRABS on $file -> $output"
    crabs --pairwise-global-alignment \
          --input "$file" \
          --amplicons "${AMP[$primer]}" \
          --output "$output" \
          --forward "${FWD[$primer]}" \
          --reverse "${REV[$primer]}" \
          --size-select 10000 \
          --percent-identity 0.95 \
          --coverage 95 \
          --threads $THREADS
  done
done

for file in split_chunk_*; do
output="pga_18S_${file}.txt"    
echo "Running CRABS on $file -> $output"                                     
crabs --pairwise-global-alignment --input "$file" --amplicons insilico_18S_test.txt --output "$output" --forward GYGGTGCATGGCCGTTSKTRGTT --reverse GTGTGYACAAAGGBCAGGGAC --size-select 10000 --percent-identity 0.95 --coverage 95 --threads 4
done                  
  
cat pga_MiFish_split_chunk_*.txt>pga_MiFish_combined.txt

crabs --dereplicate --input pga_MiFish_combined.txt --output MiFish_dereplicated.txt --dereplication-method 'unique_species'

 crabs --export --input MiFish_dereplicated.txt --output MiFish_CTPilot_lib --export-format 'qiime-fasta'                    
# Trash Section


for SAMPLE in $(ls $RAWDIR | cut -f 1 -d '_' | sort -u); do
  cutadapt
  -g mlCO1intF=GGWACWGGWTGAACWGTWTAYCCYCC -g 18SV7&8F=GYGGTGCATGGCCGTTSKTRGTT -g MiFishUF=GTCGGTAAAACTCGTGCCAGC -g   ITS-S2F=ATGCGATACTTGGTGTGAAT -g 16SV3&4F=GTGYCAGCMGCCGCGGTAA
  -G jgHCO2198=TAIACYTCIGGRTGICCRAARAAYCA -G 18SV7&8R=GTGTGYACAAAGGBCAGGGAC -G MiFishUR=CATAGTGGGGTATCTAATCCCAGTTTG   -G ITS-S3R=GACGCTTCTCCAGACTACAAT -G 16SV3&4R=GGACTACNVGGGTWTCTAAT
  -o trimmed-{name}.1.fastq.gz
  -p trimmed-{name}.2.fastq.gz
  input.1.fastq.gz input.2.fastq.gz
  "$RAWDIR"/"$SAMPLE"_R1.fastq.gz "$RAWDIR"/"$SAMPLE"_R2.fastq.gz;
done

  
