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

First we want to demultiplex our reads. These examples are using paired-end sequence data. Each read will be sorted according to primers. 

For command line history: cat ~/.zsh_history

$RAWDIR points to the libraries of stored sequence data.

```
RAWDIR_WW=/Volumes/ROSALIND/WWorth/raw/all_reads/    
```

$OUTDIR points to an output directory where trimmed reads are stored.
```
OUTDIR_WW=/Volumes/ROSALIND/WWorth/WW_trimmed
```
Demultiplexing all of the sequence data from primers using CutAdapt.

```

RAWDIR_WW=/Volumes/ROSALIND/WWorth/raw/all_reads
OUTDIR_WW=/Volumes/ROSALIND/WWorth/WW_trimmed

FORWARD_PRIMERS=("mlCO1intF=GGWACWGGWTGAACWGTWTAYCCYCC"
                "18SV7&8F=GYGGTGCATGGCCGTTSKTRGTT"
                "MiFishUF=GTCGGTAAAACTCGTGCCAGC"
                "ITS-S2F=ATGCGATACTTGGTGTGAAT"
                "16SV3&4F=GTGYCAGCMGCCGCGGTAA")

REVERSE_PRIMERS=("jgHCO2198=TAIACYTCIGGRTGICCRAARAAYCA"
                "18SV7&8R=GTGTGYACAAAGGBCAGGGAC"
                "MiFishUR=CATAGTGGGGTATCTAATCCCAGTTTG"
                "ITS-S3R=GACGCTTCTCCAGACTACAAT"
                "16SV3&4R=GGACTACNVGGGTWTCTAAT")

MARKER_NAMES=("CO1" "18S" "MiFishU" "ITS2" "16S")

for MARKER in "${MARKER_NAMES[@]}"; do
    mkdir -p "$OUTDIR_WW/$MARKER"
done

for SAMPLE in $(ls $RAWDIR_WW | grep "_R1.fastq.gz" | sed 's/_R1.fastq.gz//' | sort -u); do
  echo "Processing sample: $SAMPLE"

for i in {1..5}; do
  idx=$((i-1))
  MARKER=${MARKER_NAMES[$idx]}
  FWD=${FORWARD_PRIMERS[$idx]}
  REV=${REVERSE_PRIMERS[$idx]}

  echo "$MARKER : $FWD / $REV"

  cutadapt \
    -g "$FWD" \
    -G "$REV" \
    --cores=8 \
    --discard-untrimmed \
    -o "$OUTDIR_WW/$MARKER/${SAMPLE}_${MARKER}.1.fastq.gz" \
    -p "$OUTDIR_WW/$MARKER/${SAMPLE}_${MARKER}.2.fastq.gz" \
    "$RAWDIR_WW/${SAMPLE}_R1.fastq.gz" \
    "$RAWDIR_WW/${SAMPLE}_R2.fastq.gz"
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
## create manifest file
pwd > "$LIB"_pwd.tmptxt
find . -name "*.gz" | sort -u | cut -d '/' -f 2 > "$LIB"_filenames.tmptxt
cut -f 1 -d "_" "$LIB"_filenames.tmptxt > "$LIB"_col1.tmptxt
wc -l "$LIB"_col1.tmptxt | cut -f 1 -d ' ' > "$LIB"_lines.tmptxt
seq $(echo $(cat "$LIB"_lines.tmptxt)) | xargs -Iz echo $(cat "$LIB"_pwd.tmptxt) > "$LIB"_dirpath.tmptxt
paste "$LIB"_dirpath.tmptxt "$LIB"_filenames.tmptxt -d "/" > "$LIB"_col2.tmptxt
for i in $(cat "$LIB"_filenames.tmptxt); do
if [[ $i == *_1.fastq.gz ]];
then
  echo "forward"
else
  echo "reverse"
fi;
done > "$LIB"_col3.tmptxt
paste "$LIB"_col1.tmptxt "$LIB"_col2.tmptxt "$LIB"_col3.tmptxt -d "," > "$LIB"_manifest.tmptxt
echo 'sample-id,absolute-filepath,direction' | cat - "$LIB"_manifest.tmptxt > ../$(echo "$LIB").manifest.file
rm *.tmptxt

```

Import sequence data [O'Rourke](https://github.com/devonorourke/nhguano/blob/master/docs/sequence_processing.md)

```
qiime tools import \\
--type 'SampleData[PairedEndSequencesWithQuality]' \\
--input-path /Volumes/ROSALIND/eDNA_Pilots/WWorth/ww_manifest.tsv \\
--output-path demux.qza \\
--input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \\
--i-data demux-paired \\
--p-n 50000 \\
--o-visualization demux-sumry.qzv
# But how does QIIME demultiplex everything without the primers???
```
DADA2 Denoise Paired-End

```
## generate repseqs
qiime dada2 denoise-paired \
  --p-n-threads 8
  --p-trunc-len-f 181
  --p-trunc-len-r 181 \
  --i-demultiplexed-seqs "$LIB".demux.qza
  --o-denoising-stats "$LIB".denoisingStats.qza \
  --o-table "$LIB".raw_table.qza
  --o-representative-sequences "$LIB".raw_repSeqs.qza \

## generate summary visualization
qiime metadata tabulate \
  --m-input-file "$LIB".denoisingStats.qza
  --o-visualization "$LIB".denoisingStats.qzv  

```









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

  
