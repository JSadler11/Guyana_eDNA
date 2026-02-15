Connect to Wesleyan VPN or be on Wes Wifi

```
ssh jsadler@cottontail2.wesleyan.edu

```
Enter password
```
which sinfo
```
Should yield
```
/usr/local/slurm/bin/sinfo
```
Submit job:

```
#!/bin/bash
#SBATCH --job-name=meta #Name of Job
#SBATCH --partition=mw128 #Cluster to run from
#SBATCH --output=meta.stdout #File name for standard out
#SBATCH --error=meta.stderr #File name for standard error
#SBATCH -n 48 #Number of processors on node
#SBATCH --mem 120G
#SBATCH --exclusive #Ensures that all processors are on the same node

scratch=/localscratch/meta

rm -rf  $scratch
mkdir -p $scratch
cd
cp meta_dft.gjf $scratch/

cd $scratch
export OMP_NUM_THREADS=48
gdv meta_dft.gjf
```
Basic pipeline structure (NEEDS MODIFYING)

```
#!/bin/bash
#SBATCH --job-name=qiime2_amplicon
#SBATCH --partition=mw128  # Adjust to your cluster partition
#SBATCH --output=qiime2.stdout
#SBATCH --error=qiime2.stderr
#SBATCH -n 16  # QIIME2 can use multiple cores for some steps
#SBATCH --mem 64G  # Memory depends on dataset size
#SBATCH --time=24:00:00  # Set appropriate time limit

# Load QIIME2 (adjust based on your cluster's module system)
module load qiime2  # or: source activate qiime2

# Set working directory
WORKDIR=$HOME/qiime2_analysis
mkdir -p $WORKDIR
cd $WORKDIR

# 1. Import paired-end sequences
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.tsv \
  --input-format PairedEndFastqManifestPhred33V2 \
  --output-path demux-paired-end.qza

# 2. Quality control visualization
qiime demux summarize \
  --i-data demux-paired-end.qza \
  --o-visualization demux-paired-end.qzv

# 3. Denoise with DADA2 (adjust trim/trunc parameters based on quality)
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux-paired-end.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 250 \
  --p-trunc-len-r 250 \
  --p-n-threads 16 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza

# 4. Generate feature table summary
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file metadata.tsv

# 5. Generate representative sequences summary
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

echo "QIIME2 amplicon processing complete!"
```

**Before running, need to create:**
Manifest.tsv
```
sample-id	forward-absolute-filepath	reverse-absolute-filepath
Sample1	/path/to/Sample1_R1.fastq.gz	/path/to/Sample1_R2.fastq.gz
Sample2	/path/to/Sample2_R1.fastq.gz	/path/to/Sample2_R2.fastq.gz
```
