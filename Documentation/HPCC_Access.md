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

#SBATCH --job-name=core
#SBATCH --partition=mw128
#SBATCH --output=out_core_%j.stdout
#SBATCH --error=err_core_%j.stderr
#SBATCH -n 48
#SBATCH --mem=120G
#SBATCH --exclusive #This means to completely use this node

node=$(hostname)

scratch=/localscratch/core_$node
local_dir="/zfshomes/sacharya01/cis_tilda/"

echo "These is from core excitations"

rm -rf $scratch
mkdir -p $scratch
cd $local_dir
cp *.chk $scratch/
cp core.py $scratch/

cd $scratch
export OMP_NUM_THREADS=48

starttime=$(date '+%Y-%m-%d %H:%M:%S')

python core.py --chk polygly_24a.chk # Instead of python file you could put your Qiime methods here.

endtime=$(date '+%Y-%m-%d %H:%M:%S')

echo
printf '%.0s' {1..30}
printf '%.0s' {1..30}
echo
printf "%-25s %-25s\n" "Start Time" "End Time"
printf "%-25s %-25s\n" "$starttime" "$endtime"
echo 
printf '%.0s' {1..30}
printf '%.0s' {1..30}
echo

echo "Current node :: $node"
echo "Current jobs jobid :: ${SLURM_JOB_ID}"

cp -r  $scratch "$local_dir/outputs/gly24e_core_${SLURM_JOB_ID}"
```
Basic pipeline structure (NEEDS MODIFYING)

```
#!/bin/bash
#SBATCH --job-name=qiime2_amplicon
#SBATCH --partition=mw128  # Adjust to your cluster partition
#SBATCH --output=qiime2.stdout # Always want to be setting output to the command line, aka stdout
#SBATCH --error=qiime2.stderr
#SBATCH -n 48  # QIIME2 can use multiple cores for some steps
#SBATCH --mem 120G  # Memory depends on dataset size

scratch_dir="/localscratch/qiime"
local_dir=$pwd

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
