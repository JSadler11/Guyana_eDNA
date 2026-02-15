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
