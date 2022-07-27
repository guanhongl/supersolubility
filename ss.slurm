#!/bin/bash
#SBATCH --job-name=ssorder
#SBATCH --partition=shared
## 3 day max run time for public partitions, except 4 hour max runtime for the sandbox partition
#SBATCH --time=01-12:00:00 ## time format is DD-HH:MM:SS


#SBATCH --cpus-per-task=33
#SBATCH --mem=8G ## max amount of memory per node you require
##SBATCH --core-spec=0 ## Uncomment to allow jobs to request all cores on a node    
#SBATCH --nodes=1
#SBATCH --ntasks=1

#SBATCH --error=output/ssorder-%A.err ## %A - filled with jobid
#SBATCH --output=output/ssorder-%A.out ## %A - filled with jobid

## Useful for remote notification
##SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE,TIME_LIMIT_80
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE
#SBATCH --mail-user=guanhong@hawaii.edu

## All options and environment variables found on schedMD site: http://slurm.schedmd.com/sbatch.html

##module purge 
##./build.sh
export LD_LIBRARY_PATH="/home/guanhong/flint/lib/"
make clean
make
./ss 11 32

