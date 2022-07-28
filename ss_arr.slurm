#!/bin/bash
#SBATCH --job-name=ssorder_arr
#SBATCH --array=1-9
#SBATCH --partition=shared
## 3 day max run time for public partitions, except 4 hour max runtime for the sandbox partition
#SBATCH --time=02-00:00:00 ## time format is DD-HH:MM:SS


#SBATCH --cpus-per-task=21
#SBATCH --mem=8G ## max amount of memory per node you require
##SBATCH --core-spec=0 ## Uncomment to allow jobs to request all cores on a node    
#SBATCH --nodes=1
#SBATCH --ntasks=1

#SBATCH --error=output/ssorder-%A_%a.err ## %A - filled with jobid
#SBATCH --output=output/ssorder-%A_%a.out ## %A - filled with jobid

## Useful for remote notification
##SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE,TIME_LIMIT_80
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE
#SBATCH --mail-user=guanhong@hawaii.edu

## All options and environment variables found on schedMD site: http://slurm.schedmd.com/sbatch.html

##module purge 
##./build.sh
min=$(( SLURM_ARRAY_TASK_ID * 100000000000 + 1 ))
max=$(( ( SLURM_ARRAY_TASK_ID + 1 ) * 100000000000 ))
export LD_LIBRARY_PATH="/home/guanhong/flint/lib/"
make clean
make
./ss $min $max 20

