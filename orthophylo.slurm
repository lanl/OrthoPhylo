#!/bin/bash
# --------------------------------------------------------------
### PART 1: Requests resources to run your job.
# --------------------------------------------------------------
### Optional. Set the job name
#SBATCH --job-name=Brucella.8.19.22
### Optional. Set the output filename.
### SLURM reads %x as the job name and %j as the job ID
#SBATCH --output=slurm_out/%x-%j.out
### REQUIRED. Specify the PI group for this job
#SBATCH --account=earlm
### REQUIRED. Set the partition for your job.
#SBATCH --partition=tb
### REQUIRED. Set the number of cores that will be used for this job.
#SBATCH --nodes=1
#SBATCH --ntasks=30
### REQUIRED. Set the memory required for this job.
#SBATCH --mem=200gb
### REQUIRED. Specify the time required for this job, hhh:mm:ss
#SBATCH --time=200:00:00

# Use 1st positional arg to specify running test ($1="TESTER")
bash orthophylo.sh $1
