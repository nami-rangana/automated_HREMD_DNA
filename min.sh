#!/bin/bash 
#SBATCH --job-name=minimization 
#SBATCH --output=minimization.out 
#SBATCH --error=minimization.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=t.desilva@ufl.edu
#SBATCH --distribution=cyclic:cyclic
#SBATCH --mem-per-cpu=4000
#SBATCH --partition=gpu
#SBATCH --time=48:00:00 
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-gpu=1 
#SBATCH --gpus-per-task=1 

date;hostname;pwd 

cd $SLURM_SUBMIT_DIR

module purge; module load cuda/12.2.2  gcc/12.2.0  openmpi/4.1.5 amber/22

pmemd.cuda -O -i min.in -p system.parm7 -c dna.rst7 -r min.ncrst -o min.out -ref dna.rst7
