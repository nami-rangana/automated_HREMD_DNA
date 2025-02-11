#!/bin/bash 
#SBATCH --job-name=equlibration 
#SBATCH --output=equlibration.out 
#SBATCH --error=equlibration.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=t.desilva@ufl.edu
#SBATCH --distribution=cyclic:cyclic
#SBATCH --mem-per-cpu=4000
#SBATCH --partition=gpu
#SBATCH --time=72:00:00 
#SBATCH --nodes=7 
#SBATCH --ntasks=56 
#SBATCH --ntasks-per-node=8 
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:8

date;hostname;pwd 

cd $SLURM_SUBMIT_DIR

module purge; module load cuda/12.2.2  gcc/12.2.0  openmpi/4.1.5 amber/22

srun --mpi=pmix_v3 $AMBERHOME/bin/pmemd.cuda_SPFP.MPI -ng 56 -groupfile eq.groupfile
