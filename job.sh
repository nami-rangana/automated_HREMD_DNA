#!/bin/bash                                                                 
#SBATCH --job-name=jupyter_cpu   # Job name
#SBATCH --mail-type=ALL,END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=t.desilva@ufl.edu
#SBATCH --ntasks=1                # Number of MPI ranks
#SBATCH --cpus-per-task=1           # Number of cores per MPI rank
##SBATCH --gpus-per-task=1
#SBATCH --partition=bigmem
##SBATCH --nodes=1                  # Number of nodes
#SBATCH --mem-per-cpu=36gb          # Memory per processor
#SBATCH --time=48:00:00              # Time limit hrs:min:sec
#SBATCH --output=job-info.log     # Standard output and error log
pwd; hostname; date

source ~/.bashrc

export PATH=~/.conda/envs/myEnv/bin:$PATH

jupyter notebook --ip $(hostname) --no-browser
