#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import numpy as np
import subprocess
import time

amber = 'module purge; module load cuda/12.2.2  gcc/12.2.0  openmpi/4.1.5 amber/22'
common_script = """
#SBATCH --mail-type=ALL
#SBATCH --mail-user=t.desilva@ufl.edu
#SBATCH --distribution=cyclic:cyclic
#SBATCH --mem-per-cpu=4000
#SBATCH --partition=gpu
"""

os.system('rm delta.* delta_tortion_eq.* disang.* md.* eq.* min.* system.parm7 tleap.in leap.log logfile.* reus.* equlibration.* minimization.* dna.rst7 mdinfo* > /dev/null 2>&1')


# In[2]:


pdb = 'dna.pdb' # Name of the pdb file

nmr_Rest = 'DELTA' # NMR restraint keyword for amberMD

num_Rep = 56 # number of H-replicas
rep_per_node = 8 # number of H-replicas to run on a single node
if (num_Rep%rep_per_node) == 0:
    num_nodes = int(num_Rep/rep_per_node)
else:
    print('Number of replicas should be divisable from number of H-replicas to run on a single node')
    exit(1)

start = 55 # starting restraint value
increment = 2.0 # increment in restraint value
f_const = 150.0 # Force constant (kcal/mol/rad^2)


# In[3]:


try:
    f = open("system.parm7", 'r')
    f.close()
except:
    with open("tleap.in", mode='w') as file:
        file.write('source leaprc.DNA.bsc1 \nsource leaprc.water.spce \n') # define forcefiled and water model
        file.write('dna = loadpdb \"' + pdb + '\"\n') # load pdb
        file.write('solvateBox dna SPCBOX 15.0 \naddIons dna K+ 0 \n') # solvate and ion balance
        file.write('saveamberparm dna system.parm7 dna.rst7 \nquit\n') # save files

os.system(amber+'; tleap -f tleap.in >> leap.log') # make parameter file


# In[4]:


with open("leap.log", "r") as file:
    # Read all lines into a list
    lines = file.readlines()

    # Check if the last line starts with the desired string
    last_line = lines[-1].strip() if lines else ""  # Get the last line, or an empty string if file is empty
    if not last_line.startswith("Exiting LEaP: Errors = 0"):
        print("Error in generating parameter files!.")
        exit(1)
    else:
        print('Parmeter files generated successfully...')


# #### Generate NMR restraint files

# In[5]:


def get_seq(pdb_f):
    with open(pdb_f, 'r') as file:
        resNo = 0
        res = []
        for line in file:
            if line.startswith('ATOM'):
                if resNo != int(line.split()[4]):
                    res.append(line.split()[3])
                    resNo = int(line.split()[4])
    return res


# In[6]:


seq = get_seq(pdb)


# In[7]:


#original code
'''
for i in range(num_Rep):
    with open('delta.'+str(i+1), 'w') as delta:
        for index,res in enumerate(seq):
            delta.write('{:>3} {:>3} {:>8} {:8.2f} {:8.2f}\n'.format(index+1, res, nmr_Rest, start, start))

    os.system(amber+'; makeANG_RST -pdb '+ pdb +' -con '+'delta.'+str(i+1)+' > disang.'+str(i+1))
    os.system('sed -i "s/rk2 =   2.0, rk3 =   2.0/rk2 =   '+str(f_const)+', rk3 =   '+str(f_const)+'/g" disang.'+str(i+1))
    os.system('sed -i "s/r1 =.*/r1 = '+str(start-4)+'.0, r2 = '+str(start)+'.0, r3 = '+str(start)+'.0, r4 = '+str(start+4)+'.0,/g" disang.'+str(i+1))
    start += increment
'''

## Modify the original code to capture more in 94-102: after observing the reaction coordinate
for i in range(num_Rep):
    with open('delta.'+str(i+1), 'w') as delta:
        for index,res in enumerate(seq):
            delta.write('{:>3} {:>3} {:>8} {:8.2f} {:8.2f}\n'.format(index+1, res, nmr_Rest, start, start))

    os.system(amber+'; makeANG_RST -pdb '+ pdb +' -con '+'delta.'+str(i+1)+' > disang.'+str(i+1))
    os.system('sed -i "s/rk2 =   2.0, rk3 =   2.0/rk2 =   '+str(f_const)+', rk3 =   '+str(f_const)+'/g" disang.'+str(i+1))
    os.system('sed -i "s/r1 =.*/r1 = '+str(float(start-4))+', r2 = '+str(float(start))+', r3 = '+str(float(start))+', r4 = '+str(float(start+4))+',/g" disang.'+str(i+1))
    if 105>start>94:
        start += 1.25
    else:
        start += increment


# #### Generate group files for equlibration and REMD

# In[8]:


with open('eq.groupfile', 'w') as eq:
    for i in range(num_Rep):
        eq.write('-O -rem 0 -i eq.'+str(i+1)+'.in -o eq.'+str(i+1)+'.out -c min.ncrst -r eq.'+str(i+1)+'.rst -x eq.'+str(i+1)+'.nc -p system.parm7\n')

with open('md.groupfile', 'w') as md:
    for i in range(num_Rep):
        md.write('-O -rem 3 -i md.'+str(i+1)+'.in -o md.'+str(i+1)+'.out -c eq.'+str(i+1)+'.rst -r md.'+str(i+1)+'.rst -x md.'+str(i+1)+'.nc -p system.parm7\n')


# #### Generate amber input files

# In[9]:


with open("min.in", mode='w') as f:
    f.write('minimization \n &cntrl \n  imin=1, maxcyc=3000, ncyc=500,\n')
    f.write('  cut=8.0, ntb=1, ntp=0,\n')
    f.write('  ntc=2, ntf=2,\n /\n')

for i in range(num_Rep):
    with open('eq.'+str(i+1)+'.in', 'w') as f:
        f.write('Equilibration NVT 200 ps \n') 
        f.write(' &cntrl \n')  
        f.write('  imin=0, irest=0, ntx=1,\n')
        f.write('  nstlim=100000, dt=0.002,\n')
        f.write('  ntt=3, gamma_ln=1.0,\n')
        f.write('  temp0=300.0,\n')
        f.write('  ntc=2, ntf=2, nscm=1000,\n')
        f.write('  ntb=1, ntp=0,\n')
        f.write('  cut=8.0,\n')
        f.write('  ntpr=100, ntwx=1000, ntwr=100000,\n')
        f.write('  nmropt=1, ioutfm=1,\n')
        f.write('/\n &wt type=\'DUMPFREQ\', istep1=50,\n')
        f.write('/\n &wt type=\'END\'\n')
        f.write('/\nDISANG=disang.'+str(i+1)+' \n')
        f.write('DUMPAVE=delta_tortion_eq.'+str(i+1)+'\n')

    with open('md.'+str(i+1)+'.in', 'w') as f:
        f.write('REUS NVT 1ps*10000: 10ns \n') 
        f.write(' &cntrl \n')  
        f.write('  imin=0, irest=1, ntx=5,\n')
        f.write('  nstlim=500, dt=0.002, numexchg=10000,\n')
        f.write('  ntt=3, gamma_ln=1.0,\n')
        f.write('  temp0=300.0,\n')
        f.write('  ntc=2, ntf=2, nscm=1000,\n')
        f.write('  ntb=1, ntp=0,\n')
        f.write('  cut=8.0,\n')
        f.write('  ntpr=100, ntwx=1000, ntwr=100000,\n')
        f.write('  nmropt=1, ioutfm=1,\n')
        f.write('/\n &wt type=\'DUMPFREQ\', istep1=50,\n')
        f.write('/\n &wt type=\'END\'\n')
        f.write('/\nDISANG=disang.'+str(i+1)+' \n')
        f.write('DUMPAVE=delta_tortion.'+str(i+1)+'\n')


# #### minimization, equlibration and REUS 

# In[10]:


## make batch files
with open("min.sh", mode='w') as f:
    f.write('#!/bin/bash \n#SBATCH --job-name=minimization \n#SBATCH --output=minimization.out \n#SBATCH --error=minimization.err')
    f.write(common_script)
    f.write('#SBATCH --time=48:00:00 \n') # change time 
    f.write('#SBATCH --nodes=1 \n#SBATCH --ntasks=1 \n#SBATCH --cpus-per-gpu=1 \n#SBATCH --gpus-per-task=1 \n\n')
    f.write('date;hostname;pwd \n\n')
    f.write('cd $SLURM_SUBMIT_DIR\n\n')
    f.write(amber+'\n\n')
    f.write('pmemd.cuda -O -i min.in -p system.parm7 -c dna.rst7 -r min.ncrst -o min.out -ref dna.rst7\n')

with open("eq.sh", mode='w') as f:
    f.write('#!/bin/bash \n#SBATCH --job-name=equlibration \n#SBATCH --output=equlibration.out \n#SBATCH --error=equlibration.err')
    f.write(common_script)
    f.write('#SBATCH --time=72:00:00 \n') # change time
    f.write('#SBATCH --nodes='+str(num_nodes)+' \n#SBATCH --ntasks='+str(num_Rep)+' \n#SBATCH --ntasks-per-node='+str(rep_per_node)+' \n#SBATCH --cpus-per-task=1\n')
    f.write('#SBATCH --gres=gpu:'+str(rep_per_node)+'\n\n')
    f.write('date;hostname;pwd \n\n')
    f.write('cd $SLURM_SUBMIT_DIR\n\n')
    f.write(amber+'\n\n')
    f.write('srun --mpi=pmix_v3 $AMBERHOME/bin/pmemd.cuda_SPFP.MPI -ng '+str(num_Rep)+' -groupfile eq.groupfile\n')

with open("md.sh", mode='w') as f:
    f.write('#!/bin/bash \n#SBATCH --job-name=reus \n#SBATCH --output=reus.out \n#SBATCH --error=reus.err')
    f.write(common_script)
    f.write('#SBATCH --time=72:00:00 \n') # change time
    f.write('#SBATCH --nodes='+str(num_nodes*2)+' \n#SBATCH --ntasks='+str(num_Rep*2)+' \n#SBATCH --ntasks-per-node='+str(rep_per_node)+' \n#SBATCH --cpus-per-task=2\n')
    f.write('#SBATCH --gres=gpu:'+str(rep_per_node)+'\n\n')
    f.write('date;hostname;pwd \n\n')
    f.write('cd $SLURM_SUBMIT_DIR\n\n')
    f.write(amber+'\n\n')
    f.write('srun --mpi=pmix_v3 $AMBERHOME/bin/pmemd.cuda_SPFP.MPI -ng '+str(num_Rep)+' -groupfile md.groupfile\n')


# #### run scripts 

# In[11]:


os.system('JOBID_min=$(sbatch --parsable min.sh); JOBID_eq=$(sbatch --parsable --dependency=afterok:$JOBID_min eq.sh); sbatch --dependency=afterok:$JOBID_eq md.sh')


# In[12]:


# os.system('sbatch md.sh')


# In[ ]:




