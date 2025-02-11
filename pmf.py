#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from matplotlib import pyplot as plt
import pickle
import os, glob
import pymbar
from pymbar import timeseries
from scipy.interpolate import make_interp_spline
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

os.system('mkdir plots')


# In[3]:


num_Rep = len(glob.glob('./delta_tortion.*'))

avg_delta = []
plt.figure(figsize=(20, 3))
for rep in range(num_Rep):
    avg_delta.append(np.mean(np.loadtxt('delta_tortion.'+str(rep+1))[:,1:], axis =1))
    plt.hist(avg_delta[rep], bins=20, density=True, alpha=1, histtype='step', linewidth=1)
    # plt.hist(avg_delta[rep], bins=20, density=True, alpha=0.3)
    plt.xlabel('<delta> [degrees]')
    plt.ylabel('Density')
    
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlim(57,160)
plt.ylim(0,1)

major_locator = MultipleLocator(5)
minor_locator = AutoMinorLocator(5)

plt.gca().xaxis.set_major_locator(major_locator)
plt.gca().xaxis.set_minor_locator(minor_locator)

plt.savefig('plots/umbrella_plot.png', dpi=300, bbox_inches='tight')


# In[4]:


plt.figure(figsize=(10, 4))
for i in range(num_Rep):
    plt.scatter(np.arange(avg_delta[0].shape[0])/10000+1,avg_delta[i], s=2)
    plt.ylabel('<delta> [degrees]', fontsize=12, labelpad=10)
    plt.xlabel('time [ns]', fontsize=12, labelpad=10)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

major_locator = MultipleLocator(1)
minor_locator = AutoMinorLocator(4)

plt.gca().xaxis.set_major_locator(major_locator)
plt.gca().xaxis.set_minor_locator(minor_locator)

plt.xlim(1,10)
plt.savefig('plots/scatter_simulations.png', dpi=300, bbox_inches='tight')


# In[5]:


# get umbrella center given the replica number
def get_umbrella_center(rep_no):
    with open('disang.'+str(rep_no)) as f:
        for line in f:
            if line.startswith('	  r1'):
                x = line.replace(',',' ')
                return float(x.split()[5])
                break


# #### PyMBAR

# In[6]:


kB = 1.381e-23 * 6.022e23 / 1000.0 # Boltzmann constant in kJ/mol/K

temperature = 300.0 # assume a single temperature -- can be overridden with data from center.dat 
# Parameters
K = num_Rep # number of umbrellas
N_max = avg_delta[0].shape[0] # maximum number of snapshots/simulation
beta = 1.0 / (kB * temperature) # inverse temperature of simulations (in 1/(kJ/mol))
delta_min = get_umbrella_center(1) # min for PMF
delta_max = get_umbrella_center(num_Rep) # max for PMF
nbins = 30 # number of bins for 1D PMF

delta0_k = np.array([get_umbrella_center(i+1) for i in range(K)]) # chi0_k[k] is the spring center location (in deg) for umbrella simulation k
delta_kn = np.array([avg_delta]).reshape(K,N_max) # chi_kn[k,n] is the torsion angle (in deg) for snapshot n from umbrella simulation k


# > Used force constant in units kJ/mol/rad**2

# In[7]:


with open('reus_run.py') as f:
    for line in f:
        if line.startswith('f_const'):
            f_const = float(line.split()[2])* 4.184 * (np.pi / 180) ** 2 # convert force constant into kJ/mol/deg**2
            break

#--------------------------------------
# Allocate storage for simulation data
#--------------------------------------

N_k = np.zeros([K], dtype=int) # N_k[k] is the number of snapshots from umbrella simulation k
u_kn = np.zeros([K, N_max]) # u_kn[k,n] is the reduced potential energy without umbrella restraints of snapshot n of umbrella simulation k
g_k = np.zeros([K])


# In[8]:


for k in range(K):
    delta_radians = delta_kn[k,:] / (180.0 / np.pi)
    g_cos = timeseries.statistical_inefficiency(np.cos(delta_radians))
    g_sin = timeseries.statistical_inefficiency(np.sin(delta_radians))
    # print(f"g_cos = {g_cos:.1f} | g_sin = {g_sin:.1f}")
    g_k[k] = max(g_cos, g_sin)
    # print(f"Correlation time for set {k:5d} is {g_k[k]:10.3f}")
    indices = timeseries.subsample_correlated_data(delta_radians, g=g_k[k])
    N_k[k] = len(indices)
    u_kn[k, 0 : N_k[k]] = u_kn[k, indices]
    delta_kn[k, 0 : N_k[k]] = delta_kn[k, indices]

N_max = np.max(N_k)  # shorten the array size


# In[9]:


# compute bin centers

bin_center_i = np.zeros([nbins])
bin_edges = np.linspace(delta_min, delta_max, nbins + 1)
for i in range(nbins):
    bin_center_i[i] = 0.5 * (bin_edges[i] + bin_edges[i + 1])


# In[10]:


# u_kln[k,l,n] is the reduced potential energy of snapshot n from umbrella simulation k evaluated at umbrella l
u_kln = np.zeros([K, K, N_max])

N = np.sum(N_k)
delta_n = pymbar.utils.kn_to_n(delta_kn, N_k=N_k)


# > change retrained energy function accordingly

# In[11]:


# Evaluate reduced energies in all umbrellas
print("Evaluating reduced potential energies...")
for k in range(K):
    for n in range(N_k[k]):
        # Compute minimum-image torsion deviation from umbrella center l
        ddelta = delta_kn[k, n] - delta0_k
        for l in range(K):
            if abs(ddelta[l]) > 180.0:
                ddelta[l] = 360.0 - abs(dcdelta[l])

        # Compute energy of snapshot n from simulation k in umbrella potential l
        # u_kln[k, :, n] = u_kn[k, n] + beta * 24 * f_const * ddelta**2  
        u_kln[k, :, n] = u_kn[k, n] + beta * f_const * ddelta**2  


# In[12]:


# initialize free energy profile with the data collected
fes = pymbar.FES(u_kln, N_k, verbose=True)


# In[14]:


# Compute free energy profile in unbiased potential (in units of kT).
kde_parameters = {}
kde_parameters["bandwidth"] = 0.5 * ((delta_max - delta_min) / nbins)
fes.generate_fes(u_kn, delta_n, fes_type="kde", kde_parameters=kde_parameters)
results = fes.get_fes(bin_center_i, reference_point="from-lowest")
# Write out free energy profile from KDE
center_f_i = results["f_i"]

with open('pmf.dat', 'w') as f:
    f.write("free energy profile (in units of kcal/mol), from KDE\n")
    f.write(f"{'bin':>8s} {'f':>8s}"+'\n')
    for i in range(nbins):
        f.write(f"{bin_center_i[i]:8.1f} {center_f_i[i]* kB * temperature / 4.184:8.3f}"+'\n')


# In[15]:


plt.figure(figsize=(5, 4))

# Assuming you have defined these variables elsewhere
x = bin_center_i
y = center_f_i * kB * temperature / 4.184

x_smooth = np.linspace(bin_center_i.min(), x.max(), 200)
spl = make_interp_spline(x, y)
y_smooth = spl(x_smooth)

smooth = plt.plot(x_smooth, y_smooth, color='red', linewidth=0.6)
plt.scatter(x, y, s=10, marker='o')

plt.xlabel('<delta> [degrees]', fontsize=12, labelpad=10)
plt.ylabel('PMF [kcal/mol]', fontsize=12, labelpad=10)

# Set major ticks every 5 units and minor ticks every 1 unit
major_locator = MultipleLocator(1)
minor_locator = AutoMinorLocator(2)

plt.gca().yaxis.set_major_locator(major_locator)
plt.gca().yaxis.set_minor_locator(minor_locator)

plt.grid(which='both', linestyle='--', linewidth=0.5)
# plt.minorticks_on()

plt.xlim(55,160)
plt.ylim(-1,7)

plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)

plt.savefig('plots/pmf.png', dpi=300, bbox_inches='tight')


# In[ ]:


print('complete')

