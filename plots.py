# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/home/dirge3141/.spyder2/.temp.py
"""

# The main purpose of this code is to get summary statistics for the 
# simulation results as functions of the ratios of key parameters

# import packages for script
import pickle
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

def get_load_data(genotypes,abundances,trait_avgs,num_pts):
    pop_load = np.zeros((num_pts,))
    
    for i in range(num_pts):
        numb_of_genotypes = len(abundances[i][0])
        L1 = np.amax((genotypes[i]-np.array([trait_avgs[i]+np.array([[1,0]]) for j in range(numb_of_genotypes)])).dot(np.array([[s1],[s2]])))
        L2 = np.amax((genotypes[i]-np.array([trait_avgs[i]+np.array([[0,1]]) for j in range(numb_of_genotypes)])).dot(np.array([[s1],[s2]])))
        pop_load[i] = max([L1,L2])
    return pop_load

# parameters of simulation
N=1e11; s1=1e-2; s2=1e-2; U1=1e-4; U2=1e-4; trans_pts_end=15
data_name = '_N-10p11_c1-0d01_c2-0d01_U1-1x10pn4_U2-1x10pn4_exp1'
times = []; genotypes = []; abundances = []

# get simulation data and store genotypes as lists since they vary in dimensions over time
data_file=open('./Documents/kgrel2d/data/pythondata/times'+data_name+'.dat')
times = data_file.read().splitlines()
times = np.array(map(float,times))
num_pts = len(times)
data_file.close()
data_file=open('./Documents/kgrel2d/data/pythondata/genotypes'+data_name+'.dat')
genotypes = data_file.read().splitlines()
data_file.close()
data_file=open('./Documents/kgrel2d/data/pythondata/abundances'+data_name+'.dat')
abundances = data_file.read().splitlines()
data_file.close()
del data_file

# clean up mathematica data format
for i in range(num_pts):
    genotypes[i]='genotypes[i]=np.array(['+genotypes[i].replace('\t',',')+'])'
    genotypes[i]=genotypes[i].replace('{','[')
    genotypes[i]=genotypes[i].replace('}',']')
    exec(genotypes[i])
    abundances[i]='abundances[i]=np.array([['+abundances[i].replace('\t',',')+']])'
    exec(abundances[i])

# compute data for use in plots
trait_avgs = np.zeros((num_pts,2))
rel_fitness = genotypes
frequencies = abundances
covariance = np.zeros(np.shape(times))
variances = np.zeros((num_pts,2))
tot_variance = np.zeros(np.shape(times))

for i in range(num_pts):
    frequencies[i] = (1/np.sum(frequencies[i]))*frequencies[i]
    trait_avgs[i] = abundances[i].dot(genotypes[i])[0]
    rel_fitness[i] = genotypes[i]-np.array([trait_avgs[i] for j in range(len(genotypes[i]))])
    variances[i] = (frequencies[i].dot(((rel_fitness[i])**2)))[0]
    covariance[i] = frequencies[i].dot(rel_fitness[i][:,0]*rel_fitness[i][:,1])
    tot_variance[i] = variances[i][0]+variances[i][1]+2*covariance[i]

del rel_fitness
pop_load = get_load_data(genotypes,abundances,trait_avgs,num_pts)
mid_pt_times = 0.5*(times[2:-1]+times[1:-2])
spd_of_evol = np.transpose(np.array([((times[2:-1]-times[1:-2])**(-1)),((times[2:-1]-times[1:-2])**(-1))]))*(trait_avgs[2:-1]-trait_avgs[1:-2])

# dump data into pickle file
pickle_file_name = './Documents/kgrel2d/data/pythondata/pythondata'+data_name+'.pickle'
pickle_file = open(pickle_file_name,'wb') 
pickle.dump([times,genotypes,abundances,trait_avgs,variances,covariance,tot_variance,pop_load,mid_pt_times,spd_of_evol],pickle_file,pickle.HIGHEST_PROTOCOL)
pickle_file.close()

# load data from pickle file
pickle_file_name = './Documents/kgrel2d/data/pythondata/pythondata'+data_name+'.pickle'
pickle_file = open(pickle_file_name,'rb') 
[times,genotypes,abundances,trait_avgs,variances,covariance,tot_variance,pop_load,mid_pt_times,spd_of_evol] = pickle.load(pickle_file)
pickle_file.close()

# paper figures
# figure 1: variances, covariance and rate of adaptation
fig1,(ax1,ax2)=plt.subplots(2,1,figsize=[16,32],sharey=True)
ax1.plot(times[trans_pts_end:-1],(s1**2)*variances[trans_pts_end:-1,0],c="blue")
ax1.plot(times[trans_pts_end:-1],(s2**2)*variances[trans_pts_end:-1,1],c="red")
ax1.plot(times[trans_pts_end:-1],s1*s2*covariance[trans_pts_end:-1],c="black")
ax1.set_title('Variances-Covariances vs Time',size=18)
ax1.set_xlabel('Time (generations)',size=18)
ax1.ticklabel_format(style='sci',axis='x')
ax1.ticklabel_format(style='sci',axis='y')
#----------------------------------------------------------------------------
ax2.plot(mid_pt_times[trans_pts_end-1:-1],spd_of_evol[trans_pts_end-1:-1,0],c = "blue",linestyle="--")
ax2.plot(mid_pt_times[trans_pts_end-1:-1],spd_of_evol[trans_pts_end-1:-1,1],c = "red",linestyle="--")
#ax2.ticklabel_format(style='sci',axis='x')
#ax2.ticklabel_format(style='sci',axis='y')
ax2.set_title('Rate of Adaptation vs Time',size=18)
plt.show()
fig1.savefig('./Documents/kgrel2d/figures/fig1'+data_name+'.pdf')

# figure 2: substituational load resulting from mutation-selection balance
fig2,ax1=plt.subplots(1,1,figsize=[16,16])
ax1.plot(times[trans_pts_end:-1],pop_load[trans_pts_end:-1])
ax1.set_xlabel("Time (generations)")
ax1.ticklabel_format(style='sci',axis='x',scilimits=(0,0))
ax1.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
ax1.set_ylabel("Mean Fitness Class")
plt.show()
fig2.savefig('./Documents/kgrel2d/figures/fig2'+data_name+'.pdf')