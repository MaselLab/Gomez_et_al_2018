# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 13:27:18 2017

@author: kgomez81
"""

# import packages for script
import pickle
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

# set parameters of simulation and create required variables
N=1e9; s1=1e-2; s2=1e-2; U1=1e-5; U2=1e-5;
#vU_thry = s1*s1*(2*np.log(N*s1)-np.log(s1/(1*U1)))/((np.log(s1/(1*U1)))**2)
#v2U_thry = 0.5*s1*s1*(2*np.log(N*s1)-np.log(s1/(2*U1)))/((np.log(s1/(2*U1)))**2)
#tau_est = 0.5*s1/v2U_thry
data_name = '_N-10p09_c1-0d01_c2-0d01_U1-1x10pn5_U2-1x10pn5_exp2'
[times,genotypes,abundances] = [[],[],[]]

def get_sample_window(times,start_time,end_time):
    [num_pts,start_indx,end_indx] = [len(times),0,0]    
    for i in range(len(times)):
        if times[start_indx] <= start_time:
            start_indx = start_indx + 1
        if times[end_indx] <= end_time:
            end_indx = end_indx + 1
    return [start_indx,end_indx]

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# get simulation data and store genotypes as lists since they vary in dimensions over time
data_file=open('./Documents/kgrel2d/data/pythondata/times'+data_name+'.dat')
times = data_file.read().splitlines()
times = np.array(map(float,times))
data_file.close()
data_file=open('./Documents/kgrel2d/data/pythondata/genotypes'+data_name+'.dat')
genotypes = data_file.read().splitlines()
data_file.close()
data_file=open('./Documents/kgrel2d/data/pythondata/abundances'+data_name+'.dat')
abundances = data_file.read().splitlines()
data_file.close()
del data_file
num_pts = len(times)

# clean up mathematica data format
for i in range(num_pts):
    genotypes[i]='genotypes[i]=np.array(['+genotypes[i].replace('\t',',')+'])'
    genotypes[i]=genotypes[i].replace('{','[')
    genotypes[i]=genotypes[i].replace('}',']')
    exec(genotypes[i])
    abundances[i]='abundances[i]=np.array([['+abundances[i].replace('\t',',')+']])'
    exec(abundances[i])
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------    

# compute data for use in plots
[rel_fitness,frequencies] = [genotypes[:],abundances[:]]
[trait_avgs,variances] = [np.zeros((num_pts,2)) for j in range(2)]
[covariance,pop_load,dcov_dt] = [np.zeros(np.shape(times)) for j in range(3)]

for i in range(num_pts):
    numb_genotypes = len(abundances[i][0])    
    frequencies[i] = (1/np.sum(frequencies[i]))*frequencies[i]
    trait_avgs[i] = frequencies[i].dot(genotypes[i])[0]
    rel_fitness[i] = rel_fitness[i]-np.array([trait_avgs[i] for j in range(numb_genotypes)])
    rel_fitness[i] = rel_fitness[i]*np.array([[s1,s2] for j in range(numb_genotypes)])
    variances[i] = (frequencies[i].dot(((rel_fitness[i])**2)))[0]
    covariance[i] = frequencies[i].dot(rel_fitness[i][:,0]*rel_fitness[i][:,1])
    dcov_dt[i] = frequencies[i].dot(rel_fitness[i][:,0]**2*rel_fitness[i][:,1]+rel_fitness[i][:,1]**2*rel_fitness[i][:,0])
    
    L1 = np.amax((rel_fitness[i]-np.array([[s1,0] for j in range(numb_genotypes)])).dot(np.array([[1],[1]])))
    L2 = np.amax((rel_fitness[i]-np.array([[0,s2] for j in range(numb_genotypes)])).dot(np.array([[1],[1]])))
    pop_load[i] = max([L1,L2])
    
# figure for covariances
plt.plot(times,covariance)
