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

# parameters of simulation
N=1e11; s1=1e-2; s2=1e-2; U1=1e-4; U2=1e-4; kill_pts=10
times = []; genotypes = []; abundances = []

# get simulation data and store genotypes as lists since they vary in dimensions over time
data_file=open('./Documents/kgrel2d/data/pythondata/times_N-10p11_c1-0d01_c2-0d01_U1-1x10pn4_U2-1x10pn4_exp1.dat')
times = data_file.read().splitlines()
times = np.array(map(float,times))
data_file.close()

data_file=open('./Documents/kgrel2d/data/pythondata/genotypes_N-10p11_c1-0d01_c2-0d01_U1-1x10pn4_U2-1x10pn4_exp1.dat')
genotypes = data_file.read().splitlines()
data_file.close()

data_file=open('./Documents/kgrel2d/data/pythondata/abundances_N-10p11_c1-0d01_c2-0d01_U1-1x10pn4_U2-1x10pn4_exp1.dat')
abundances = data_file.read().splitlines()
data_file.close()

num_pts = len(times)

for i in range(num_pts):
    genotypes[i]='genotypes[i]=np.array(['+genotypes[i].replace('\t',',')+'])'
    genotypes[i]=genotypes[i].replace('{','[')
    genotypes[i]=genotypes[i].replace('}',']')
    exec(genotypes[i])
    abundances[i]='abundances[i]=np.array([['+abundances[i].replace('\t',',')+']])'
    exec(abundances[i])

# compute the variances and covariances from the data
trait_avgs = np.zeros((len(times),2))
relfitness = genotypes
frequencies = abundances
covariance = np.zeros(np.shape(times))
variances = np.zeros((len(times),2))

for i in range(num_pts):
    frequencies[i]=(1/np.sum(frequencies[i]))*frequencies[i]
    trait_avgs[i] = abundances[i].dot(genotypes[i])[0]
    relfitness[i] = genotypes[i]-np.array([trait_avgs[i] for j in range(len(genotypes[i]))])
    variances[i]=(frequencies[i].dot(((relfitness[i])**2)))[0]
    covariance[i]=frequencies[i].dot(relfitness[i][:,0]*relfitness[i][:,1])

midtimes=0.5*(times[2:-1]+times[1:-2])
spdofevol=((times[2:-1]-times[1:-2])**(-1))*()

for i in range(num_pts-kill_pts)
    
# paper figures
# figure 1: variances, covariance and rate of adaptation
fig1,axes1=plt.subplots(1,1,figsize=[16,16])
plt.plot(times[kill_pts:-1],(s1**2)*variances[10:-1,0],c="blue")
plt.plot(times[kill_pts:-1],(s2**2)*variances[10:-1,1],c="red")
plt.plot(times[kill_pts:-1],s1*s2*covariance[10:-1],c="black")
plt.xlabel("time (generations)")
axes1.ticklabel_format(style='sci', useOffset=False)
axes1.set_ylabel("Variance-Covariance (gen$^{-2}$)",size=14)

axes2=axes1.twinx()

plt.show()
fig1.savefig('./Documents/kgrel2d/figures/fig1')

# CALCULATE  THE GEOMETRIC MEANS OF THE DATA
# get all the base data for model with no competitive mutations
#------------------------------------------------------------------------------

plt.plot(times[10:-1],(s1**2)*variances[10:-1,0]+s1*s2*covariance[10:-1],c="blue")
plt.plot(times[10:-1],(s2**2)*variances[10:-1,1]+s1*s2*covariance[10:-1],c="red")
plt.show()

# BURN A SET OF THE TRANSITIONS DATA
#------------------------------------------------------------------------------


#-------------------- The plots of the abs fit stoch chain --------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

fig2,axes2=plt.subplots(1,1,figsize=[8,8])
subfig2=plt.subplot(111)
for i in range(len(data_slices)-2):
    plt.plot(data_slices[i+2][1],data_slices[i+2][2],label="UrUa="+str(data_slices[i+2][0]))
#plt.plot(simgen1[:,0],simgen1[:,2],c="blue")
#plt.plot(simgen1[:,0],simgen1[:,3],c="red")
#plt.plot(simgen1[:,0],mvgavg1,c="black",linewidth=2)
#plt.xlabel("generations")
#plt.ylabel("abs fitness class")
#plt.subplot(312)
#plt.plot(simgen1[:,0],simgen1[:,4],c="blue")
#plt.plot(simgen1[:,0],simgen1[:,5],c="blue")
#plt.plot(simgen1[:,0],simgen1[:,6],c="red")
#plt.xlabel("generations")
#plt.ylabel("rel fitness class")
#plt.subplot(313)
#plt.plot(simgen1[:,0],simgen1[:,7],c="red")
plt.legend(loc=4)
plt.xlabel("Ratio of Rel. & Abs. Fitness Increments")
plt.ylabel("Mean Fitness Class")
plt.show()
fig2.savefig('./UrUaSrSaCurves')