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
        L1 = np.amax((genotypes[i]-np.array([trait_avgs[i]-np.array([[1,0]]) for j in range(numb_of_genotypes)])).dot(np.array([[s1],[s2]])))
        L2 = np.amax((genotypes[i]-np.array([trait_avgs[i]-np.array([[0,1]]) for j in range(numb_of_genotypes)])).dot(np.array([[s1],[s2]])))
        pop_load[i] = max([L1,L2])
    return pop_load

def get_sample_window(times,start_time,end_time):
    num_pts = len(times)
    start_indx = 0
    end_indx = 0
    for i in range(len(times)):
        if times[start_indx] <= start_time:
            start_indx = start_indx + 1
        if times[end_indx] <= end_time:
            end_indx = end_indx + 1
    return [start_indx,end_indx]
    
def get_characteristic_time(covariance,times,percent_change):
    times_distr = []
    max_cov_change = percent_change*np.abs(np.mean(covariance))
    indx1 = 0
    indx2 = 0
    max_loop = len(covariance)-1
    
    while (indx1 < max_loop):
        indx2 = indx1        
        while (indx2 <= max_loop):
            if (np.abs(covariance[indx1]-covariance[indx2])<= max_cov_change):
                indx2 = indx2 + 1
            else:
                times_distr = times_distr+[times[indx2]-times[indx1]]
                indx2 = max_loop + 1                
        indx1 = indx1 + 1
        
    return times_distr
    
# parameters of simulation
N=1e9; s1=1e-2; s2=1e-2; U1=1e-5; U2=1e-5;
vU_thry = s1*s1*(2*np.log(N*s1)-np.log(s1/(1*U1)))/((np.log(s1/(1*U1)))**2)
v2U_thry = 0.5*s1*s1*(2*np.log(N*s1)-np.log(s1/(2*U1)))/((np.log(s1/(2*U1)))**2)
tau_est = 0.5*s1/v2U_thry
data_name = '_N-10p09_c1-0d01_c2-0d01_U1-1x10pn5_U2-1x10pn5_exp1'
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
    rel_fitness[i] = (genotypes[i]-np.array([trait_avgs[i] for j in range(len(genotypes[i]))]))*np.array([[s1,s2] for j in range(len(genotypes[i]))])
    variances[i] = (frequencies[i].dot(((rel_fitness[i])**2)))[0]
    covariance[i] = frequencies[i].dot(rel_fitness[i][:,0]*rel_fitness[i][:,1])
    tot_variance[i] = variances[i][0]+variances[i][1]+2*covariance[i]

del rel_fitness

pop_load = get_load_data(genotypes,abundances,trait_avgs,num_pts)
mid_pt_times = 0.5*(times[2:-1]+times[1:-2])
spd_of_evol = np.transpose(np.array([((times[2:-1]-times[1:-2])**(-1)),((times[2:-1]-times[1:-2])**(-1))]))*(trait_avgs[2:-1]-trait_avgs[1:-2])

# select indices for plots
[start_indx,end_indx] = get_sample_window(times,1e4,2e4)
unit_array = np.ones(np.shape(times[start_indx:end_indx])) 
var1_avg = np.mean(variances[start_indx:end_indx,0])
var2_avg = np.mean(variances[start_indx:end_indx,1])
cov_avg = np.mean(covariance)

# compute characteristic times of covariance and rate of adaptation
times_distr_cov = get_characteristic_time(covariance,times,0.1)
times_distr_roa = get_characteristic_time(variances[:,0]+covariance,times,0.1)
time_scale_cov = np.median(times_distr_cov)
time_scale_roa = np.median(times_distr_roa)

# dump data into pickle file
pickle_file_name = './Documents/kgrel2d/data/pythondata/pythondata'+data_name+'.pickle'
pickle_file = open(pickle_file_name,'wb') 
pickle.dump([times,genotypes,abundances,trait_avgs,variances,covariance,tot_variance,pop_load
            ,mid_pt_times,spd_of_evol,var1_avg,var2_avg,cov_avg,unit_array,vU_thry,v2U_thry,
            num_pts,times_distr_cov,times_distr_roa],pickle_file,pickle.HIGHEST_PROTOCOL)
pickle_file.close()

# load data from pickle file
pickle_file_name = './Documents/kgrel2d/data/pythondata/pythondata'+data_name+'.pickle'
pickle_file = open(pickle_file_name,'rb') 
[times,genotypes,abundances,trait_avgs,variances,covariance,tot_variance,pop_load
        ,mid_pt_times,spd_of_evol,var1_avg,var2_avg,cov_avg,unit_array,vU_thry
        ,v2U_thry,num_pts,times_distr_cov,times_distr_roa] = pickle.load(pickle_file)
pickle_file.close()

# paper figures

# figure 1: comparison of variances, covariance and rate of adaptation
fig1,ax1=plt.subplots(nrows=1,ncols=1,figsize=[8,8])
ax12=ax1.twinx()
ax1.plot(times[start_indx:end_indx],var1_avg*unit_array,c="blue",label='var($r_1$)='+str(round(np.mean(variances[start_indx:end_indx,0]),7)),linewidth=1.0)
ax1.plot(times[start_indx:end_indx],cov_avg*unit_array,c="green",label='cov($r_1$,$r_2$)='+str(round(np.mean(covariance[start_indx:end_indx]),7)),linewidth=1.0)
ax1.axhline(linewidth=0.5, color = 'k')
ax1.set_ylabel('Variances-Covariances',fontsize=18)
ax1.set_xlabel('Time (Generations)',fontsize=18)
ax1.set_ylim((-3e-4,4e-4))
ax1.set_xlim((1e4,2e4))
ax1.tick_params(axis='both',labelsize=14)
ax1.ticklabel_format(style='sci',axis='both',scilimits=(0,0))
ax1.legend()
ax12.plot(times[start_indx:end_indx],variances[start_indx:end_indx,0],c="blue",linestyle="--")
ax12.plot(times[start_indx:end_indx],covariance[start_indx:end_indx],c="green",linestyle="--")
ax12.set_ylim((-3e-4,4e-4))
ax12.tick_params(axis='both',labelsize=14)
ax12.ticklabel_format(style='sci',axis='both',scilimits=(0,0))
plt.show()
fig1.savefig('./Documents/kgrel2d/figures/fig1'+data_name+'.pdf')

# figure 2: rates of adaptation and variance+covariance
fig2,ax2=plt.subplots(nrows=1,ncols=1,figsize=[8,8])
ax22=ax2.twinx()
ax2.plot(times[start_indx:end_indx],vU_thry*unit_array,c="black",label='v(U)='+str(round(vU_thry,7)),linewidth=2)
ax2.plot(times[start_indx:end_indx],v2U_thry*unit_array,c="black",label='$2^{-1}$v(2U)='+str(round(v2U_thry,7)),linewidth=2,linestyle="-.")
ax2.plot(times[start_indx:end_indx],(var1_avg+cov_avg)*unit_array,c="orange",label='v($r_1$)='+str(round(var1_avg+cov_avg,7)),linewidth=2)
ax2.plot(times[start_indx:end_indx],var1_avg*unit_array,c="blue",label='var($r_1$)='+str(round(np.mean(variances[start_indx:end_indx,0]),7)),linewidth=2)
ax2.plot(times[start_indx:end_indx],cov_avg*unit_array,c="green",label='cov($r_1$,$r_2$)='+str(round(np.mean(covariance[start_indx:end_indx]),7)),linewidth=2)
ax2.axhline(linewidth=0.5, color = 'k')
ax2.set_ylabel('Rate of Adaptation',fontsize=18)
ax2.set_xlabel('Time (Generations)',fontsize=18)
ax2.set_ylim((-1e-4,3e-4))
ax2.set_xlim((1e4,2e4))
ax2.tick_params(axis='both',labelsize=14)
ax2.ticklabel_format(style='sci',axis='both',scilimits=(0,0))
ax2.legend()
ax22.plot(times[start_indx:end_indx],variances[start_indx:end_indx,0]+covariance[start_indx:end_indx],c="orange",linestyle=":")
ax22.set_ylim((-1e-4,3e-4))
ax22.tick_params(axis='both',labelsize=14)
ax22.ticklabel_format(style='sci',axis='both',scilimits=(0,0))
plt.show()
fig2.savefig('./Documents/kgrel2d/figures/fig2'+data_name+'.pdf')
#----------------------------------------------------------------------------

# figure 3: covariance time scale
sample_time1 = 1.3e4
sample_time2 = 1.55e4
y_range = 0.1*np.abs(cov_avg)
[start_indx2,end_indx2] = get_sample_window(times,sample_time1,sample_time1+0.5*time_scale_cov)
[start_indx3,end_indx3] = get_sample_window(times,sample_time2,sample_time2+0.5*time_scale_cov)

fig3,ax3=plt.subplots(1,1,figsize=[8,8])
ax3.plot(times[start_indx:end_indx],covariance[start_indx:end_indx],c='green',label='cov($r_1$,$r_2$)',linestyle=":",linewidth=1.5)
ax3.plot(times[start_indx:end_indx],cov_avg*unit_array,c="green",label='cov($r_1$,$r_2$)='+str(round(np.mean(covariance[start_indx:end_indx]),7)),linewidth=1.0)
ax3.plot(times[start_indx2:end_indx2],covariance[start_indx2:end_indx2],c='red',label='cov($r_1$,$r_2$)',linewidth=3.0)
ax3.plot(times[start_indx3:end_indx3],covariance[start_indx3:end_indx3],c='blue',label='cov($r_1$,$r_2$)',linewidth=3.0)
ax3.set_xlabel('Time (generations)',fontsize=18)
ax3.set_ylim((-3e-4,0e-4))
ax3.ticklabel_format(style='sci',axis='both',scilimits=(0,0))
ax3.set_ylabel("Covariance",fontsize=18)

ax32 = fig3.add_axes([0.45,0.18,0.33,0.15])
ax32.plot(times[start_indx2:end_indx2],covariance[start_indx2:end_indx2],c='red')
ax32.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
ax32.set_ylim((np.min(covariance[start_indx2:end_indx2]),np.min(covariance[start_indx2:end_indx2])+y_range))
ax32.tick_params(axis='both', which='major', labelsize=8)

ax33 = fig3.add_axes([0.45,0.40,0.33,0.15])
ax33.plot(times[start_indx3:end_indx3],covariance[start_indx3:end_indx3],c='blue')
ax33.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
ax33.set_ylim((np.min(covariance[start_indx3:end_indx3]),np.min(covariance[start_indx3:end_indx3])+y_range))
ax33.tick_params(axis='both', which='major', labelsize=8)

plt.show()
fig3.savefig('./Documents/kgrel2d/figures/fig3'+data_name+'.pdf')