# -*- coding: utf-8 -*-
"""
Created on Mon Jan 12 13:27:18 2017
@author: Kevin Gomez (Masel Lab)
Script for plots used in "Modeling How Evolution in One Trait is Affected by Adapatation in Another"
"""

#---------------CHANGE WORKING DIRECTORY FIRST---------------------------------
#cd C://Users/dirge/Documents/kgrel2d/

import matplotlib.ticker as mtick
import pickle
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
import matplotlib.mlab as mlab
import plotfunctions as pltfun

# ************************************************************************************
# ************************************************************************************
# ************************************************************************************
#------------------FIGURE: MeanVarianceCovarianceNsU------------------
# ************************************************************************************
# ************************************************************************************
# ************************************************************************************

# Need to change this
#1. add more points
#2. decrease space between plots
#3. change dots to one dot in legend

# -----------------------------------------------------------------------------------
# load time series data of distrStats from plotdata.py output
#pickle_file_name = './data/pythondata/2dwave_data_time_avg_stats_mm-01.pickle'     #old mathematica data
pickle_file_name = './data/2dwave_data_time_avg_stats_ml-01.pickle'      #new matlab data
pickle_file = open(pickle_file_name,'rb') 
[var, cov, vUthry, v2Uthry, varp, covp, vUthryp, v2Uthryp, NsUparam] = pickle.load(pickle_file)
pickle_file.close()

NsUparam = np.asarray(NsUparam)
num_exp = len(NsUparam)

#indices for partitioning parameter array
#end points had to be manually calibrated with additional data
[start1,start2,start3] = [0,51,102]         
[end1,end2,end3] = [51,102,153]

# used a different order than in simulations
[N_min,N_max] = [min(NsUparam[start1:end1,0]),max(NsUparam[start1:end1,0])]
[s_min,s_max] = [min(NsUparam[start2:end2,1]),max(NsUparam[start2:end2,1])]
[U_min,U_max] = [min(NsUparam[start3:end3,2]),max(NsUparam[start3:end3,2])]

#change in v1, variance and covariance as function of N
xN = np.asarray([N_min + i*(N_max-N_min)/100 for i in range(101)])
pN = np.asarray(NsUparam[start1:end1,0])
vN = np.asarray([pltfun.get_vNsU_perChg(xN[i],NsUparam[start1,1],NsUparam[start1,2],2) for i in range(101)])
VN = np.asarray(varp[start1:end1])
CN = -np.asarray(covp[start1:end1])
VarN = np.asarray([pltfun.get_normlzd_thry_indv_var(xN[i],NsUparam[start1,1],NsUparam[start1,2]) for i in range(101)])
CovN = np.asarray([-pltfun.get_normlzd_thry_cov(xN[i],NsUparam[start1,1],NsUparam[start1,2]) for i in range(101)])

#change in v1, variance and covariance as function of s
xs = np.asarray([s_min + i*(s_max-s_min)/100 for i in range(101)])
ps = np.asarray(NsUparam[start2:end2,1])
vs = np.asarray([pltfun.get_vNsU_perChg(NsUparam[start2,0],xs[i],NsUparam[start2,2],2) for i in range(101)])
Vs = np.asarray(varp[start2:end2])
Cs = -np.asarray(covp[start2:end2])
Vars = np.asarray([pltfun.get_normlzd_thry_indv_var(NsUparam[start2,0],xs[i],NsUparam[start2,2]) for i in range(101)])
Covs = np.asarray([-pltfun.get_normlzd_thry_cov(NsUparam[start2,0],xs[i],NsUparam[start2,2]) for i in range(101)])

#change in v1, variance and covariance as function of U
xU = np.asarray([U_min + i*(U_max-U_min)/100 for i in range(101)])
pU = np.asarray(NsUparam[start3:end3,2])
vU = np.asarray([pltfun.get_vNsU_perChg(NsUparam[start3,0],NsUparam[start3,1],xU[i],2) for i in range(101)])
VU = np.asarray(varp[start3:end3])
CU = -np.asarray(covp[start3:end3])
VarU = np.asarray([pltfun.get_normlzd_thry_indv_var(NsUparam[start3,0],NsUparam[start3,1],xU[i]) for i in range(101)])
CovU = np.asarray([-pltfun.get_normlzd_thry_cov(NsUparam[start3,0],NsUparam[start3,1],xU[i]) for i in range(101)])

# figure for change in selection coefficient
fig=plt.figure(figsize=[24,8])
ax=plt.subplot(131)
ax.plot(xs,vs,c="blue",linewidth=3.0,linestyle = '-')
ax.plot(xs,Vars,c="cyan",linewidth=3.0,linestyle = '-')
ax.plot(xs,Covs,c="magenta",linewidth=3.0,linestyle = '-')
ax.scatter(ps,Vs-Cs,c="blue",label='$v_{1}$ sim',s=40.0,marker = 'o')        
ax.scatter(ps,Vs,c="cyan",label='$\sigma_1^2$ sim',s=40.0,marker = 'o')
ax.scatter(ps,Cs,c="magenta",label='$|\sigma_{1,2}|$ sim',s=40.0,marker = 'o')
ax.plot([],[],c="black",label='theory',linewidth=3.0,linestyle = '-')
ax.legend(framealpha=1,frameon=True,loc=2, ncol=2,fontsize=18,numpoints=1,scatterpoints = 1)
ax.set_ylabel(r'Multiples of v(U,N,s)',fontsize=20)
ax.set_xlabel(r'Selection coefficient',fontsize=20)
ax.tick_params(labelsize=18,pad=15)
ax.set_ylim((0,3.5))
ax.set_xlim((0.75*s_min,1.25*s_max))
ax.set_xscale('log')
ax.xaxis.set_tick_params(which='both',length=5)
ax.yaxis.set_tick_params(which='both',length=5)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.annotate('A',xy=(0.91,0.93),xycoords='axes fraction',fontsize=22,weight='bold')

# figure for change in mutation rate
ax=plt.subplot(132)
ax.plot(xU,vU,c="blue",label='$v_{1}$ theory',linewidth=3.0,linestyle = '-')
ax.plot(xU,VarU,c="cyan",label='$\sigma_{1}^2$ theory',linewidth=3.0,linestyle = '-')
ax.plot(xU,CovU,c="magenta",label='$\sigma_{1,2}$ theory',linewidth=3.0,linestyle = '-')
ax.scatter(pU,VU-CU,c="blue",label='$v_{1}$ sim',s=40.0,marker = 'o')        
ax.scatter(pU,VU,c="cyan",label='$\sigma_1^2$',s=40.0,marker = 'o')
ax.scatter(pU,CU,c="magenta",label='$|\sigma_{1,2}|$',s=40.0,marker = 'o')
ax.set_xlabel(r'Mutation rate',fontsize=20)
ax.set_ylim((0,3.5))
ax.set_xlim((0.75*U_min,1.25*U_max))
ax.set_xscale('log')
ax.tick_params(labelsize=18,pad=15)
locs,labels = plt.yticks()
plt.yticks(locs,[])
ax.xaxis.set_tick_params(which='both',length=5)
ax.yaxis.set_tick_params(which='both',length=5)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.annotate('B',xy=(0.91,0.93),xycoords='axes fraction',fontsize=22,weight='bold')

# figure for change in population size
ax=plt.subplot(133)
ax.plot(xN,vN,c="blue",label='$v_{1}$ theory',linewidth=3.0,linestyle = '-')
ax.plot(xN,VarN,c="cyan",label='$\sigma_{1}^2$ theory',linewidth=3.0,linestyle = '-')
ax.plot(xN,CovN,c="magenta",label='$\sigma_{1,2}$ theory',linewidth=3.0,linestyle = '-')
ax.scatter(pN,VN-CN,c="blue",label='$v_{1}$ sim',s=40.0,marker = 'o')        
ax.scatter(pN,VN,c="cyan",label='$\sigma_1^2$',s=40.0,marker = 'o')
ax.scatter(pN,CN,c="magenta",label='$|\sigma_{1,2}|$',s=40.0,marker = 'o')                
ax.set_xlabel(r'Population size',fontsize=20)
ax.set_ylim((0,3.5))
ax.set_xlim((0.75*N_min,1.25*N_max))
ax.set_xscale('log')
ax.tick_params(labelsize=18,pad=15)
locs,labels = plt.yticks()
plt.yticks(locs,[])
ax.xaxis.set_tick_params(which='both',length=5)
ax.yaxis.set_tick_params(which='both',length=5)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.annotate('C',xy=(0.91,0.93),xycoords='axes fraction',fontsize=22,weight='bold')

fig.subplots_adjust(wspace=0.1)
fig.subplots_adjust(bottom=0.25)
plt.tight_layout
#plt.savefig('./figures/MeanVarianceCovarianceNsU.pdf',bbox_inches='tight')     #figure with old mathematica data
plt.savefig('./figures/MeanVarianceCovarianceNsU-updated.pdf',bbox_inches='tight')

