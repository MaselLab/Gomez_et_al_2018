# -*- coding: utf-8 -*-
"""
Created on Mon Jan 12 13:27:18 2017
@author: Kevin Gomez (Masel Lab)
Script for plots used in "Modeling How Evolution in One Trait is Affected by Adapatation in Another"
"""

#---------------CHANGE WORKING DIRECTORY FIRST---------------------------------
#cd C://Users/dirge/Documents/kgrel2d/

from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import multivariate_normal
from numpy import inf
import matplotlib.ticker as mtick
import pickle
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
import matplotlib.mlab as mlab
import plotfunctions as pltfun


# ************************************************************************************
# ************************************************************************************
#               General format parameters
#
#   Label font size = 20, 24 (tex)
#   tick font size = 18
#   tick padding = 15
#   legend = 20 (tex)
#
# ************************************************************************************
# ************************************************************************************
# ************************************************************************************
#--------------FIGURE: FluctuationsStabilityG----------------------------
# ************************************************************************************
# ************************************************************************************
# ************************************************************************************

# ----------------plots for angle of G and eigenvalues--------------------------------
# load time series data of distrStats from plotdata.py output
#pickle_file_name = './data/2dwave_data_time_series_stab_mm-01.pickle'        #old mathematica data
pickle_file_name = './data/2dwave_data_time_series_stab_ml-01.pickle'         #new matlab data
pickle_file = open(pickle_file_name,'rb') 
[times,fit_var,fit_cov,vU_thry,v2U_thry,lambda1,lambda2,Gang,parameters] = pickle.load(pickle_file)
pickle_file.close()

# compute means of the time series data and store as constant arry
[N,s,U] = parameters 
rate_adpt_t1 = np.asarray(fit_var[:,0])+np.asarray(fit_cov[:])
var1_avg = np.mean(fit_var[:,0])*np.ones(np.shape(times))
var2_avg = np.mean(fit_var[:,1])*np.ones(np.shape(times))
cov_avg = np.mean(fit_cov[:])*np.ones(np.shape(times))
rate_adpt_t1_avg = var1_avg + cov_avg

# ------------------------------------------------------------------------------------
# plots for figure 4a - the variance and covariance, and trait 1 rate of adaptation

lsize1=14
lsize2=14

my_xticks = [5000+i*1000 for i in range(11)]
my_xlabel = ['5', '', '7', '', '9', '', '11', '', '13', '', '15']
my_yticks = [(i-5) for i in range(15)]
my_ylabel = ['', '-4', '', '-2', '', '0', '', '2', '', '4', '', '6','', '8', '']

fig = plt.figure(figsize=[5,8.5])

ax=plt.subplot(311)
ax.plot(times,(1/vU_thry)*rate_adpt_t1,c="black",label=r'$v_1$',linestyle="-",linewidth=2.0)        
ax.plot(times,(1/vU_thry)*fit_var[:,0],c="black",label=r'$\sigma_1^2$',linestyle="--",linewidth=2.0)
ax.plot(times,(1/vU_thry)*fit_cov[:],c="black",label=r'$\sigma_{1,2}$',linestyle=":",linewidth=2.0)
ax.set_ylabel(r'$\sigma_1^2$ and $\sigma_{1,2}$ / $v(U,N,s)$',fontsize=lsize1)
ax.yaxis.set_label_coords(-0.12,0.5)
ax.set_xlim((5e3,1.5e4))
ax.set_ylim((-5.5,9.5))
ax.tick_params(labelbottom='on',labelleft='on',labelright='off',axis='both',labelsize=lsize2)
plt.xticks(my_xticks,[])
plt.yticks(my_yticks,my_ylabel)
ax.grid(b='off', which='both', axis='both')
ax.axhline(linewidth=0.5, color = 'k')    
ax.legend(loc='upper center',ncol=3,fontsize=lsize1,handletextpad=0.3,columnspacing =1.3)
plt.annotate('A',xy=(0.92,0.89),xycoords='axes fraction',fontsize=lsize1,weight='bold') 
plt.tight_layout()
#plt.axvspan(5.35e3,6.4e3, color='gray', alpha=0.2)
#plt.axvspan(1.26e4,1.45e4, color='gray', alpha=0.2)


# ------------------------------------------------------------------------------------

my_yticks = [(i-1)/3.0 for i in range(4)]
my_ylabel = [r'-$15^o$',r'$0^o$',r'$15^o$', r'$30^o$']

ax=plt.subplot(312)

ax.plot(times,2*Gang,c="black",linewidth=2.0,linestyle="-",label='Angle')
ax.axhline(linewidth=0.5, color = 'k')
ax.set_ylabel(r'Angle of $2^{nd}$ eigenvector',fontsize=lsize1,color="black")
ax.yaxis.set_label_coords(-0.12,0.5)
ax.set_xlim((5e3,1.5e4))
ax.set_ylim((-0.50,0.80))
ax.tick_params(labelbottom='on',labelleft='on',labelright='off',axis='both',labelsize=lsize2)
ax.legend(loc='upper center',ncol=2,fontsize=lsize1,handletextpad=0.3,columnspacing =1.3)
plt.xticks(my_xticks,[])
plt.yticks(my_yticks,my_ylabel)
ax.grid(b='off', which='both', axis='both')
ax.axhline(linewidth=0.5, color = 'k')    
plt.annotate('B',xy=(0.92,0.89),xycoords='axes fraction',fontsize=lsize1,weight='bold') 

plt.tight_layout()
#plt.axvspan(5.35e3,6.4e3, color='gray', alpha=0.2)
#plt.axvspan(1.26e4,1.45e4, color='gray', alpha=0.2)

# ------------------------------------------------------------------------------------
# plots for figure 4b - The eigenvalues of the G matrix

my_yticks = [i for i in range(5)]
my_ylabel = ['0', '1', '2', '3', '4']

snapshot = [0.6985e4, 0.7565e4, 0.8198e4, 0.8790e4, 0.9486e4, 1.0162e4]

tl1 = np.asarray([snapshot[0] for i in range(11)])
tl2 = np.asarray([snapshot[1] for i in range(11)])
tl3 = np.asarray([snapshot[2] for i in range(11)])
tl4 = np.asarray([snapshot[3] for i in range(11)])
tl5 = np.asarray([snapshot[4] for i in range(11)])
tl6 = np.asarray([snapshot[5] for i in range(11)])
yl = np.asarray([5*i/10 for i in range(11)])

ax=plt.subplot(313)
ax.plot(times,(1/np.mean(lambda1))*lambda1,c="black",linewidth=2.0,linestyle="-",label='$\lambda_1$ / $\overline{\lambda}_1$')
ax.plot(times,(1/np.mean(lambda2))*lambda2,c="black",linewidth=2.0,linestyle=":",label='$\lambda_2$ / $\overline{\lambda}_2$')
ax.plot(tl1,yl,c="blue",linewidth=2.0,linestyle="--")
ax.plot(tl2,yl,c="green",linewidth=2.0,linestyle="--")
ax.plot(tl3,yl,c="red",linewidth=2.0,linestyle="--")
ax.plot(tl4,yl,c="darkcyan",linewidth=2.0,linestyle="--")
ax.plot(tl5,yl,c="magenta",linewidth=2.0,linestyle="--")
ax.plot(tl6,yl,c="darkorange",linewidth=2.0,linestyle="--")
ax.yaxis.set_label_coords(-0.12,0.5)
ax.set_ylabel(r'Normalized eignvalues of G',fontsize=lsize1)
ax.axhline(linewidth=0.5, color = 'k')
ax.set_xlim((5e3,1.5e4))
ax.set_ylim((0,4.5))
ax.legend(loc='upper center',ncol=2,fontsize=lsize1,handletextpad=0.3,columnspacing =1.3)
ax.grid(b='off', which='both', axis='both')
ax.set_xlabel('Time (Thousands of Generations)',fontsize=lsize1)
ax.tick_params(labelbottom='on',labelleft='on',labelright='off',axis='both',labelsize=lsize2)
plt.xticks(my_xticks,my_xlabel)
plt.yticks(my_yticks,my_ylabel)

#plt.axvspan(5.35e3,6.4e3, color='gray', alpha=0.2)
#plt.axvspan(1.26e4,1.45e4, color='gray', alpha=0.2)
#ax.tick_params(labelsize=18)

plt.annotate('C',xy=(0.92,0.89),xycoords='axes fraction',fontsize=lsize1,weight='bold')
#plt.axvspan(3, 6, color='grey', alpha=0.1)
plt.tight_layout()


#fig.savefig('./figures/FluctuationsStabilityG.pdf',bbox_inches='tight')             #old mathematica data
fig.savefig('./figures/FluctuationsStabilityG-updated.pdf',bbox_inches='tight')     #new matlab data
