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
#--------------FIGURE: CORR BTWN COVARIANCES FRONT MEAN-----------------------------
# ************************************************************************************
# ************************************************************************************
# ************************************************************************************

data_name = '_N-10p09_c1-0d01_c2-0d01_U1-1x10pn5_U2-1x10pn5_exp1'
folder_location = ''
[N,U,s] = [10**9, 2*10**(-5),10**(-2)]

# compute the right offset and construct nose covariance 
tau_q = ((np.log(s/U))**2)/(s*(2*np.log(N*s)-np.log(s/U)))
q = (2*np.log(N*s))/(np.log(s/U))

# load covariance data of front from pickle file
# what is lead_cov data?? lead cov stores the covariance of each const fit line
pickle_file_name = './'+folder_location+'data/pythondata/covdata'+data_name+'.pickle'
pickle_file = open(pickle_file_name,'rb') 
[times,times2,tau_fix_avg,lead_cov,nose_cov,fit_cov,mean_fix_time,
             t_off,t_cov,new_times,new_covs,new_ncovs] = pickle.load(pickle_file)
pickle_file.close()

xline = np.asarray(np.ones([11,1]))
yline = np.asarray([0+i*max(t_cov)/10 for i in range(11)])


my_xticks = [i*0.25 for i in range(8)]
my_xlabel = ['0', '', '0.5', '', '1', '', '1.5', '']
my_yticks = [i*0.1 for i in range(11)]
my_ylabel = ['0', '', '0.2', '', '0.4', '', '0.6', '', '0.8', '', '1']

# figure 6: cross-covariance between bulk and nose
fig,ax = plt.subplots(figsize=[8,8])
ax.plot((1/tau_fix_avg)*np.asarray(t_off),np.asarray(t_cov),c="black",linestyle="-",linewidth=3.0)
ax.plot(xline,yline,c="black",linestyle="--",linewidth=3.0)
ax.set_ylabel(r'Front-Bulk $\sigma_{1,2}$ Correlation',fontsize=28)
ax.set_xlabel(r'Time offset (Multiples of $\tau_{SW}$)',fontsize=28)
ax.axhline(linewidth=0.5, color = 'k')        
ax.set_ylim((0,1))
ax.set_xlim((0,1.3e3/tau_fix_avg))        
ax.xaxis.set_tick_params(which='both',length=5)
ax.yaxis.set_tick_params(which='both',length=5)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.annotate('(b)',xy=(0.87,0.90),xycoords='axes fraction',fontsize=32)
plt.xticks(my_xticks,my_xlabel)
plt.yticks(my_yticks,my_ylabel)

ax.tick_params(labelbottom='on',labelleft='on',labelright='off',axis='both',labelsize=26)
fig.savefig('./figures/CorrelationNoseBulkCovariance.pdf',bbox_inches='tight')

# figure 6b: Time displaced Covariance for Poster
fig,ax = plt.subplots(1,1,figsize=[8,8])
ax.plot(times,fit_cov[:],c="black",linestyle="-",linewidth=1.0)
ax.plot(times2,nose_cov[:],c="red",linestyle="-",linewidth=1.0)
ax.set_ylabel('Covariance',fontsize=18)
ax.set_xlabel('Time (Generations)',fontsize=18)
ax.axhline(linewidth=0.5, color = 'k')        
ax.set_xlim((5e3,1.5e4))        
ax.tick_params(labelbottom='off',labelleft='off',labelright='off',axis='both',labelsize=14)
ax.xaxis.set_tick_params(which='both',length=5)
ax.yaxis.set_tick_params(which='both',length=5)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.legend()
plt.show()
fig.savefig('./figures/fig5d.pdf')
