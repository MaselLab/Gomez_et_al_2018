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
#------------------FIGURE: ReductionMultipleTraits.-------------------------------
# ************************************************************************************
# ************************************************************************************
# ************************************************************************************

#plot fig 2 trait no. effect
[N,s,U] = [1e7,2e-2,1e-5]
per_decr_vrate1 = [1-1*pltfun.get_vNsU_perChg(N,s,U,i+1) for i in range(10)]        
traitNo = [np.log(i+1) for i in range(10)]
my_xticks = [-0.25]+[log(i+1) for i in range(11)]
my_xlabel = ['', '1', '2', '', '4', '', '', '', '8', '', '', '']
my_yticks = [0+0.20*i for i in range(6)]
my_ylabel = [str(0+0.20*i) for i in range(len(my_yticks))]

fig, ax = plt.subplots(1,1,figsize=[8,8])
#fig2g.subplots_adjust(bottom=0.25)
fig.subplots_adjust(left=0.15)
ax.scatter(traitNo,per_decr_vrate1,c="black",label=r's=$0.02$',s=60,marker="o")        
ax.axhline(linewidth=0.5, color = 'k')      
ax.set_ylabel('Adaptation rate / v(U,N,s)',fontsize=20)
ax.set_xlabel('Number of traits (log scale)',fontsize=20)
plt.xticks(my_xticks,my_xlabel)
plt.yticks(my_yticks,my_ylabel)
ax.tick_params(labelbottom='on',labelleft='on',labelright='off',axis='both',labelsize=18)
ax.grid(b='off', which='both', axis='both')
ax.set_ylim((0,1.1))
ax.set_xlim(-0.25,np.log(11)) 
#ax.legend(loc=1,ncol=1,fontsize=20,frameon=True,scatterpoints = 1)        
#ax2g.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=2,fontsize=20)        
fig.savefig('./figures/vReductionMultipleTraits.pdf')
