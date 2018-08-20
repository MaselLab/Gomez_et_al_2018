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
#---------------FIGURE 1: 2dDistributionSample.-----------------------------------
# ************************************************************************************
# ************************************************************************************
# ************************************************************************************

# figure 2: plot of sampled two dimensional distribution from simulated data
[N,s,U] = [1e9,1e-2,1e-5]
[sim_start,sim_end,snapshot] = [1e4,4e4,1.313e4]

# load time series data of distrStats from plotdata.py output
pickle_file_name = './data/pythondata/timesGenosAbund_N-10p09_c1-0d01_c2-0d01_U1-1x10pn5_U2-1x10pn5_exp1.pickle'
pickle_file = open(pickle_file_name,'rb') 
[times,genotypes,abundances] = pickle.load(pickle_file)
pickle_file.close()

# change the index to plot a different distribution
snapshot_indx = pltfun.get_sample_window(times,sim_start,sim_start)[0]
genotypes = genotypes[snapshot_indx]
abundances = abundances[snapshot_indx]    

#get width of box to enclose 2d distribution
fit_clss_width = np.max([np.max(genotypes[:,0])-np.min(genotypes[:,0])+5,np.max(genotypes[:,1])-np.min(genotypes[:,1])+5])
                         
# use time index indx = 13630
box_dim = [[fit_clss_width,2],[fit_clss_width,2]]
[distr_grid,class_xlabels,class_ylabels] = pltfun.get_2D_distr(genotypes,abundances,box_dim)
distr_grid = np.log10(N*distr_grid)
distr_grid[distr_grid == -inf] = 0
            
# plot figure 1 with general with constructed discrete gaussian
fig1a, ax1a = plt.subplots(1,1,figsize=[10,8])
fit_distr_2d = ax2.pcolormesh(distr_grid.transpose(),cmap=plt.cm.gray_r)
cbar = plt.colorbar(fit_distr_2d)
ax1a.axis('tight')        
ax1a.set_xticks(np.arange(distr_grid.shape[1])+0.5)
ax1a.set_yticks(np.arange(distr_grid.shape[0])+0.5)        
ax1a.set_xticklabels(class_xlabels[1:],rotation=90)
ax1a.set_yticklabels(class_ylabels[1:])        
ax1a.set_xlabel('Beneficial mutaitons trait 1',fontsize=18,labelpad=20)
ax1a.set_ylabel('Beneficial mutaitons trait 2',fontsize=18,labelpad=10)
ax1a.tick_params(axis='both',labelsize=14)        
cbar.ax.text(2.5,0.65,'Log$_{10}$ of Abundances',rotation=90,fontsize=18)
#ax1.scatter([45],[33])
#fit_line = plt.plot([i+40 for i in range(13)],[-i+38 for i in range(13)],ls="--", c=".3")
fig1a.savefig('./figures/2dDistributionSample.pdf')

del times, genotypes, abundances, fit_clss_width, class_xlabels
del class_ylabels, fit_distr_2d, cbar, fig1a, ax1a, N, s, U