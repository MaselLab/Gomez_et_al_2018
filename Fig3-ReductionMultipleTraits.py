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

#------------------------first panel--------------------------------------------
#plot fig 2 trait no. effect
[N,s,U] = [1e7,2e-2,1e-5]
per_decr_vrate1 = [pltfun.get_vNsU_perChg(N,s,U,i+1) for i in range(10)]        
traitNo = [np.log(i+1) for i in range(10)]
my_xticks = [-0.25]+[log(i+1) for i in range(11)]
my_xlabel = ['', '1', '2', '', '4', '', '', '', '8', '', '', '']
my_yticks = [0+0.20*i for i in range(6)]
my_ylabel = [str(0+0.20*i) for i in range(len(my_yticks))]

fig = plt.figure(figsize=[5,9])

#fig2g.subplots_adjust(bottom=0.25)
ax = plt.subplot(211)

ax.scatter(traitNo,per_decr_vrate1,c="black",label=r's=$0.02$',s=30,marker="o")        
ax.axhline(linewidth=0.5, color = 'k')      
ax.set_ylabel('Adaptation rate / v(U,N,s)',fontsize=14)
ax.set_xlabel('Number of traits (log scale)',fontsize=14)
ax.yaxis.set_label_coords(-0.14,0.5)
plt.xticks(my_xticks,my_xlabel)
plt.yticks(my_yticks,my_ylabel)
ax.tick_params(labelbottom='on',labelleft='on',labelright='off',axis='both',labelsize=12)
ax.grid(b='off', which='both', axis='both')
ax.set_ylim((0,1.1))
ax.set_xlim(-0.25,np.log(11)) 
#ax.legend(loc=1,ncol=1,fontsize=20,frameon=True,scatterpoints = 1)        
#ax2g.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=2,fontsize=20) 

#--------------------Plot reduction in v1 as a function of U and s-------------
# -------------------Second panel-------------------------------------------


# Basic contour plot
ax=plt.subplot(212)

N=1e9
s = np.asarray([10**(-3)*10**(i/50.0) for i in range(101)])
u = np.asarray([10**(-6)*10**(i/50.0) for i in range(101)])

S, U = np.meshgrid(s, u)
VP = pltfun.get_vNsU_perChg(N,s,U,2)

class nf(float):
    def __repr__(self):
        str = '%.3f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.3f' % self.__float__()
        else:
            return '%.3f' % self.__float__()

my_sticks = [np.log10(10**(-3)*10**(i/10.0)) for i in range(21)]
my_slabel = [r'$10^{-3}$','','','','','','','','','',r'$10^{-2}$','','','','','','','','','',r'$10^{-1}$']
my_uticks = [np.log10(10**(-6)*10**(i/10.0)) for i in range(21)]
my_ulabel = [r'$10^{-6}$','','','','','','','','','',r'$10^{-5}$','','','','','','','','','',r'$10^{-4}$']
my_clevels = np.asarray([(0+5*i)/100.0 for i in range(21)])
CS = ax.contour(np.log10(S), np.log10(U),VP,levels=my_clevels,colors='black')
CS.levels = [100*nf(val) for val in CS.levels]

# Label levels with specially formatted floats
if plt.rcParams["text.usetex"]:
    fmt = r'%i \%%'
else:
    fmt = '%i %%'

ax.clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=13)
ax.set_ylabel(r'Mutation rate',fontsize=14)
ax.set_xlabel(r'Selection coefficient',fontsize=14)

plt.xticks(my_sticks,my_slabel,fontsize=14)
plt.yticks(my_uticks,my_ulabel,fontsize=14)
#fig.text(0.95, 0.5, 'Reduction in v1 (multiples of v(U,N,s))', ha='center', va='center', rotation='vertical',fontsize=20)
       
fig.savefig('./figures/vReductionMultipleTraits-Updated.pdf',bbox_inches='tight')
