# -*- coding: utf-8 -*-
"""
Created on Mon Jan 12 13:27:18 2017
@author: Kevin Gomez (Masel Lab)
Script for plots used in "Modeling How Evolution in One Trait is Affected by Adapatation in Another"
"""

#---------------CHANGE WORKING DIRECTORY FIRST---------------------------------
#cd C://Users/dirge/Documents/kgrel2d/

#-------------------------------------------------------------------------------------------
#---------------FIGURE 0: PLOT OF MARGINAL TRAIT DISTRIBUTIONS------------------------------
#-------------------------------------------------------------------------------------------
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

# figure 4: plot of trait or fitness distribution versus normal distribution
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

[mu, var] = pltfun.get_trait_mean_var(genotypes,abundances,traitno)
[trait_bins,trait_tot] = pltfun.get_1D_proj(genotypes,abundances,traitno)
        
# plot data for figure 2a and 2b 
fig0, ax0 = plt.subplots(1,1,figsize=[8,8])
ax4.hist(trait_bins,len(trait_bins),normed=1,facecolor='green',alpha=0.75,weights=trait_tot)
x = np.asarray([np.min(trait_bins)+(np.max(trait_bins)-np.min(trait_bins))*i/100.0 for i in range(101)])
y = mlab.normpdf(x, mu, np.sqrt(var))
ax0.plot(x, y, 'r--', linewidth=1)
ax0.set_title('1D distribution trait '+str(traitno))
fig0.savefig('./figures/fig0.pdf')

# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------------------
#---------------FIGURE 1: SAMPLE OF 2D DISTRIBUTION-----------------------------------------
#-------------------------------------------------------------------------------------------
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
ax1a.set_xlabel('Beneficial Mutaitons Trait 1',fontsize=18,labelpad=20)
ax1a.set_ylabel('Beneficial Mutaitons Trait 2',fontsize=18,labelpad=10)
ax1a.tick_params(axis='both',labelsize=14)        
cbar.ax.text(2.5,0.65,'Log$_{10}$ of Abundances',rotation=90,fontsize=18)
#ax1.scatter([45],[33])
#fit_line = plt.plot([i+40 for i in range(13)],[-i+38 for i in range(13)],ls="--", c=".3")
fig1a.savefig('./figures/fig1a.pdf')

del times, genotypes, abundances, fit_clss_width, class_xlabels
del class_ylabels, fit_distr_2d, cbar, fig1a, ax1a, N, s, U

#-----------------------------------------------------------------------------------------
#------------------FIGURE 2: %CHANGE IN V WRT PARAMETERS----------------------------------
#-----------------------------------------------------------------------------------------
import matplotlib.ticker as mtick
import pickle
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
import matplotlib.mlab as mlab
import plotfunctions as pltfun

#generate figure 2 plots
        
# load time series data of distrStats from plotdata.py output
pickle_file_name = './data/pythondata/sumdata_exp6.pickle'
pickle_file = open(pickle_file_name,'rb') 
[var, cov, vUthry, v2Uthry, varp, covp, vUthryp, v2Uthryp, NsUparam] = pickle.load(pickle_file)
pickle_file.close()

del var, cov, vUthry, v2Uthry, vUthryp, v2Uthryp
NsUparam = np.asarray(NsUparam)
num_exp = len(NsUparam)

#indices for partitioning parameter array
[start1,start2,start3] = [0,num_exp/3,2*num_exp/3]         
[end1,end2,end3] = [num_exp/3,2*num_exp/3,num_exp]

#-------------------------------------------------------------------------------------------
# figure for change in selection coefficient
NsUparam[start1:end1,1]=np.round(NsUparam[start1:end1,1],4)
[min_par,max_par] = [min(NsUparam[start1:end1,1]),max(NsUparam[start1:end1,1])]
[N,s,U] = NsUparam[start1]
param_set = [min_par + i*(max_par-min_par)/100 for i in range(101)]
param_vchg = [pltfun.get_vNsU_perChg(N,param_set[i],U,2) for i in range(101)]

my_xticks = [0.1e-2, 0.3e-2, 0.5e-2, 0.7e-2, 0.9e-2, 1.1e-2]
my_xlabel = ['0.1', '0.3', '0.5', '0.7', '0.9', '1.1']

my_yticks = [20*i for i in range(9)]
my_ylabel = ['']+[str(20*(i+1))+'%' for i in range(9)]

fig2d, ax2d = plt.subplots(1,1,figsize=[8,8])
#fig2d.subplots_adjust(bottom=0.25)
fig2d.subplots_adjust(left=0.15)        
ax2d.plot(np.asarray(param_set),100*np.asarray(param_vchg),c="black",label='$\Delta$v$_{1}$',linewidth=3.0,linestyle = '-')                
ax2d.scatter(NsUparam[start1:end1,1],100*np.asarray(varp[start1:end1])-100,c="black",label='$\Delta \sigma_1^2$',s=40,marker = 'o')
ax2d.scatter(NsUparam[start1:end1,1],-100*np.asarray(covp[start1:end1]),c="black",label='|$\sigma_{12}$|',s=40,marker = 'D')        
ax2d.axhline(linewidth=0.5, color = 'k')     
ax2d.set_ylabel(r'Scaled Variance-Covariance',fontsize=18)
ax2d.set_xlabel(r'Selection Coefficient',fontsize=18)
ax2d.tick_params(axis='both',labelsize=18)
plt.xticks(my_xticks,my_xlabel)
plt.yticks(my_yticks,my_ylabel)
plt.annotate(r'$\times 10^{-2}$',xy=(505,35),xycoords='figure points',fontsize=20)
plt.annotate('(a)',xy=(115,485),xycoords='figure points',fontsize=20)
ax2d.grid(b='off', which='both', axis='both')
#ax2d.ticklabel_format(style='plain',axis='x')
ax2d.set_ylim((0,160))
ax2d.set_xlim((min_par,max_par))
ax2d.legend(loc=4, ncol=3,fontsize=16,frameon=True)
#ax2d.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=2,fontsize=20)
fig2d.savefig('./figures/fig2d.pdf')

#-------------------------------------------------------------------------------------------
# figure for change in beneficial mutation rate
[min_par,max_par] = [min(NsUparam[start2:end2,2]),max(NsUparam[start2:end2,2])]
[N,s,U] = NsUparam[start2]
param_set = [min_par + i*(max_par-min_par)/100.0 for i in range(101)]
param_vchg = [pltfun.get_vNsU_perChg(N,s,param_set[i],2) for i in range(101)]

my_xticks = [0.1e-5, 0.3e-5, 0.5e-5, 0.7e-5, 0.9e-5, 1.1e-5]
my_xlabel = ['0.1', '0.3', '0.5', '0.7', '0.9', '1.1']

my_yticks = [20*i for i in range(9)]
my_ylabel = ['']+[str(20*(i+1))+'%' for i in range(9)]

fig2e, ax2e = plt.subplots(1,1,figsize=[8,8])
#fig2e.subplots_adjust(bottom=0.25)        
fig2e.subplots_adjust(left=0.15)
ax2e.plot(np.asarray(param_set),100*np.asarray(param_vchg),c="black",label='$\Delta$v$_{1}$',linewidth=3.0,linestyle = '-')
ax2e.scatter(NsUparam[start2:end2,2],100*np.asarray(varp[start2:end2])-100,c="black",label='$\Delta\sigma_1^2$',s=40,marker = 'o')
ax2e.scatter(NsUparam[start2:end2,2],-100*np.asarray(covp[start2:end2]),c="black",label='$|\Delta\sigma_{12}$|',s=40,marker = 'D')
ax2e.axhline(linewidth=0.5, color = 'k')
ax2e.set_ylabel(r'Scaled Variance-Covariance',fontsize=18)
ax2e.set_xlabel(r'Mutation Rate',fontsize=18)
ax2e.tick_params(axis='both',labelsize=18)
plt.xticks(my_xticks,my_xlabel)
plt.yticks(my_yticks,my_ylabel)
plt.annotate(r'$\times 10^{-5}$',xy=(505,35),xycoords='figure points',fontsize=20)
plt.annotate('(b)',xy=(115,485),xycoords='figure points',fontsize=20)
ax2e.grid(b='off', which='both', axis='both')
#ax2e.ticklabel_format(style='plain',axis='x')
ax2e.set_ylim((0,160))
ax2e.set_xlim((min_par,max_par))  
ax2e.legend(loc=4,ncol=3,fontsize=16,frameon=True)        
#ax2e.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=2,fontsize=20)
fig2e.savefig('./figures/fig2e.pdf')

#-------------------------------------------------------------------------------------------
# figure for change in population size
NsUparam[start3:end3,0]=np.round(NsUparam[start3:end3,0],8)
[min_par,max_par] = [min(NsUparam[start3:end3,0]),max(NsUparam[start3:end3,0])]
[N,s,U] = NsUparam[start3]
param_set = [min_par + i*(max_par-min_par)/100 for i in range(101)]
param_vchg = [pltfun.get_vNsU_perChg(param_set[i],s,U,2) for i in range(101)]

my_xticks = [1.0e8, 7.0e8, 13.0e8, 19.0e8, 25.0e8]
my_xlabel = ['1.0', '7.0', '13.0', '19.0', '25.0']

my_yticks = [20*i for i in range(9)]
my_ylabel = ['']+[str(20*(i+1))+'%' for i in range(9)]

fig2f, ax2f = plt.subplots(1,1,figsize=[8,8])
#fig2f.subplots_adjust(bottom=0.25)
fig2f.subplots_adjust(left=0.15)
ax2f.plot(np.asarray(param_set),100*np.asarray(param_vchg),c="black",label='$\Delta$v$_{1}$',linewidth=3.0,linestyle = '-')
ax2f.scatter(NsUparam[start3:end3,0],100*np.asarray(varp[start3:end3])-100,c="black",label='$\Delta\sigma_1^2$',s=40,marker = 'o')
ax2f.scatter(NsUparam[start3:end3,0],-100*np.asarray(covp[start3:end3]),c="black",label='|$\Delta\sigma_{12}$|',s=40,marker = 'D')                
ax2f.axhline(linewidth=0.5, color = 'k')           
ax2f.set_ylabel(r'Scaled Variance-Covariance',fontsize=18)
ax2f.set_xlabel(r'Population Size',fontsize=18)
ax2f.tick_params(axis='both',labelsize=18)
plt.xticks(my_xticks,my_xlabel)
plt.yticks(my_yticks,my_ylabel)
plt.annotate(r'$\times 10^{8}$',xy=(505,35),xycoords='figure points',fontsize=20)
plt.annotate('(c)',xy=(115,485),xycoords='figure points',fontsize=20)
ax2f.grid(b='off', which='both', axis='both')
#ax2f.ticklabel_format(style='plain',axis='x')
ax2f.set_ylim((0,160))
ax2f.set_xlim((min_par,max_par)) 
ax2f.legend(loc=4, ncol=3,fontsize=16,frameon=True)        
#ax2f.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=2,fontsize=20)        
fig2f.savefig('./figures/fig2f.pdf')

#-----------------------------------------------------------------------------------------
#------------------FIGURE 2d: %CHANGE IN V WRT TRAIT NO.----------------------------------
#-----------------------------------------------------------------------------------------
import matplotlib.ticker as mtick
import pickle
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
import matplotlib.mlab as mlab
import plotfunctions as pltfun

#plot fig 2 trait no. effect
[N,s,U] = [1e7,2e-2,1e-5]

per_decr_vrate1 = [100*pltfun.get_vNsU_perChg(N,s,U,i+2) for i in range(9)]        
per_decr_vrate2 = [100*pltfun.get_vNsU_perChg(N,s/4,U,i+2) for i in range(9)]        
per_decr_vrate3 = [100*pltfun.get_vNsU_perChg(N,s/10,U,i+2) for i in range(9)]        

my_xticks = [i+1 for i in range(9)]
my_xlabel = [i+1 for i in range(9)]
my_yticks = [20+10*i for i in range(7)]
my_ylabel = [str(20+10*i)+'%' for i in range(len(my_yticks))]

fig2g, ax2g = plt.subplots(1,1,figsize=[8,8])
#fig2g.subplots_adjust(bottom=0.25)
fig2g.subplots_adjust(left=0.15)
ax2g.scatter(my_xticks,per_decr_vrate1,c="black",label=r's=$2 \times 10^{-2}$',s=60,marker="o")        
ax2g.scatter(my_xticks,per_decr_vrate2,c="black",label=r's=$5 \times 10^{-3}$',s=70,marker="*")
ax2g.scatter(my_xticks,per_decr_vrate3,c="black",label=r's=$2 \times 10^{-3}$',s=50,marker="D")                
ax2g.axhline(linewidth=0.5, color = 'k')      
ax2g.set_ylabel('Scaled Decrease in Rate if Adapt.',fontsize=18)
ax2g.set_xlabel('Number of Traits Added',fontsize=18)
#ax2g.tick_params(axis='both',labelsize=18)
plt.xticks(my_xticks,my_xlabel)
plt.yticks(my_yticks,my_ylabel)
plt.annotate('(d)',xy=(115,485),xycoords='figure points',fontsize=20) 
ax2g.grid(b='off', which='both', axis='both')
ax2g.set_ylim((20,80))
ax2g.set_xlim(0,10) 
ax2g.legend(loc=4,ncol=1,fontsize=16,frameon=True)        
#ax2g.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=2,fontsize=20)        
fig2g.savefig('./figures/fig2g.pdf')

#-------------------------------------------------------------------------------------------
#------------------FIGURE 2: %CHANGE IN V WRT PARAMETERS SHARED AXES----------------------
#-------------------------------------------------------------------------------------------

import matplotlib.ticker as mtick
import pickle
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
import matplotlib.mlab as mlab
import plotfunctions as pltfun

# load time series data of distrStats from plotdata.py output
pickle_file_name = './data/pythondata/sumdata_exp6.pickle'
pickle_file = open(pickle_file_name,'rb') 
[var, cov, vUthry, v2Uthry, varp, covp, vUthryp, v2Uthryp, NsUparam] = pickle.load(pickle_file)
pickle_file.close()

del var, cov, vUthry, v2Uthry, vUthryp, v2Uthryp
NsUparam = np.asarray(NsUparam)
num_exp = len(NsUparam)

#indices for partitioning parameter array
[start1,start2,start3] = [0,num_exp/3,2*num_exp/3]         
[end1,end2,end3] = [num_exp/3,2*num_exp/3,num_exp]

[s_min,s_max] = [min(NsUparam[start1:end1,1]),max(NsUparam[start1:end1,1])]
[U_min,U_max] = [min(NsUparam[start2:end2,2]),max(NsUparam[start2:end2,2])]
[N_min,N_max] = [min(NsUparam[start3:end3,0]),max(NsUparam[start3:end3,0])]

#change in v1, variance and covariance as function of s
xs = np.asarray([s_min + i*(s_max-s_min)/100 for i in range(101)])
ps = np.asarray(NsUparam[start1:end1,1])
vs = np.asarray([pltfun.get_vNsU_perChg(NsUparam[start1,0],xs[i],NsUparam[start1,2],2) for i in range(101)])
Vs = np.asarray(varp[start1:end1])
Cs = -np.asarray(covp[start1:end1])

#change in v1, variance and covariance as function of U
xU = np.asarray([U_min + i*(U_max-U_min)/100 for i in range(101)])
pU = np.asarray(NsUparam[start2:end2,2])
vU = np.asarray([pltfun.get_vNsU_perChg(NsUparam[start2,0],NsUparam[start2,1],xU[i],2) for i in range(101)])
VU = np.asarray(varp[start2:end2])
CU = -np.asarray(covp[start2:end2])

#change in v1, variance and covariance as function of N
xN = np.asarray([N_min + i*(N_max-N_min)/100 for i in range(101)])
pN = np.asarray(NsUparam[start3:end3,0])
vN = np.asarray([pltfun.get_vNsU_perChg(xN[i],NsUparam[start3,1],NsUparam[start3,2],2) for i in range(101)])
VN = np.asarray(varp[start3:end3])
CN = -np.asarray(covp[start3:end3])

# figure for change in selection coefficient
fig=plt.figure(figsize=[24,8])
ax=plt.subplot(131)
ax.plot(xs,vs,c="black",label='$\Delta v_{1}$',linewidth=3.0,linestyle = '-')                
ax.scatter(ps,Vs,c="white",label='$\sigma_1^2$',s=40.0,marker = 'D')
ax.scatter(ps,Cs,c="black",label='$|\sigma_{12}|$',s=40.0,marker = 'o')        
ax.set_ylabel(r'Multiples of v(U,N,s)',fontsize=20)
ax.set_xlabel(r'Selection Coefficient',fontsize=20)
ax.tick_params(labelsize=20)
ax.set_ylim((0,2.5))
ax.set_xlim((0.75*s_min,1.25*s_max))
ax.set_xscale('log')
ax.xaxis.set_tick_params(which='both',length=5)
ax.yaxis.set_tick_params(which='both',length=5)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.annotate('(a)',xy=(100,395),xycoords='figure points',fontsize=20)

# figure for change in mutation rate
ax=plt.subplot(132)
ax.plot(xU,vU,c="black",label='$\Delta v_{1}$',linewidth=3.0,linestyle = '-')
ax.scatter(pU,VU,c="white",label='$\sigma_1^2$',s=40.0,marker = 'D')
ax.scatter(pU,CU,c="black",label='$|\sigma_{12}|$',s=40.0,marker = 'o')
ax.set_xlabel(r'Mutation Rate',fontsize=20)
ax.set_ylim((0,2.5))
ax.set_xlim((0.75*U_min,1.25*U_max))
ax.set_xscale('log')
ax.tick_params(labelsize=20)
locs,labels = plt.yticks()
plt.yticks(locs,[])
ax.xaxis.set_tick_params(which='both',length=5)
ax.yaxis.set_tick_params(which='both',length=5)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.annotate('(b)',xy=(495,395),xycoords='figure points',fontsize=20)

# figure for change in population size
ax=plt.subplot(133)
ax.plot(xN,vN,c="black",label='$\Delta v_{1}$',linewidth=3.0,linestyle = '-')
ax.scatter(pN,VN,c="white",label='$\sigma_1^2$',s=40.0,marker = 'D')
ax.scatter(pN,CN,c="black",label='$|\sigma_{12}|$',s=40.0,marker = 'o')                
ax.legend(loc=1, ncol=1,fontsize=16)
ax.set_xlabel(r'Population Size',fontsize=18)
ax.set_ylim((0,2.5))
ax.set_xlim((0.75*N_min,1.25*N_max))
ax.set_xscale('log')
ax.tick_params(labelsize=20)
locs,labels = plt.yticks()
plt.yticks(locs,[])
ax.xaxis.set_tick_params(which='both',length=5)
ax.yaxis.set_tick_params(which='both',length=5)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.annotate('(c)',xy=(890,395),xycoords='figure points',fontsize=20)

fig.subplots_adjust(wspace=0.2)
fig.subplots_adjust(bottom=0.25)
plt.tight_layout

plt.savefig('./figures/fig2.pdf',bbox_inches='tight')

#-------------------------------------------------------------------------------------------
#--------------FIGURE 3: FLUCTUATIONS IN VAR AND COV----------------------------------------
#-------------------------------------------------------------------------------------------
import pickle
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
import matplotlib.mlab as mlab
import plotfunctions as pltfun

# figure 3: plot of rate of adaptation, variances, covariance and their means
[N,s,U] = [1e9,1e-2,1e-5]
[sim_start,sim_end,snapshot] = [5e3,4e4,1.313e4]

# load time series data of distrStats from plotdata.py output
pickle_file_name = './data/pythondata/distrStats_N-10p09_c1-0d01_c2-0d01_U1-1x10pn5_U2-1x10pn5_exp1.pickle'
pickle_file = open(pickle_file_name,'rb') 
[times,mean_fit,fit_var,fit_cov,pop_load,dcov_dt,vU_thry,v2U_thry] = pickle.load(pickle_file)
pickle_file.close()

del mean_fit,pop_load,vU_thry,v2U_thry

# select interval of simulation data that will be used for plot
# reduce loaded data to subset corresponding to selected interval
[start_indx,end_indx] = pltfun.get_sample_window(times,sim_start,sim_end)

times = times[start_indx:end_indx]
fit_var = fit_var[start_indx:end_indx]
fit_cov = fit_cov[start_indx:end_indx]

# compute means of the time series data and store as constant arry
rate_adpt = np.asarray(fit_var[:,0])+np.asarray(fit_var[:,1])+2*np.asarray(fit_cov[:])
var1_avg = np.mean(fit_var[:,0])*np.ones(np.shape(times))
var2_avg = np.mean(fit_var[:,1])*np.ones(np.shape(times))
cov_avg = np.mean(fit_cov[:])*np.ones(np.shape(times))
rate_adpt_avg = var1_avg + var2_avg + 2*cov_avg

# plot data for figure 3a and 3b 
fig3a, ax3a = plt.subplots(1,1,figsize=[8,8])
#fig3a.subplots_adjust(bottom=0.2)
ax3b = plt.twinx(ax3a)
ax3a.plot(times,rate_adpt_avg,c="black",label='var$_{fitness}$',linewidth=2.0,linestyle = '-')                
ax3a.plot(times,var1_avg,c="black",label='var($r_1|r_2$)',linewidth=3.0,linestyle = '--')
ax3a.plot(times,cov_avg,c="black",label='cov($r_1$,$r_2$)',linewidth=3.0,linestyle = ':')
ax3a.axhline(linewidth=0.5, color = 'k')
ax3a.set_ylabel('Fitness Variances & Covariance',fontsize=18)
ax3a.set_xlabel('Time (Generations)',fontsize=18)
ax3a.set_ylim((-3e-4,4e-4))
ax3a.set_xlim((1e4,2e4))
ax3a.tick_params(labelbottom='off',labelleft='off',labelright='off',axis='both',labelsize=18)
ax3a.grid(b='off', which='both', axis='both')
ax3a.ticklabel_format(style='plain',axis='both',scilimits=(0,0))
ax3a.legend(loc=4,ncol=2,fontsize=14)
#ax3a.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),fancybox=True, shadow=True, ncol=3,fontsize=20)

#time dependent trajectories of variance and covariance
ax3b.plot(times,rate_adpt,c="black",linestyle="-",linewidth=3.0)        
ax3b.plot(times,fit_var[:,0],c="black",linestyle="--",linewidth=3.0)
ax3b.plot(times,fit_cov[:],c="black",linestyle=":",linewidth=3.0)
ax3b.axhline(linewidth=0.5, color = 'k')        
ax3b.set_ylim((-3e-4,4e-4))
ax3b.set_xlim((1e4,2e4))        
ax3b.ticklabel_format(style='plain',axis='both',scilimits=(0,0))
ax3b.tick_params(labelbottom='off',labelleft='off',labelright='off',axis='both',labelsize=18)

fig3a.savefig('./figures/fig3.pdf')

del start_indx, end_indx
del times, fit_var, fit_cov, rate_adpt
del rate_adpt_avg, var1_avg, var2_avg,cov_avg
del fig3a, ax3a, ax3b, N, s, U

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

import pickle
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
import matplotlib.mlab as mlab
import plotfunctions as pltfun

# figure 3: plot of rate of adaptation, variances, covariance and their means
[N,s,U] = [1e9,1e-2,1e-5]
[sim_start,sim_end,snapshot] = [5e3,4e4,1.313e4]

# dump data into a pickle files
pickle_file_name = './'+folder_location+'data/pythondata/covdata'+data_name+'.pickle'
pickle_file = open(pickle_file_name,'rb') 
[times,lead_cov,nose_cov] = pickle.load(pickle_file)
pickle_file.close()

# dump data into a pickle files
pickle_file_name = './'+folder_location+'data/pythondata/distrStats'+data_name+'.pickle'
pickle_file = open(pickle_file_name,'rb') 
[times,mean_fit,fit_var,fit_cov,pop_load,dcov_dt,vU_thry,v2U_thry] = pickle.load(pickle_file)
pickle_file.close()

#nose2_cov = [lead_cov[i][-5][2] for i in range(len(lead_cov))]

[N,U,s] = [10**9, 2*10**(-5),10**(-2)]
tau_q = ((np.log(s/U))**2)/(s*(2*np.log(N*s)-np.log(s/U)))
q = (2*np.log(N*s))/(np.log(s/U))
times2 = [times[i]+np.floor(q*tau_q) for i in range(len(times))]

[t_off,t_cov,new_times,new_covs,new_ncovs]= get_cov_cov(times,nose_cov,fit_cov,N,s,U)

# Time displaced Covariance for Poster
fig1,ax1a = plt.subplots(1,1,figsize=[8,8])
ax1a.plot(times,fit_cov[:],c="black",linestyle="-",linewidth=1.0)
#ax1a.plot(sub_times,sub_nose_cov[:],c="red",linestyle="-",linewidth=1.0)
ax1a.plot(times2,nose_cov[:],c="red",linestyle="-",linewidth=1.0)
ax1a.set_ylabel('Covariance',fontsize=18)
ax1a.set_xlabel('Time (Generations)',fontsize=18)
ax1a.axhline(linewidth=0.5, color = 'k')        
#ax1a.set_ylim((-3e-4,4e-4))
ax1a.set_xlim((5e3,1.5e4))        
#ax1a.ticklabel_format(style='sci',axis='both',scilimits=(0,0))
ax1a.tick_params(labelbottom='off',labelleft='off',labelright='off',axis='both',labelsize=14)
ax1a.legend()
plt.show()

# Time Correlation
fig2,ax2a = plt.subplots(1,1,figsize=[8,8])
ax2a.plot(t_off,t_cov,c="red",linestyle="-",linewidth=3.0)
#ax2a.axvline(x=np.floor((q-.5)*tau_q),c="blue",linewidth=3.0)
#ax2a.plot(times2,nose_cov[:],c="red",linestyle="-",linewidth=1.0)
ax2a.set_ylabel('Temporal Correlations',fontsize=18)
ax2a.set_xlabel('Time (Generations)',fontsize=18)
ax2a.axhline(linewidth=0.5, color = 'k')        
ax2a.set_ylim((0,1.1*max(t_cov)))
ax2a.set_xlim((0,1.3e3))        
#ax2a.ticklabel_format(style='sci',axis='both',scilimits=(0,0))
ax2a.tick_params(labelbottom='off',labelleft='off',labelright='off',axis='both',labelsize=14)
ax2a.legend()
plt.show()
