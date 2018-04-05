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
#---------------FIGURE 1: SAMPLE OF 2D DISTRIBUTION-----------------------------------
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

# ************************************************************************************
# ************************************************************************************
# ************************************************************************************
#------------------FIGURE: %CHANGE IN V WRT PARAMETERS SHARED AXES------------------
# ************************************************************************************
# ************************************************************************************
# ************************************************************************************

# Need to change this
#1. add more points
#2. decrease space between plots
#3. change dots to one dot in legend

# load time series data of distrStats from plotdata.py output
pickle_file_name = './data/pythondata/sumdata_exp7.pickle'
pickle_file = open(pickle_file_name,'rb') 
[var, cov, vUthry, v2Uthry, varp, covp, vUthryp, v2Uthryp, NsUparam] = pickle.load(pickle_file)
pickle_file.close()

del var, cov, vUthry, v2Uthry, vUthryp, v2Uthryp
NsUparam = np.asarray(NsUparam)
num_exp = len(NsUparam)

#indices for partitioning parameter array
#end points had to be manually calibrated with additional data
[start1,start2,start3] = [0,41,81]         
[end1,end2,end3] = [40,80,133]

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
ax.plot(xs,vs,c="black",label='$v_{1}$ theory',linewidth=3.0,linestyle = '-')
ax.scatter(ps,Vs-Cs,c="black",label='$v_{1}$ sim',s=40.0,marker = 'o')        
ax.scatter(ps,Vs,c="white",label='$\sigma_1^2$',s=40.0,marker = 'D')
ax.scatter(ps,Cs,c="white",label='$|\sigma_{1,2}|$',s=40.0,marker = 'o')
ax.legend(loc=2, ncol=2,fontsize=18,numpoints=1,scatterpoints = 1)
ax.set_ylabel(r'Multiples of v(U,N,s)',fontsize=20)
ax.set_xlabel(r'Selection coefficient',fontsize=20)
ax.tick_params(labelsize=18,pad=15)
ax.set_ylim((0,2.5))
ax.set_xlim((0.75*s_min,1.25*s_max))
ax.set_xscale('log')
ax.xaxis.set_tick_params(which='both',length=5)
ax.yaxis.set_tick_params(which='both',length=5)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.annotate('(a)',xy=(0.9,0.93),xycoords='axes fraction',fontsize=20)

# figure for change in mutation rate
ax=plt.subplot(132)
ax.plot(xU,vU,c="black",label='$v_{1}$ theory',linewidth=3.0,linestyle = '-')
ax.scatter(pU,VU-CU,c="black",label='$v_{1}$ sim',s=40.0,marker = 'o')        
ax.scatter(pU,VU,c="white",label='$\sigma_1^2$',s=40.0,marker = 'D')
ax.scatter(pU,CU,c="white",label='$|\sigma_{1,2}|$',s=40.0,marker = 'o')
ax.set_xlabel(r'Mutation rate',fontsize=20)
ax.set_ylim((0,2.5))
ax.set_xlim((0.75*U_min,1.25*U_max))
ax.set_xscale('log')
ax.tick_params(labelsize=18,pad=15)
locs,labels = plt.yticks()
plt.yticks(locs,[])
ax.xaxis.set_tick_params(which='both',length=5)
ax.yaxis.set_tick_params(which='both',length=5)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.annotate('(b)',xy=(0.9,0.93),xycoords='axes fraction',fontsize=20)

# figure for change in population size
ax=plt.subplot(133)
ax.plot(xN,vN,c="black",label='$v_{1}$ theory',linewidth=3.0,linestyle = '-')
ax.scatter(pN,VN-CN,c="black",label='$v_{1}$ sim',s=40.0,marker = 'o')        
ax.scatter(pN,VN,c="white",label='$\sigma_1^2$',s=40.0,marker = 'D')
ax.scatter(pN,CN,c="white",label='$|\sigma_{1,2}|$',s=40.0,marker = 'o')                
ax.set_xlabel(r'Population size',fontsize=20)
ax.set_ylim((0,2.5))
ax.set_xlim((0.75*N_min,1.25*N_max))
ax.set_xscale('log')
ax.tick_params(labelsize=18,pad=15)
locs,labels = plt.yticks()
plt.yticks(locs,[])
ax.xaxis.set_tick_params(which='both',length=5)
ax.yaxis.set_tick_params(which='both',length=5)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.annotate('(c)',xy=(0.9,0.93),xycoords='axes fraction',fontsize=20)

fig.subplots_adjust(wspace=0.1)
fig.subplots_adjust(bottom=0.25)
plt.tight_layout
plt.savefig('./figures/MeanVarianceCovarianceNsU.pdf',bbox_inches='tight')

# ************************************************************************************
# ************************************************************************************
# ************************************************************************************
#------------------FIGURE: %CHANGE IN V WRT TRAIT NO.-------------------------------
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

# ************************************************************************************
# ************************************************************************************
# ************************************************************************************
#--------------FIGURE: FLUCTUATIONS IN VAR-COV AND EVALS----------------------------
# ************************************************************************************
# ************************************************************************************
# ************************************************************************************

# ----------------plots for angle of G and eigenvalues--------------------------------
# load time series data of distrStats from plotdata.py output
pickle_file_name = './data/pythondata/Gstability_exp1.pickle'
pickle_file = open(pickle_file_name,'rb') 
[times,fit_var,fit_cov,pop_load,dcov_dt,vU_thry,v2U_thry,var_diff,
     n1,trG,detG,Gmatr,Xmatr,lambda1,lambda2,Gvec,Gval,Gang,parameters] = pickle.load(pickle_file)
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

fig = plt.figure(figsize=[6,14])
my_yticks = [(i-4) for i in range(11)]
my_ylabel = ['-4', '', '-2', '', '0', '', '2', '', '4', '', '6']

ax=plt.subplot(311)
ax.plot(times,(1/vU_thry)*rate_adpt_t1_avg,c="black",label=r'$v_1$',linewidth=2.0,linestyle = '-')                
ax.plot(times,(1/vU_thry)*var1_avg,c="black",label=r'$\sigma_1^2$',linewidth=2.0,linestyle = '--')
ax.plot(times,(1/vU_thry)*cov_avg,c="black",label=r'$\sigma_{1,2}$',linewidth=2.0,linestyle = ':')
ax.plot(times,(1/vU_thry)*rate_adpt_t1,c="black",linestyle="-",linewidth=2.0)        
ax.plot(times,(1/vU_thry)*fit_var[:,0],c="black",linestyle="--",linewidth=2.0)
ax.plot(times,(1/vU_thry)*fit_cov[:],c="black",linestyle=":",linewidth=2.0)
ax.yaxis.set_label_coords(-0.05,0.5)
ax.axhline(linewidth=0.5, color = 'k')
ax.set_ylabel(r'$\sigma_1^2$ and $\sigma_{1,2}$ / $v(U,N,s)$',fontsize=lsize1)
ax.legend(loc='upper center',ncol=3,fontsize=14)
ax.tick_params(labelbottom='on',labelleft='on',labelright='off',axis='both',labelsize=lsize2)
ax.set_ylim((-5,7))
ax.set_xlim((5e3,1.5e4))
plt.xticks(my_xticks,[])
plt.yticks(my_yticks,my_ylabel)
ax.grid(b='off', which='both', axis='both')
ax.axhline(linewidth=0.5, color = 'k')    
plt.annotate('(a)',xy=(0.90,0.90),xycoords='axes fraction',fontsize=lsize1)         

# ------------------------------------------------------------------------------------
# plots for figure 4b - The eigenvalues of the G matrix

my_yticks = [i for i in range(5)]
my_ylabel = ['0', '1', '2', '3', '4']

tl1 = np.asarray([0.870e4 for i in range(11)])
tl2 = np.asarray([0.920e4 for i in range(11)])
tl3 = np.asarray([0.990e4 for i in range(11)])
tl4 = np.asarray([1.050e4 for i in range(11)])
tl5 = np.asarray([1.120e4 for i in range(11)])
tl6 = np.asarray([1.200e4 for i in range(11)])
yl = np.asarray([4*i/10 for i in range(11)])

ax=plt.subplot(312)
ax.plot(times,(1/np.mean(lambda1))*lambda1,c="black",linewidth=2.0,linestyle=":",label='$\lambda_1$ / $\overline{\lambda}_1$')
ax.plot(times,(1/np.mean(lambda2))*lambda2,c="black",linewidth=2.0,linestyle="-",label='$\lambda_2$ / $\overline{\lambda}_2$')
ax.plot(tl1,yl,c="blue",linewidth=2.0,linestyle="--")
ax.plot(tl2,yl,c="green",linewidth=2.0,linestyle="--")
ax.plot(tl3,yl,c="red",linewidth=2.0,linestyle="--")
ax.plot(tl4,yl,c="darkcyan",linewidth=2.0,linestyle="--")
ax.plot(tl5,yl,c="magenta",linewidth=2.0,linestyle="--")
ax.plot(tl6,yl,c="darkorange",linewidth=2.0,linestyle="--")
ax.yaxis.set_label_coords(-0.05,0.5)
#ax.annotate('',xy=(1.05e4,3.9), xycoords='data',xytext=(1.05e4,4.9),arrowprops=dict(facecolor='black', shrink=0.01, width=1, headwidth=8))
#ax.annotate('',xy=(1.12e4,2.1), xycoords='data',xytext=(1.12e4,3.1),arrowprops=dict(facecolor='black', shrink=0.01, width=1, headwidth=8))
#ax.annotate('',xy=(1.2e4,0.2), xycoords='data',xytext=(1.2e4,-0.8),arrowprops=dict(facecolor='black', shrink=0.01, width=1, headwidth=8))
ax.set_ylabel(r'Normalized eignvalues of G',fontsize=lsize1)
ax.axhline(linewidth=0.5, color = 'k')
ax.set_xlim((5e3,1.5e4))
ax.set_ylim((0,4))
ax.legend(loc='upper center',ncol=2,fontsize=lsize1)
ax.grid(b='off', which='both', axis='both')
ax.tick_params(labelbottom='on',labelleft='on',labelright='off',axis='both',labelsize=lsize2)
plt.xticks(my_xticks,[])
plt.yticks(my_yticks,my_ylabel)
plt.annotate('(b)',xy=(0.90,0.90),xycoords='axes fraction',fontsize=lsize1)

#ax.tick_params(labelsize=18)

# ------------------------------------------------------------------------------------
# plots for figure 4c - The angle of the G matrix wrt selection gradient

#my_yticks = [(i-4) for i in range(9)]
#my_ylabel = [r'$-90^o$','', r'$-45^o$', '', r'$0^o$', '', r'$45^o$', '', r'$95^o$']
my_yticks = [(i-2)/2.0 for i in range(5)]
my_ylabel = [r'-$45^o$','', r'$0^o$','', r'$45^o$']

ax=plt.subplot(313)
lg1 = ax.plot(times,(pi/2)*Gang,c="black",linewidth=2.0,linestyle="-",label='Angle')
ax.axhline(linewidth=0.5, color = 'k')
ax.set_ylabel(r'Angle of $1^{st}$ eigenvector',fontsize=lsize1)
ax.yaxis.set_label_coords(-0.05,0.5)
ax.set_xlim((5e3,1.5e4))
ax.set_ylim((-1,1))
ax.grid(b='off', which='both', axis='both')
ax.tick_params(labelbottom='on',labelleft='on',labelright='off',axis='both',labelsize=lsize2)
ax.set_xlabel('Time (Thousands of Generations)',fontsize=lsize1)
ax.legend(loc=3,ncol=2,fontsize=lsize1)
plt.xticks(my_xticks,my_xlabel)
plt.yticks(my_yticks,my_ylabel)

my_yticks = [(i-4) for i in range(9)]
my_ylabel = ['-4', '', '-2', '', '0', '', '2', '', '4']

ax2 = plt.twinx(ax)
#lg2 = ax2.plot(times,(np.mean(fit_cov))*fit_cov**(-1),c="black",linewidth=2.0,linestyle=":",label='$k$ / $\sigma_{1,2}$')
lg2 = ax2.plot(times,(1/vU_thry)*fit_cov[:],c="black",linestyle=":",linewidth=2.0,label='$\sigma_{1,2}$')
ax2.set_ylabel(r'$\sigma_{1,2}$ / $v(U,N,s)$',fontsize=lsize1)
#ax2.yaxis.set_label_coords(1.05,0.5)
ax2.set_xlim((5e3,1.5e4))
ax2.set_ylim((-4,4))
ax2.tick_params(labelsize=lsize2)
plt.yticks(my_yticks,my_ylabel)

lg = lg1+lg2
labs = [l.get_label() for l in lg]
ax.legend(lg, labs, loc='upper center',ncol=2,fontsize=lsize1)

plt.annotate('(c)',xy=(0.90,0.90),xycoords='axes fraction',fontsize=lsize1)
plt.axvspan(5.1e3,6.3e3, color='gray', alpha=0.2)
plt.axvspan(1.19e4,1.25e4, color='gray', alpha=0.2)
plt.axvspan(1.40e4,1.46e4, color='gray', alpha=0.2)
#plt.axvspan(3, 6, color='grey', alpha=0.1)
plt.tight_layout()

fig.savefig('./figures/FluctuationsStabilityG.pdf',bbox_inches='tight')

# ************************************************************************************
# ************************************************************************************
# ************************************************************************************
# Figure: 2d distribution at peaks and valleys for comparison
# ************************************************************************************
# ************************************************************************************
# ************************************************************************************

# figure 2: plot of sampled two dimensional distribution from simulated data
[N,s,U] = [1e9,1e-2,1e-5]
[sim_start,sim_end] = [1e4,4e4]
snapshot = [0.870e4,0.920e4,0.990e4,1.050e4,1.120e4,1.200e4]

# load time series data of distrStats from plotdata.py output
pickle_file_name = './data/pythondata/timesGenosAbund_N-10p09_c1-0d01_c2-0d01_U1-1x10pn5_U2-1x10pn5_exp1.pickle'
pickle_file = open(pickle_file_name,'rb') 
[times,genotypes_all,abundances_all] = pickle.load(pickle_file)
pickle_file.close()

# ------------------------------------------------------------------------------------                                   
# 2D distribution sample at time 8,700 generations

fig=plt.figure(figsize=[12,8])

ax=plt.subplot(231)
snapshot_indx = pltfun.get_sample_window(times,snapshot[0],snapshot[0])[0]
genotypes = genotypes_all[snapshot_indx]
abundances = abundances_all[snapshot_indx]  
fit_clss_width = np.max([np.max(genotypes[:,0])-np.min(genotypes[:,0])+5,np.max(genotypes[:,1])-np.min(genotypes[:,1])+5])
box_dim = [[fit_clss_width,2],[fit_clss_width,2]]
[distr_grid,class_xlabels,class_ylabels,hhf_points] = pltfun.get_2D_distr(genotypes,abundances,box_dim)
distr_grid = np.log10(N*distr_grid)
distr_grid[distr_grid == -inf] = 0
fit_distr_2d = ax.pcolormesh(distr_grid.transpose(),cmap=plt.cm.gray_r,vmin=0, vmax=9*np.log10(10))

# mark classes at the high fitness front
ax.scatter(hhf_points[:,0],hhf_points[:,1],c="pink",marker="s",linewidth='0',s=100)

# get line data for pre-high fitness front
[xl,yl] = get_hifit_front_line(genotypes,40,box_dim)
ax.plot(xl,yl,c="black",linestyle="-")

#cbar = plt.colorbar(fit_distr_2d)
ax.axis('tight')        
ax.set_xticks(np.arange(distr_grid.shape[1])+0.5)
ax.set_yticks(np.arange(distr_grid.shape[0])+0.5)        
#ax.set_xticklabels(class_xlabels[1:],rotation=90)
#ax.set_yticklabels(class_ylabels[1:])
ax.set_xticklabels([])
ax.set_yticklabels([])
#ax.set_xlabel('Beneficial mutaitons trait 1',fontsize=24,labelpad=20)
#ax.set_ylabel('Beneficial mutaitons trait 2',fontsize=24,labelpad=10)
ax.tick_params(axis='both',labelsize=14)        
#cbar.ax.text(2.5,0.65,'Log$_{10}$ of Abundances',rotation=90,fontsize=18)
plt.annotate('t = 8,700',xy=(0.1,0.05),xycoords='axes fraction',fontsize=16,color="Blue")
plt.annotate('(a)',xy=(0.85,0.90),xycoords='axes fraction',fontsize=16) 

# ---------------------------------------------------------------------------------
# 2D distribution sample at time 9,200 generations

ax=plt.subplot(232)
snapshot_indx = pltfun.get_sample_window(times,snapshot[1],snapshot[1])[0]
genotypes = genotypes_all[snapshot_indx]
abundances = abundances_all[snapshot_indx]
fit_clss_width = np.max([np.max(genotypes[:,0])-np.min(genotypes[:,0])+5,np.max(genotypes[:,1])-np.min(genotypes[:,1])+5])
box_dim = [[fit_clss_width,2],[fit_clss_width,2]]
[distr_grid,class_xlabels,class_ylabels,hhf_points] = pltfun.get_2D_distr(genotypes,abundances,box_dim)
distr_grid = np.log10(N*distr_grid)
distr_grid[distr_grid == -inf] = 0
fit_distr_2d = ax.pcolormesh(distr_grid.transpose(),cmap=plt.cm.gray_r,vmin=0, vmax=9*np.log10(10))

# mark classes at the high fitness front
ax.scatter(hhf_points[:,0],hhf_points[:,1],c="pink",marker="s",linewidth='0',s=100)

# get line data for pre-high fitness front
[xl,yl] = get_hifit_front_line(genotypes,40,box_dim)
ax.plot(xl,yl,c="black",linestyle="-")

#cbar = plt.colorbar(fit_distr_2d)
ax.axis('tight')        
ax.set_xticks(np.arange(distr_grid.shape[1])+0.5)
ax.set_yticks(np.arange(distr_grid.shape[0])+0.5)        
#ax.set_xticklabels(class_xlabels[1:],rotation=90)
#ax.set_yticklabels(class_ylabels[1:])
ax.set_xticklabels([])
ax.set_yticklabels([])    
#ax.set_xlabel('Beneficial mutaitons trait 1',fontsize=24,labelpad=20)
#ax.set_ylabel('Beneficial mutaitons trait 2',fontsize=24,labelpad=10)
ax.tick_params(axis='both',labelsize=14)        
#cbar.ax.text(2.5,0.65,'Log$_{10}$ of Abundances',rotation=90,fontsize=18)
plt.annotate('t = 9,200',xy=(0.1,0.05),xycoords='axes fraction',fontsize=16,color="green")
plt.annotate('(b)',xy=(0.85,0.90),xycoords='axes fraction',fontsize=16)

# ---------------------------------------------------------------------------------
# 2D distribution sample at time 9,900 generations

ax=plt.subplot(233)
snapshot_indx = pltfun.get_sample_window(times,snapshot[2],snapshot[2])[0]
genotypes = genotypes_all[snapshot_indx]
abundances = abundances_all[snapshot_indx]  
fit_clss_width = np.max([np.max(genotypes[:,0])-np.min(genotypes[:,0])+5,np.max(genotypes[:,1])-np.min(genotypes[:,1])+5])
box_dim = [[fit_clss_width,2],[fit_clss_width,2]]
[distr_grid,class_xlabels,class_ylabels,hhf_points] = pltfun.get_2D_distr(genotypes,abundances,box_dim)
distr_grid = np.log10(N*distr_grid)
distr_grid[distr_grid == -inf] = 0
fit_distr_2d = ax.pcolormesh(distr_grid.transpose(),cmap=plt.cm.gray_r,vmin=0, vmax=9*np.log10(10))

# mark classes at the high fitness front
ax.scatter(hhf_points[:,0],hhf_points[:,1],c="pink",marker="s",linewidth='0',s=100)

# get line data for pre-high fitness front
[xl,yl] = get_hifit_front_line(genotypes,40,box_dim)
ax.plot(xl,yl,c="black",linestyle="-")

#cbar = plt.colorbar(fit_distr_2d)
ax.axis('tight')        
ax.set_xticks(np.arange(distr_grid.shape[1])+0.5)
ax.set_yticks(np.arange(distr_grid.shape[0])+0.5)        
#ax.set_xticklabels(class_xlabels[1:],rotation=90)
#ax.set_yticklabels(class_ylabels[1:])
ax.set_xticklabels([])
ax.set_yticklabels([])       
#ax.set_xlabel('Beneficial mutaitons trait 1',fontsize=24,labelpad=20)
#ax.set_ylabel('Beneficial mutaitons trait 2',fontsize=24,labelpad=10)
ax.tick_params(axis='both',labelsize=14)        
#cbar.ax.text(2.5,0.65,'Log$_{10}$ of Abundances',rotation=90,fontsize=18)
plt.annotate('t = 9,900',xy=(0.1,0.05),xycoords='axes fraction',fontsize=16,color="red")
plt.annotate('(c)',xy=(0.85,0.90),xycoords='axes fraction',fontsize=16)

# ---------------------------------------------------------------------------------
# 2D distribution sample at time 10,500 generations

ax=plt.subplot(234)
snapshot_indx = pltfun.get_sample_window(times,snapshot[3],snapshot[3])[0]
genotypes = genotypes_all[snapshot_indx]
abundances = abundances_all[snapshot_indx]  
fit_clss_width = np.max([np.max(genotypes[:,0])-np.min(genotypes[:,0])+5,np.max(genotypes[:,1])-np.min(genotypes[:,1])+5])
box_dim = [[fit_clss_width,2],[fit_clss_width,2]]
[distr_grid,class_xlabels,class_ylabels,hhf_points] = pltfun.get_2D_distr(genotypes,abundances,box_dim)
distr_grid = np.log10(N*distr_grid)
distr_grid[distr_grid == -inf] = 0
fit_distr_2d = ax.pcolormesh(distr_grid.transpose(),cmap=plt.cm.gray_r,vmin=0, vmax=9*np.log10(10))

# mark classes at the high fitness front
ax.scatter(hhf_points[:,0],hhf_points[:,1],c="pink",marker="s",linewidth='0',s=100)

# get line data for pre-high fitness front
[xl,yl] = get_hifit_front_line(genotypes,40,box_dim)
ax.plot(xl,yl,c="black",linestyle="-")

#cbar = plt.colorbar(fit_distr_2d)
ax.axis('tight')        
ax.set_xticks(np.arange(distr_grid.shape[1])+0.5)
ax.set_yticks(np.arange(distr_grid.shape[0])+0.5)        
#ax.set_xticklabels(class_xlabels[1:],rotation=90)
#ax.set_yticklabels(class_ylabels[1:])
ax.set_xticklabels([])
ax.set_yticklabels([])
#ax.set_xlabel('Beneficial mutaitons trait 1',fontsize=24,labelpad=20)
#ax.set_ylabel('Beneficial mutaitons trait 2',fontsize=24,labelpad=10)
ax.tick_params(axis='both',labelsize=14)        
#cbar.ax.text(2.5,0.65,'Log$_{10}$ of Abundances',rotation=90,fontsize=18)
plt.annotate('t = 10,500',xy=(0.1,0.05),xycoords='axes fraction',fontsize=16,color="darkcyan")
plt.annotate('(d)',xy=(0.85,0.90),xycoords='axes fraction',fontsize=16) 

# ---------------------------------------------------------------------------------
# 2D distribution sample at time 11,200 generations
ax=plt.subplot(235)
snapshot_indx = pltfun.get_sample_window(times,snapshot[4],snapshot[4])[0]
genotypes = genotypes_all[snapshot_indx]
abundances = abundances_all[snapshot_indx]
fit_clss_width = np.max([np.max(genotypes[:,0])-np.min(genotypes[:,0])+5,np.max(genotypes[:,1])-np.min(genotypes[:,1])+5])
box_dim = [[fit_clss_width,2],[fit_clss_width,2]]
[distr_grid,class_xlabels,class_ylabels,hhf_points] = pltfun.get_2D_distr(genotypes,abundances,box_dim)
distr_grid = np.log10(N*distr_grid)
distr_grid[distr_grid == -inf] = 0
fit_distr_2d = ax.pcolormesh(distr_grid.transpose(),cmap=plt.cm.gray_r,vmin=0, vmax=9*np.log10(10))

# mark classes at the high fitness front
ax.scatter(hhf_points[:,0],hhf_points[:,1],c="pink",marker="s",linewidth='0',s=100)

# get line data for pre-high fitness front
[xl,yl] = get_hifit_front_line(genotypes,40,box_dim)
ax.plot(xl,yl,c="black",linestyle="-")

#cbar = plt.colorbar(fit_distr_2d)
ax.axis('tight')        
ax.set_xticks(np.arange(distr_grid.shape[1])+0.5)
ax.set_yticks(np.arange(distr_grid.shape[0])+0.5)        
#ax.set_xticklabels(class_xlabels[1:],rotation=90)
#ax.set_yticklabels(class_ylabels[1:])
ax.set_xticklabels([])
ax.set_yticklabels([])    
#ax.set_xlabel('Beneficial mutaitons trait 1',fontsize=24,labelpad=20)
#ax.set_ylabel('Beneficial mutaitons trait 2',fontsize=24,labelpad=10)
ax.tick_params(axis='both',labelsize=14)        
#cbar.ax.text(2.5,0.65,'Log$_{10}$ of Abundances',rotation=90,fontsize=18)
plt.annotate('t = 11,200',xy=(0.1,0.05),xycoords='axes fraction',fontsize=16,color="magenta")
plt.annotate('(e)',xy=(0.85,0.90),xycoords='axes fraction',fontsize=16)

# ---------------------------------------------------------------------------------
# 2D distribution sample at time 12,000 generations
ax=plt.subplot(236)
snapshot_indx = pltfun.get_sample_window(times,snapshot[5],snapshot[5])[0]
genotypes = genotypes_all[snapshot_indx]
abundances = abundances_all[snapshot_indx]  
fit_clss_width = np.max([np.max(genotypes[:,0])-np.min(genotypes[:,0])+5,np.max(genotypes[:,1])-np.min(genotypes[:,1])+5])
box_dim = [[fit_clss_width,2],[fit_clss_width,2]]
[distr_grid,class_xlabels,class_ylabels,hhf_points] = pltfun.get_2D_distr(genotypes,abundances,box_dim)
distr_grid = np.log10(N*distr_grid)
distr_grid[distr_grid == -inf] = 0     
fit_distr_2d = ax.pcolormesh(distr_grid.transpose(),cmap=plt.cm.gray_r,vmin=0, vmax=9*np.log10(10))

# mark classes at the high fitness front
ax.scatter(hhf_points[:,0],hhf_points[:,1],c="pink",marker="s",linewidth='0',s=100)

# get line data for pre-high fitness front
[xl,yl] = get_hifit_front_line(genotypes,40,box_dim)
ax.plot(xl,yl,c="black",linestyle="-")

#cbar = plt.colorbar(fit_distr_2d)
ax.axis('tight')        
ax.set_xticks(np.arange(distr_grid.shape[1])+0.5)
ax.set_yticks(np.arange(distr_grid.shape[0])+0.5)        
#ax.set_xticklabels(class_xlabels[1:],rotation=90)
#ax.set_yticklabels(class_ylabels[1:])
ax.set_xticklabels([])
ax.set_yticklabels([])       
#ax.set_xlabel('Beneficial mutaitons trait 1',fontsize=24,labelpad=20)
#ax.set_ylabel('Beneficial mutaitons trait 2',fontsize=24,labelpad=10)
ax.tick_params(axis='both',labelsize=14)        
#cbar.ax.text(2.5,0.65,'Log$_{10}$ of Abundances',rotation=90,fontsize=18)
plt.annotate('t = 12,000',xy=(0.1,0.05),xycoords='axes fraction',fontsize=16,color="darkorange")
plt.annotate('(f)',xy=(0.85,0.90),xycoords='axes fraction',fontsize=16)

#------------------------------------------------------------------------------------

fig.text(0.51, 0.02, 'Beneficial mutations trait 1', ha='center', va='center',fontsize=14)
fig.text(0.1, 0.5, 'Beneficial mutations trait 2', ha='center', va='center', rotation='vertical',fontsize=14)
fig.text(0.95, 0.5, 'Log10 of abundances', ha='center', va='center', rotation='vertical',fontsize=14)
fig.subplots_adjust(right=0.89)
cbar_ax = fig.add_axes([0.90, 0.10, 0.02, 0.8])
fig.colorbar(fit_distr_2d, cax=cbar_ax)

fig.subplots_adjust(wspace=0.05)
#fig.subplots_adjust(bottom=0.25)
#plt.tight_layout()

fig.savefig('./figures/2dBulkBreakup.pdf',bbox_inches='tight')
