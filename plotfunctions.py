# -*- coding: utf-8 -*-
"""
Created on Tue Aug 08 18:08:03 2017
@author: Kevin Gomez (Masel Lab)
Library of functions used in plots.py and plotdata.py
"""


#--------FUNCTIONS REQUIRE PACKAGES LISTED:------------------------------------
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import multivariate_normal
from numpy import inf
import pickle
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
import matplotlib.mlab as mlab
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
def get_sample_window(times,start_time,end_time):
# returns: indeces of times that correspond to start_time and end_time
 
    [num_pts,start_indx,end_indx] = [len(times),0,0]
    
    for i in range(num_pts):
        if times[start_indx] <= start_time:
            start_indx = start_indx + 1
        if times[end_indx] <= end_time:
            end_indx = end_indx + 1
    
    return [start_indx,end_indx]

# -----------------------------------------------------------------------------
def get_trait_mean_var(genotypes,abundances,traitno):
    
    if ((traitno == 1) | (traitno == 2)):
        mean = (genotypes[:,traitno-1].dot(abundances[0]))/sum(abundances[0])
        var = (((genotypes[:,traitno-1]-mean*np.ones(np.shape(genotypes[:,traitno-1])))**2).dot(abundances[0]))/sum(abundances[0])
    if (traitno == 0):
        mean1 = (genotypes[:,0].dot(abundances[0]))/sum(abundances[0])
        mean2 = (genotypes[:,1].dot(abundances[0]))/sum(abundances[0])
        means_arry = np.asarray([[mean1,mean2] for i in range(len(genotypes[:,0]))])
        mean = mean1 + mean2
        var = (abundances.dot((((genotypes - means_arry)**2).dot(np.ones([2,1]))))[0][0])/sum(abundances[0])
    return [mean, var]
    
# -----------------------------------------------------------------------------
def get_1D_proj(genotypes,abundances,traitno):

    if ((traitno == 1) | (traitno == 2)):
        trait_min = np.min(genotypes[:,traitno-1])
        trait_max = np.max(genotypes[:,traitno-1])
    
        trait_classes = [trait_min+i-3 for i in range(trait_max-trait_min+1+6)]
        trait_totals = [0 for i in range(trait_max-trait_min+1+6)]
    
        for i in range(len(genotypes[:,traitno-1])):
            indx = genotypes[i,traitno-1]-trait_min+3
            trait_totals[indx] = trait_totals[indx]+abundances[0][i] 
    
    if (traitno == 0):
        genotype_fitnesses = genotypes[:,0]+genotypes[:,1]
        trait_min = np.min(genotype_fitnesses)
        trait_max = np.max(genotype_fitnesses)
    
        trait_classes = [trait_min+i-3 for i in range(trait_max-trait_min+1+6)]
        trait_totals = [0 for i in range(trait_max-trait_min+1+6)]
    
        for i in range(len(genotype_fitnesses)):
            indx = genotype_fitnesses[i]-trait_min+3
            trait_totals[indx] = trait_totals[indx]+abundances[0][i] 
    
    return [trait_classes, trait_totals]

# -----------------------------------------------------------------------------
def generate_2D_discrete_gauss(distr_grid, min_fit,mean_fit, cov_matrix, N):
# distr_grid =  square array with 1's where the subpop is nonzero
# mean_fit_class = (i,j) in the array that will be the mean of the gaussian
# returns: two lists genotypes & abundances corresponding to desired distribution
    
    [genotypes,abundances] = [[],[]]    
    
    for i in range(len(distr_grid[:,0])):
        for j in range(len(distr_grid[0,:])):
            if(distr_grid[i,j] != 0):
                genotypes = genotypes+[[i+min_fit[0],j+min_fit[1]]]
                abundances = abundances + [multivariate_normal.pdf([i,j],mean_fit,cov_matrix)]
    
    abundances = [(tot_pop_size/sum(abundances))*abundances[i] for i in range(len(abundances))]
    
    return [genotypes,abundances]

# -----------------------------------------------------------------------------
def get_2D_distr(genotypes,abundances,box_dim):
# box_dim = gives array data to bound distr correponding to genotypes & abund.
#           [[width1,margin1],[width2,margin2]]
# returns: an array whose elements are the abundances of the fit classes

    tot_pop_size = sum(abundances[0])  # be careful with your sums of arrays!!!    
    dim1_data = [np.min(genotypes[:,0]),np.max(genotypes[:,0])]
    dim2_data = [np.min(genotypes[:,1]),np.max(genotypes[:,1])]
    
    if((box_dim[0][0] < dim1_data[1]-dim1_data[0]) | (box_dim[1][0] < dim2_data[1]-dim2_data[0])):
        print "Error with box dimensions"
        end()
    
    my_distr = np.zeros([box_dim[0][0],box_dim[1][0]])
    
    for i in range(len(genotypes)):
        indx1 = genotypes[i][0] - dim1_data[0] + box_dim[0][1]
        indx2 = genotypes[i][1] - dim2_data[0] + box_dim[1][1]
        my_distr[indx1,indx2] = max([abundances[0][i]/tot_pop_size,1/tot_pop_size])
    
    xlabels = [(dim1_data[0] - box_dim[0][1]-1 + i) for i in range(box_dim[0][0]+1)]
    ylabels = [(dim2_data[0] - box_dim[1][1]-1 + i) for i in range(box_dim[1][0]+1)]
    
    return [my_distr,xlabels,ylabels]

# -----------------------------------------------------------------------------
def get_cov_by_fitness_line(genotypes,abundances):
    
    mean_fit1 = get_trait_mean_var(genotypes,abundances,1)[0]
    mean_fit2 = get_trait_mean_var(genotypes,abundances,2)[0]
    num_genotypes = len(abundances[0])
    
    fit1D = [genotypes[i,0]+genotypes[i,1] for i in range(num_genotypes)]
    fit1Dshrt = list(set(fit1D))
    fit1Dcovs = []    
    tempcov = 0
    tempfreq = 0
    popsize = sum(abundances[0])
    
    for i in range(len(fit1Dshrt)):
        tempcov = 0
        tempfreq = 0
        for j in range(num_genotypes):
            if(fit1D[j]==fit1Dshrt[i]):
                tempfreq += abundances[0][j]/popsize
                tempcov += (genotypes[j,0]-mean_fit1)*(genotypes[j,1]-mean_fit2)*(abundances[0][j]/popsize)
        fit1Dcovs = fit1Dcovs +[[fit1Dshrt[i],tempfreq,tempcov/tempfreq]]
    
    # should return [[fit_i,p_i,cov_i,] for i = min_fit,...,max_fit]
    return fit1Dcovs

# -----------------------------------------------------------------------------
def generate_figure(figNum,data_name,folder_location,sim_start,sim_end,pop_param,fname,traitno,num_exp):
# figure 1: representation of two-dimensional distribution  (no data required)
    [N,s1,s2,U1,U2] = pop_param
    plt.ioff()
    
    # figure 1: general two dimensional distribution formed as gaussian
    if(figNum == 1):
        
        # set up distribution data for the figure
        [min_fit_clss,mean_fit_clss, cov_mtrx, box_dim] = [[0,0],[7,7],[[1,0],[0,1]],[[20,7],[20,7]]]
        [cut_off, distr_grid] = [1e-4,np.zeros([box_dim[0][0],box_dim[1][0]])]
        [arry_dim1,arry_dim2]=[np.shape(distr_grid[:,0])[0],np.shape(distr_grid[0,:])[0]]
        
        # set up non-zero classes on array
        for i in range(arry_dim1):
            for j in range(arry_dim2):
                if(multivariate_normal.pdf([i,j],mean_fit_class,cov_mtrx) >= cut_off):
                    distr_grid[j,i] = 1
                    
        [genotypes,abundances] = generate_2D_discrete_gauss(distr_grid,min_fit_clss,mean_fit_clss,cov_matrx, N)
        [distr_grid,dim1_data,dim2_data] = get_2D_distr(genotypes,abundances,box_dim)
        
        # plot figure 1 with general with constructed discrete gaussian
        fig1, ax1 = plt.subplots()
        image = np.ones(np.shape(distr_grid)) - distr_grid
        
        ax1.imshow(image, cmap=plt.cm.gray, interpolation='nearest')
        ax1.set_title('dropped spines')
        
        #-------------------------------------------------------------------------
        # THIS WHOLE SECTION NEEDS TO BE FIXED, AND SHOULD LOOK LIKE FIG 2 SECTION
        #-------------------------------------------------------------------------
        # Move left and bottom spines outward by 10 points
        ax1.spines['left'].set_position(('outward', 10))
        ax1.spines['bottom'].set_position(('outward', 10))
        # Hide the right and top spines
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        # Only show ticks on the left and bottom spines
        ax1.yaxis.set_ticks_position('left')
        ax1.xaxis.set_ticks_position('bottom')
        
        fig1.savefig('./'+folder_location+'figures/'+fname+data_name+'.pdf')
        
        del min_fit_clss, mean_fit_clss, cov_mtrx, box_dim, cut_off
        del distr_grid, arry_dim1, arry_dim2, fig1, ax1, N, s1, s2, U1, U2

# figure 2: plot of sampled two dimensional distribution from simulated data        
    if(figNum == 2):
        
        # load time series data of distrStats from plotdata.py output
        pickle_file_name = './'+folder_location+'data/pythondata/timesGenosAbund'+data_name+'.pickle'
        pickle_file = open(pickle_file_name,'rb') 
        [times,genotypes,abundances] = pickle.load(pickle_file)
        pickle_file.close()
        
        # change the index to plot a different distribution
        snapshot_indx = get_sample_window(times,sim_start,sim_start)[0]
        genotypes = genotypes[snapshot_indx]
        abundances = abundances[snapshot_indx]    
        
        fit_clss_width = np.max([np.max(genotypes[:,0])-np.min(genotypes[:,0])+5,
                                 np.max(genotypes[:,1])-np.min(genotypes[:,1])+5])
                                 
        # use time index indx = 13630
        box_dim = [[fit_clss_width,2],[fit_clss_width,2]]

        [distr_grid,class_xlabels,class_ylabels] = get_2D_distr(genotypes,abundances,box_dim)

        distr_grid = np.log10(N*distr_grid)
        distr_grid[distr_grid == -inf] = 0
                    
        print genotypes
        # plot figure 1 with general with constructed discrete gaussian
        fig2, ax2 = plt.subplots(1,1,figsize=[10,8])

        fit_distr_2d = ax2.pcolormesh(distr_grid.transpose(),cmap=plt.cm.gray_r)
        cbar = plt.colorbar(fit_distr_2d)
        ax2.axis('tight')        
        ax2.set_xticks(np.arange(distr_grid.shape[1])+0.5)
        ax2.set_yticks(np.arange(distr_grid.shape[0])+0.5)        
        ax2.set_xticklabels(class_xlabels[1:],rotation=90)
        ax2.set_yticklabels(class_ylabels[1:])        
        ax2.set_xlabel('Beneficial Mutaitons Trait 1',fontsize=18,labelpad=20)
        ax2.set_ylabel('Beneficial Mutaitons Trait 2',fontsize=18,labelpad=10)
        ax2.tick_params(axis='both',labelsize=14)        
        cbar.ax.text(2.5,0.65,'Log$_{10}$ of Abundances',rotation=90,fontsize=18)

#        ax2.scatter([45],[33])
#        fit_line = plt.plot([i+40 for i in range(13)],[-i+38 for i in range(13)],ls="--", c=".3")

        fig2.savefig('./'+folder_location+'figures/'+fname+data_name+'.pdf')
        
        del times, genotypes, abundances, fit_clss_width, class_xlabels
        del class_ylabels, fit_distr_2d, cbar, fig2, ax2, N, s1, s2, U1, U2
    
# figure 3: plot of rate of adaptation, variances, covariance and their means
    
    if (figNum == 3):

        # load time series data of distrStats from plotdata.py output
        pickle_file_name = './'+folder_location+'data/pythondata/distrStats'+data_name+'.pickle'
        pickle_file = open(pickle_file_name,'rb') 
        [times,mean_fit,fit_var,fit_cov,pop_load,dcov_dt,vU_thry,v2U_thry] = pickle.load(pickle_file)
        pickle_file.close()
        
        # select interval of simulation data that will be used for plot
        # reduce loaded data to subset corresponding to selected interval
        [start_indx,end_indx] = get_sample_window(times,sim_start,sim_end)
        times = times[start_indx:end_indx]
        mean_fit = mean_fit[start_indx:end_indx]
        fit_var = fit_var[start_indx:end_indx]
        fit_cov = fit_cov[start_indx:end_indx]
        pop_load = pop_load[start_indx:end_indx]
        rate_adpt1 = fit_var[:,0] + fit_cov
        rate_adpt2 = fit_var[:,0] + fit_cov
        
        del start_indx, end_indx
        
        # compute means of the time series data and store as constant arry
        var1_avg = np.mean(fit_var[:,0])*np.ones(np.shape(times))
        var2_avg = np.mean(fit_var[:,1])*np.ones(np.shape(times))
        cov_avg = np.mean(fit_cov[:])*np.ones(np.shape(times))
        rate_adpt1_avg = var1_avg + cov_avg
        rate_adpt2_avg = var2_avg + cov_avg
        
        # plot data for figure 2a and 2b 
        fig3a, ax3a = plt.subplots(1,1,figsize=[8,8])
        ax3b = plt.twinx(ax3a)
        ax3a.plot(times,var1_avg,c="black",label='var($r_1|r_2$)=' + str(round(var1_avg[0],7)),linewidth=2.0,linestyle = '--')
        ax3a.plot(times,cov_avg,c="black",label='cov($r_1$,$r_2$)=' + str(round(cov_avg[0],7)),linewidth=2.0,linestyle = '-.')
        ax3a.plot(times,vU_thry*np.ones(np.shape(var1_avg)),c="black",label='var($r_1$)=' + str(round(vU_thry,7)),linewidth=2.0,linestyle = '-')        
        ax3a.axhline(linewidth=0.5, color = 'k')
        ax3a.set_ylabel('Fitness Variances & Covariance',fontsize=18)
        ax3a.set_xlabel('Time (Generations)',fontsize=18)
        ax3a.set_ylim((-3e-4,4e-4))
        ax3a.set_xlim((1e4,2e4))
        ax3a.tick_params(axis='both',labelsize=14)
        ax3a.ticklabel_format(style='sci',axis='both',scilimits=(0,0))
        ax3a.legend()
        
        ax3b.plot(times,fit_var[:,0],c="black",linestyle="--",linewidth=1.0)
        ax3b.plot(times,fit_cov[:],c="black",linestyle="-.",linewidth=1.0)
        ax3b.axhline(linewidth=0.5, color = 'k')        
        ax3b.set_ylim((-3e-4,4e-4))
        ax3b.set_xlim((1e4,2e4))        
        ax3b.ticklabel_format(style='sci',axis='both',scilimits=(0,0))
        ax3b.tick_params(axis='both',labelsize=14)

        fig3a.savefig('./'+folder_location+'figures/'+fname+data_name+'.pdf')
        
        del times, mean_fit, fit_var, fit_cov, pop_load, dcov_dt, vU_thry
        del v2U_thry, rate_adpt1, rate_adpt2, var1_avg, var2_avg,cov_avg
        del rate_adpt1_avg, rate_adpt2_avg, fig3a, ax3a, ax3b, N, s1, s2, U1, U2

# figure 4: plot of trait or fitness distribution versus normal distribution
    
    if (figNum == 4):

        # load time series data of distrStats from plotdata.py output
        pickle_file_name = './'+folder_location+'data/pythondata/timesGenosAbund'+data_name+'.pickle'
        pickle_file = open(pickle_file_name,'rb') 
        [times,genotypes,abundances] = pickle.load(pickle_file)
        pickle_file.close()
        
        # change the index to plot a different distribution
        snapshot_indx = get_sample_window(times,sim_start,sim_start)[0]
        genotypes = genotypes[snapshot_indx]
        abundances = abundances[snapshot_indx]   
        
        [mu, var] = get_trait_mean_var(genotypes,abundances,traitno)
        [trait_bins,trait_tot] = get_1D_proj(genotypes,abundances,traitno)
                
        # plot data for figure 2a and 2b 
        fig4, ax4 = plt.subplots(1,1,figsize=[8,8])
        ax4.hist(trait_bins, len(trait_bins), normed=1, 
                                    facecolor='green', alpha=0.75, weights=trait_tot)
        x = np.asarray([np.min(trait_bins)+(np.max(trait_bins)-np.min(trait_bins))
                                    *i/100.0 for i in range(101)])
        y = mlab.normpdf(x, mu, np.sqrt(var))
        ax4.plot(x, y, 'r--', linewidth=1)
        ax4.set_title('1D distribution trait '+str(traitno))
        fig4.savefig('./'+folder_location+'figures/'+fname+data_name+'.png')

# figure 5: plots of var-cov means wrt 
    if (figNum == 5):
        
        # load time series data of distrStats from plotdata.py output
        pickle_file_name = './'+folder_location+'data/pythondata/sumdata_exp5.pickle'
        pickle_file = open(pickle_file_name,'rb') 
        [var, cov, vUthry, v2Uthry, varp, covp, vUthryp, v2Uthryp,NsUparam] = pickle.load(pickle_file)
        pickle_file.close()
        
        [start1,start2,start3] = [1,num_exp/3+1,2*num_exp/3+1]         
        [end1,end2,end3] = [num_exp/3,2*num_exp/3,num_exp]
        NsUparam = np.asarray(NsUparam)        
        
        fig5a, ax5a = plt.subplots(1,1,figsize=[16,8])
        ax5a.plot(NsUparam[start1:end1,1],-np.asarray(covp[start1:end1])-np.asarray(varp[start1:end1])+np.asarray(vUthryp[start1:end1]),c="black",label='$v_1$',linewidth=2.0,linestyle = '-')                
        ax5a.plot(NsUparam[start1:end1,1],np.asarray(varp[start1:end1]),c="black",label='$\sigma_{1}^2$',linewidth=2.0,linestyle = ':')
        ax5a.plot(NsUparam[start1:end1,1],np.asarray(covp[start1:end1]),c="black",label='$\sigma_{12}$',linewidth=2.0,linestyle = '-.')        
        ax5a.axhline(linewidth=0.5, color = 'k')
        ax5a.set_ylabel('Scaled Variance-Covariance',fontsize=18)
        ax5a.set_xlabel('selection coeff ($10^{-3} - 2 x 10^{-2}$)',fontsize=18)
        ax5a.tick_params(axis='both',labelsize=14)
        ax5a.ticklabel_format(style='sci',axis='x',scilimits=(0,0))
        ax5a.set_ylim((-2,2))
        ax5a.set_xlim((min(NsUparam[start1:end1,1]),max(NsUparam[start1:end1,1])))  
        ax5a.legend()
        fig5a.savefig('./'+folder_location+'figures/'+fname+data_name+'d.pdf')
        
        fig5b, ax5b = plt.subplots(1,1,figsize=[16,8])
        ax5b.plot(NsUparam[start2:end2,2],-np.asarray(covp[start2:end2])-np.asarray(varp[start2:end2])+np.asarray(vUthryp[start2:end2]),c="black",label='$v_1$',linewidth=2.0,linestyle = '-')
        ax5b.plot(NsUparam[start2:end2,2],np.asarray(varp[start2:end2]),c="black",label='$\sigma_{1}^2$',linewidth=2.0,linestyle = ':')
        ax5b.plot(NsUparam[start2:end2,2],np.asarray(covp[start2:end2]),c="black",label='$\sigma_{12}$',linewidth=2.0,linestyle = '-.')
        ax5b.axhline(linewidth=0.5, color = 'k')
        ax5b.set_ylabel('Scaled Variance-Covariance',fontsize=18)
        ax5b.set_xlabel('mutation rate ($10^{-6} - 10^{-4}$)',fontsize=18)
        ax5b.tick_params(axis='both',labelsize=14)
        ax5b.ticklabel_format(style='sci',axis='x',scilimits=(0,0))
        ax5b.set_ylim((-2,2))
        ax5b.set_xlim((min(NsUparam[start2:end2,2]),max(NsUparam[start2:end2,2])))          
        ax5b.legend()
        fig5b.savefig('./'+folder_location+'figures/'+fname+data_name+'e.pdf')
        
        fig5c, ax5c = plt.subplots(1,1,figsize=[16,8])
        ax5c.plot(NsUparam[start3:end3,0],-np.asarray(covp[start3:end3])-np.asarray(varp[start3:end3])+np.asarray(vUthryp[start3:end3]),c="black",label='$v_1$',linewidth=2.0,linestyle = '-')        
        ax5c.plot(NsUparam[start3:end3,0],np.asarray(varp[start3:end3]),c="black",label='$\sigma_{1}^2$',linewidth=2.0,linestyle = ':')
        ax5c.plot(NsUparam[start3:end3,0],np.asarray(covp[start3:end3]),c="black",label='$\sigma_{12}$',linewidth=2.0,linestyle = '-.')                
        ax5c.axhline(linewidth=0.5, color = 'k')
        ax5c.set_ylabel('Scaled Variance-Covariance',fontsize=18)
        ax5c.set_xlabel('population size ($10^{6} - 10^{9}$)',fontsize=18)
        ax5c.tick_params(axis='both',labelsize=14)
        ax5c.ticklabel_format(style='sci',axis='x',scilimits=(0,0))
        ax5c.set_ylim((-2,2))
        ax5c.set_xlim((min(NsUparam[start3:end3,0]),max(NsUparam[start3:end3,0]))) 
        ax5c.legend()        
        fig5c.savefig('./'+folder_location+'figures/'+fname+data_name+'f.pdf')

    return None