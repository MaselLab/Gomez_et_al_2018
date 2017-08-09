# -*- coding: utf-8 -*-
"""
Created on Tue Aug 08 18:08:03 2017
@author: Kevin Gomez (Masel Lab)
Library of functions used in plots.py and plotdata.py
"""

# import packages for script
import scipy as sp
import numpy as np
from scipy.stats import multivariate_normal

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
# returns: an array whose elements are the abundances of the fit classes

    tot_pop_size = sum(abundances)
    
    dim1_data = [np.min([genotypes[i][0] for i in range(len(abundances))]),
                 np.max([genotypes[i][0] for i in range(len(abundances))])]
    dim2_data = [np.min([genotypes[i][1] for i in range(len(abundances))]),
                 np.max([genotypes[i][1] for i in range(len(abundances))])]
                 
    my_distr = np.zeros([box_dim[0][0],box_dim[1][0]])
        
    for i in range(len(abundances)):
        indx1 = genotypes[i][0] - dim1_data[1] + box_dim[0][1]
        indx2 = genotypes[i][1] - dim2_data[1] + box_dim[1][1]
        my_distr[indx1,indx2] = abundances[i]/tot_pop_size

    return [np.flipud(my_distr.transpose()),dim1_data,dim2_data]