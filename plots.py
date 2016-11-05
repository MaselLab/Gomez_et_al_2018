# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/home/dirge3141/.spyder2/.temp.py
"""

# The main purpose of this code is to get summary statistics for the 
# simulation results as functions of the ratios of key parameters

# IMPORT PACKAGES AND DEFINE LOCAL FUNCTIONS----------------------------------
#-----------------------------------------------------------------------------
import pickle
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

#reference of uploaded data types:
getp=dict(Ur=0,Ua=1,sr=2,sa=3,K=4,b=5,d=6,alldata=7)

# the first thing to do is the get some visualizations and store them somwhere

# GETTING PROCESSED DATA - can only due one data set at a time due to size
[amin, amax, aavg, times, data,select_data]=[[],[],[],[],[],[]]

for expNo in range(30):
    for i in range(16):
        fd = open('./AllData/PythonData/simdata'+str(i)+'.data')
        data = pickle.load(fd)
        select_data=[data[getp['Ur']],data[getp['Ua']],data[getp['sr']],data[getp['sa']], \
            data[getp['alldata']][expNo]]
        amin=[24]+[val[1] for val in select_data[4]]
        amax=[24]+[val[2] for val in select_data[4]]
        aavg=[24]+[val[3] for val in select_data[4]]
        times=[0]+[val[0] for val in select_data[4]]
        fig1,axes1=plt.subplots(1,1,figsize=[8,8])
        subfig1=plt.subplot(111)
        plt.plot(times,amin,c='blue')
        plt.plot(times,amax,c='blue')
        plt.plot(times,aavg,c='red')
        plt.ylim([0,60])
        plt.xlim([0,90000])    
        plt.title("Sim: "+str(i)+", Exp: "+str(expNo+1)+", Ur: "+str(select_data[0]) \
                +", Ua: "+str(select_data[1])+", sr: "+str(select_data[2]) \
                +", sa: "+str(select_data[3]))
        plt.xlabel("generations")
        plt.ylabel("absolute fitness")
        plt.show()
        fig1.savefig('./images/RelAbsImgs/sim'+str(i)+"exp"+str(expNo+1))

#Before we plot the means of the curves, it will be worth it to check out the
#of some of these changes

fig1,axes1=plt.subplots(1,1,figsize=[8,8])
subfig=plt.subplot(111)
for i in range(len(data_slices)-2):
    plt.plot(data_slices[i+2][1],data_slices[i+2][2],label="UrUa="+str(data_slices[i+2][0]))
#plt.plot(simgen1[:,0],simgen1[:,2],c="blue")
#plt.plot(simgen1[:,0],simgen1[:,3],c="red")
#plt.plot(simgen1[:,0],mvgavg1,c="black",linewidth=2)
#plt.xlabel("generations")
#plt.ylabel("abs fitness class")
#plt.subplot(312)
#plt.plot(simgen1[:,0],simgen1[:,4],c="blue")
#plt.plot(simgen1[:,0],simgen1[:,5],c="blue")
#plt.plot(simgen1[:,0],simgen1[:,6],c="red")
#plt.xlabel("generations")
#plt.ylabel("rel fitness class")
#plt.subplot(313)
#plt.plot(simgen1[:,0],simgen1[:,7],c="red")
plt.legend(loc=4)
plt.xlabel("Ratio of Rel. & Abs. Fitness Increments")
plt.ylabel("Mean Fitness Class")
plt.show()
fig1.savefig('./UrUaSrSaCurves')

# CALCULATE  THE GEOMETRIC MEANS OF THE DATA
# get all the base data for model with no competitive mutations
#------------------------------------------------------------------------------


# BURN A SET OF THE TRANSITIONS DATA
#------------------------------------------------------------------------------


#-------------------- The plots of the abs fit stoch chain --------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

fig2,axes2=plt.subplots(1,1,figsize=[8,8])
subfig2=plt.subplot(111)
for i in range(len(data_slices)-2):
    plt.plot(data_slices[i+2][1],data_slices[i+2][2],label="UrUa="+str(data_slices[i+2][0]))
#plt.plot(simgen1[:,0],simgen1[:,2],c="blue")
#plt.plot(simgen1[:,0],simgen1[:,3],c="red")
#plt.plot(simgen1[:,0],mvgavg1,c="black",linewidth=2)
#plt.xlabel("generations")
#plt.ylabel("abs fitness class")
#plt.subplot(312)
#plt.plot(simgen1[:,0],simgen1[:,4],c="blue")
#plt.plot(simgen1[:,0],simgen1[:,5],c="blue")
#plt.plot(simgen1[:,0],simgen1[:,6],c="red")
#plt.xlabel("generations")
#plt.ylabel("rel fitness class")
#plt.subplot(313)
#plt.plot(simgen1[:,0],simgen1[:,7],c="red")
plt.legend(loc=4)
plt.xlabel("Ratio of Rel. & Abs. Fitness Increments")
plt.ylabel("Mean Fitness Class")
plt.show()
fig2.savefig('./UrUaSrSaCurves')