# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 19:56:37 2021

@author: Admin
"""


# https://github.com/LiviaJanke/SummerProject.git
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
import math
#%%
t_gas,tt,co2,ch4,n2o=np.loadtxt('C:/Users/Admin/Documents/SummerProject/Data/cleandata.csv',skiprows=1,delimiter=',',unpack=True)
years,temp_no,temp=np.loadtxt('C:/Users/Admin/Documents/SummerProject/Data/graph.txt',skiprows=5,unpack=True)
time=years[1:]

surface=510e12
T_lw=0.2 #transmittance of the long wavelength of the atmosphere
B=2.218 #Wm^-2K^-1 constant for OLR, determined empirically
alpha_co2=5.35
ach4=0.036
an2o=0.12
#%% 70% water, considering depths of 70m for well-mixed increases, 1000kgm^-3, heat cap 
 
depth_water=70
depth_soil=1.6

def calc_mc_tot(depth_water,depth_soil):
    c_water=4180 #Jkg^-1C-1
    av_c_soil=1000
    density_water=1000
    density_soil=1300
    mc_water=0.7*surface*depth_water*density_water*c_water
    mc_soil=0.3*surface*depth_soil*density_soil*av_c_soil
    mc_tot=(mc_water+mc_soil)/3.1536e7 #in units of years
    return mc_tot
mc_tot=calc_mc_tot(70,1.6)

#%%
def forcing_CO2(data):
    avco2=data
    dfco2=[]
    for i in range(0,len(time)):
        beta=avco2[i+1]/avco2[0]
        F=alpha_co2*np.log(beta)
        dfco2.append(F)
    
    return dfco2


def f(M,N):
    f_function=0.47*np.log(1+2.01e-5*(M*N)**0.75 +5.31e-15*M*(M*N)**1.52)
    return f_function

def forcing_methane(meth_data,n2o_data):
    dfch4_list=[]
    for i in range(len(meth_data)-1):
        dfch4 = 0
        dfch4 += (meth_data[i+1])**0.5
        dfch4 -= (meth_data[0])**0.5
        dfch4 *= ach4
        dfch4 -= f(meth_data[i+1],n2o_data[i+1])-f(meth_data[0],n2o_data[0])
        dfch4_list.append(dfch4)         
    return dfch4_list
      
def forcing_n2o(meth_data,n2o_data):
    dfn2o_list=[]
    for i in range(len(meth_data)-1):
        dfn2o=0
        dfn2o+=(n2o_data[i+1])**0.5-(n2o_data[0]**0.5)
        dfn2o*=0.12
        dfn2o-=f(meth_data[0],n2o_data[i+1])-f(meth_data[0],n2o_data[0])
        dfn2o_list.append(dfn2o)
    return dfn2o_list

## by taking base temperature as 287K, and thus 1958 temperature as 287K+anomaly
# according to penas report, the relationship between olr and temp is as follows: A+B*T where A=-339.647 and B=2.218 

def temp_increase(number_years):
    equilibrium_temperature=287
    T=286.1
    anomaly=-0.09
    temperature=[]
    increase_temp=-0.09
    excess_planetary_energy=[]
    dOLR=B*anomaly 
    dF_CO2=forcing_CO2(co2)
    dF_N2O=forcing_n2o(ch4,n2o)
    dF_methane=forcing_methane(ch4,n2o)
    for i in range (0,number_years):
        dF_tot=+dF_CO2[i]+dF_N2O[i]+dF_methane[i]
        excess_planetary_energy=(dF_tot-dOLR)*surface
        dT=excess_planetary_energy/mc_tot
        anomaly+=dT
        dOLR=B*anomaly
        increase_temp+=dT
        temperature.append(increase_temp)
    return temperature

temperature=temp_increase(len(time))
plt.plot(years,temp)
plt.plot(time,temperature)

plt.show()

print(temperature[len(time)-1])

