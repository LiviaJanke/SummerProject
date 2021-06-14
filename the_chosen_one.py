# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 14:48:39 2021

@author: louis
"""

# https://github.com/LiviaJanke/SummerProject.git
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
import math
import seaborn as sns
#%%
t_gas,tt,co2,ch4,n2o=np.loadtxt('Data/cleandata.csv',skiprows=1,delimiter=',',unpack=True)
years,temp_no,temp=np.loadtxt('Data/graph.txt',skiprows=5,unpack=True)
time=years[1:]
print(len(time))
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
    for i in range(0,len(data)-1):
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
        dfn2o*=an2o
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
        dF_tot=dF_CO2[i]+dF_N2O[i]+dF_methane[i]
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
print(time[0])
#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op
import sympy as sym
import csv
#%% projection to the future
time2, n, co2_conc, ch4_conc, n2o_conc = np.loadtxt('Data/cleandata.csv', skiprows = 81, delimiter = ',', unpack = 1)




def exp_projection(time: list, conc:list, gas_name_unit: str, int_guess: list):
    x = np.linspace(2020, 2070, 51)
    def exp_func(time2, *int_guess):
        return int_guess[0]*np.e**(int_guess[1]*((time2-1960)**int_guess[2]))                                           #plots twice a month
    params, pcov = op.curve_fit(exp_func, time2, conc, int_guess)  
    plt.figure(figsize = (10, 7))
    plt.errorbar(x[::24], exp_func((x[::24]), *params), yerr = 0, fmt = 'y+', mew=2, ms=3, capsize = 2)
    plt.grid()
    plt.xlabel('Year')
    plt.ylabel('Concentration of: {}'.format(gas_name_unit))
    plt.title('Extrapolated Concentration over time of: {}'.format(gas_name_unit), fontsize = 15)
    plt.plot(x, exp_func(x, *params), 'r', linewidth = 2)

    plt.figure(figsize = (10, 7))
    plt.errorbar(time2, conc, yerr = 0, fmt = 'r+', mew=2, ms=2, capsize = 2)
    plt.grid()
    plt.xlabel('Year')
    plt.ylabel('Concentration of: {}'.format(gas_name_unit))
    plt.title('Concentration over time of: {}'.format(gas_name_unit), fontsize = 15)
    plt.plot(time2, exp_func(time2, *params), 'b', linewidth = 3)
    print("Exponential Projection Implies an increase, per year, by a factor of approx: {}".format(np.e**params[1])) 
    print(params, ' = parameters' , pcov, ' = covariance')
    
    proj_data = [x, exp_func(x, *params)]
    data = np.array([time, exp_func(time, *params)])
    return proj_data, data, params

co2_exp_proj, co2_fit_data, co2_params = exp_projection(time, co2_conc, 'CO2 (ppm)', int_guess = [2.465e-3, 5.99e-3, 1])

n2o_exp_proj, n2o_fit_data, n2o_params = exp_projection(time, n2o_conc, 'N2O (ppb)', [0, 0.005, 1])

ch4_exp_proj, ch4_fit_data, ch4_params = exp_projection(time, ch4_conc, 'CH4 (ppb)', [0, 0.005, 1])
#%%
#print(co2_exp_proj)
#print(ch4_exp_proj)
#print(n2o_exp_proj[0])

#%% projection using only ch4 n2o co2 and current trends
t_gas,tt,co2,ch4,n2o=np.loadtxt('Data/cleandata.csv',skiprows=1,delimiter=',',unpack=True)
#print(len(co2))
#mc_tot_mean=calc_mc_tot(70,1.6)
#mc_tot_min=calc_mc_tot(40,1.6)
#mc_tot_max=calc_mc_tot(100,1.6)
def temp_increase(number_years,data_co2,data_n2o,data_ch4):
    equilibrium_temperature=287
    T=286.1
    anomaly=-0.09
    temperature=[]
    increase_temp=-0.09
    excess_planetary_energy=[]
    dOLR=B*anomaly 
    dF_CO2=forcing_CO2(data_co2)
    print(len(dF_CO2))
    dF_N2O=forcing_n2o(data_ch4,data_n2o)
    print(len(dF_N2O))
    dF_methane=forcing_methane(data_ch4,data_n2o)
    print(len(dF_methane))
    for i in range (0,number_years):
        dF_tot=dF_CO2[i]+dF_N2O[i]+dF_methane[i]
        excess_planetary_energy=(dF_tot-dOLR)*surface
        dT=excess_planetary_energy/mc_tot
        anomaly+=dT
        dOLR=B*anomaly
        increase_temp+=dT
        temperature.append(increase_temp)
    return temperature

co2_proj=co2_exp_proj[1]
n2o_proj=n2o_exp_proj[1]
ch4_proj=ch4_exp_proj[1]

new_co2=np.concatenate((co2,co2_proj[1:]))
new_n2o=np.concatenate((n2o,n2o_proj[1:]))
new_ch4=np.concatenate((ch4,ch4_proj[1:]))
new_time=np.linspace(1881,2070,190)
#print(len(ch4))
#print(co2)
#print(len(co2))
#print(time3)
#print(len(time3))


temperature_projection=temp_increase(len(new_time),new_co2,new_n2o,new_ch4)

plt.plot(years,temp)
plt.plot(new_time,temperature_projection)

plt.show()
print(new_time[139])
print(new_time[len(new_time)-1])
print(temperature_projection[len(new_time)-1])
print(temperature_projection[139])
print(temperature_projection[len(new_time)-1]-temperature_projection[139])

#%% fitting including CFC 
alpha_CFC11= 0.25
alpha_CFC12 = 0.32
CFC11_years, CFC11x = np.loadtxt('Data/CFC11_1880_to_present.csv', skiprows = 1, delimiter = ',', unpack = True)
CFC12_years, CFC12x = np.loadtxt('Data/CFC12_1880_to_2021_means.csv', skiprows = 1, delimiter = ',', unpack = True)
CFC11=CFC11x[0:142]
CFC12=CFC12x[0:142]
def forcing_CFC11(CFC11_data):
    dfCFC11 =[]
    for i in range(0,len(CFC11_data)-1): # -2 because don't want 2021 
        F = alpha_CFC11 *( CFC11_data[i+1]-CFC11_data[0])
        dfCFC11.append(F)
    return dfCFC11

def forcing_CFC12(CFC12_data):
    dfCFC12 =[]
    for i in range(0,len(CFC12_data)-1):
        F = alpha_CFC12 * (CFC12_data[i+1]-CFC12_data[0])
        dfCFC12.append(F)
    return dfCFC12

def temp_increase_new(number_years,data_co2,data_n2o,data_ch4,data_cfc11,data_cfc12):
    equilibrium_temperature=287
    T=286.1
    anomaly=-0.09
    temperature=[]
    increase_temp=-0.09
    excess_planetary_energy=[]
    dOLR=B*anomaly 
    dF_CO2=forcing_CO2(data_co2)
    dF_N2O=forcing_n2o(data_ch4,data_n2o)
    dF_methane=forcing_methane(data_ch4,data_n2o)
    dF_CFC11 = forcing_CFC11(data_cfc11)
    dF_CFC12 = forcing_CFC12(data_cfc12)
    for i in range (0,number_years):
        dF_tot=dF_CO2[i]+dF_N2O[i]+dF_methane[i] + dF_CFC11[i] + dF_CFC12[i]
        excess_planetary_energy=(dF_tot-dOLR)*surface
        dT=excess_planetary_energy/mc_tot
        anomaly+=dT
        dOLR=B*anomaly
        increase_temp+=dT
        temperature.append(increase_temp)
    return temperature

temperature_old=temp_increase(len(time),co2,n2o,ch4)
temperature_new=temp_increase_new(len(time),co2,n2o,ch4,CFC11,CFC12)
plt.plot(years,temp)
plt.xlabel('Years')
plt.ylabel('Temperature Anomaly')
plt.plot(time,temperature_old, color = 'red', label = 'old model')
plt.plot(time,temperature_new, color = 'green', label = 'new model')
plt.legend()

plt.show()

print(temperature_new[len(time)-1])


#%% fitting with everything so volcanoes
t_gas,tt,co2,ch4,n2o=np.loadtxt('Data/cleandata.csv',skiprows=1,delimiter=',',unpack=True)
years,temp_no,temp=np.loadtxt('Data/graph.txt',skiprows=5,unpack=True)
sulph_years, sulphates = np.loadtxt('Data/sulphate_annual_medians_from1880_v2.csv', skiprows = 1, delimiter = ',', unpack = True)
CFC11_years, CFC11 = np.loadtxt('Data/CFC11_1880_to_present.csv', skiprows = 1, delimiter = ',', unpack = True)
CFC12_years, CFC12 = np.loadtxt('Data/CFC12_1880_to_2021_means.csv', skiprows = 1, delimiter = ',', unpack = True)
volcanic_years, volcanic_forcing = np.loadtxt('Data/volcanicforcingdata.csv', skiprows = 1, delimiter = ',', unpack = True)

aerosols=-0.5
#aerosols=np.linspace(0,-0.5,141,endpoint=True)
def temp_increase_newest(number_years,data_co2,data_n2o,data_ch4,data_cfc11,data_cfc12,data_volc):
    equilibrium_temperature=287
    T=286.1
    anomaly=-0.09
    temperature=[]
    increase_temp=-0.09
    excess_planetary_energy=[]
    dOLR=B*anomaly 
    dF_CO2=forcing_CO2(data_co2)
    dF_N2O=forcing_n2o(ch4,n2o)
    dF_methane=forcing_methane(data_ch4,data_n2o)
    dF_CFC11 = forcing_CFC11(data_cfc11)
    dF_CFC12 = forcing_CFC12(data_cfc12)
    dF_volcano = data_volc
    for i in range (0,number_years):
#        dF_tot=+dF_CO2[i]+dF_N2O[i]+dF_methane[i] + dF_sulph[i] + dF_CFC11[i] + dF_CFC12[i]
        dF_tot=+dF_CO2[i]+dF_N2O[i]+dF_methane[i] + dF_CFC11[i] + dF_CFC12[i] + dF_volcano[i]+aerosols
        excess_planetary_energy=(dF_tot-dOLR)*surface
        dT=excess_planetary_energy/mc_tot
        anomaly+=dT
        dOLR=B*anomaly
        increase_temp+=dT
        temperature.append(increase_temp)
    return temperature


temperature_newest=temp_increase_newest(len(time),co2,n2o,ch4,CFC11,CFC12,volcanic_forcing)
plt.plot(years,temp)
plt.xlabel('Years')
plt.ylabel('Temperature Anomaly')

plt.plot(time,temperature_newest, color = 'green', label = 'new model')
plt.legend()

plt.show()

print(temperature_newest[len(time)-1])


#%% fitting with cfc 
def forcing_CFC11(CFC11_data):
    dfCFC11 =[]
    for i in range(0,len(CFC11_data)-1): 
        F = alpha_CFC11 *( CFC11_data[i+1]-CFC11_data[0])
        dfCFC11.append(F)
    return dfCFC11

def forcing_CFC12(CFC12_data):
    dfCFC12 =[]
    for i in range(0,len(CFC12_data)-1):
        F = alpha_CFC12 * (CFC12_data[i+1]-CFC12_data[0])
        dfCFC12.append(F)
    return dfCFC12

t_gas,tt,co2,ch4,n2o=np.loadtxt('Data/cleandata.csv',skiprows=1,delimiter=',',unpack=True)
CFC11_years, CFC11 = np.loadtxt('Data/CFC11_1880_to_present.csv', skiprows = 1, delimiter = ',', unpack = True)
CFC12_years, CFC12 = np.loadtxt('Data/CFC12_1880_to_2021_means.csv', skiprows = 1, delimiter = ',', unpack = True)
time_data,co2_proj,n2o_proj,ch4_proj,cfc11_proj,cfc12_proj=np.loadtxt('Projections/All_Projections.csv',skiprows=1,delimiter=',',unpack=True)

timee,cfc11_proj,cfc12_proj=np.loadtxt("Data/cfc_projections.csv",skiprows=1,delimiter=',',unpack=True)



new_cfc11=np.concatenate((CFC11[0:141],cfc11_proj[1:]))
new_cfc12=np.concatenate((CFC12[0:141],cfc12_proj[1:]))
print(len(new_cfc11))

#plt.plot(new_time,new_cfc12[1:],'x')
#plt.plot(new_time,new_cfc12[1:])
#plt.plot(new_time,new_cfc11[1:])
#plt.show()
temperature_projection=temp_increase_new(len(new_time),new_co2,new_n2o,new_ch4,new_cfc11,new_cfc12)

plt.plot(years,temp)
plt.plot(new_time,temperature_projection)

plt.show()
#print(new_time[139])
#print(new_time[len(new_time)-1])
#print(temperature_projection[len(new_time)-1])
#print(temperature_projection[139])
print(temperature_projection[len(new_time)-1]-temperature_projection[139])


#%% target projections
target_time,target_co2,target_n2o,target_ch4=np.loadtxt('Projections/Target_Projections.csv',skiprows=1,delimiter=',',unpack=True)


tar_co2=np.concatenate((co2,target_co2[1:]))
tar_n2o=np.concatenate((n2o,target_n2o[1:]))
tar_ch4=np.concatenate((ch4,target_ch4[1:]))
tar_time=np.linspace(1881,2070,190)
print(len(tar_co2))
target_projection=temp_increase(len(tar_time),tar_co2,tar_n2o,tar_ch4)

plt.plot(years,temp)
plt.plot(tar_time,target_projection)

plt.show()
#print(new_time[139])
#print(new_time[len(new_time)-1])
#print(temperature_projection[len(new_time)-1])
#print(temperature_projection[139])
print(target_projection[len(tar_time)-1]-target_projection[139])

#%% Figuring out aerosols

aero=-0.5 #np.linspace(0,-0.5,141,endpoint=True)
aero_min=-0.26 #np.linspace(0,-0.26,141,endpoint=True)
aero_max=-0.82 #np.linspace(0,-0.82,141,endpoint=True)
def temp_increase_corrected(number_years,data_co2,data_n2o,data_ch4,data_cfc11,data_cfc12,aerosols):
    equilibrium_temperature=287
    T=286.1
    anomaly=-0.09
    temperature=[]
    increase_temp=-0.09
    excess_planetary_energy=[]
    dOLR=B*anomaly 
    dF_CO2=forcing_CO2(data_co2)
    dF_N2O=forcing_n2o(data_ch4,data_n2o)
    dF_methane=forcing_methane(data_ch4,data_n2o)
    dF_CFC11 = forcing_CFC11(data_cfc11)
    dF_CFC12 = forcing_CFC12(data_cfc12)
    for i in range (0,number_years):
        dF_tot=dF_CO2[i]+dF_N2O[i]+dF_methane[i] + dF_CFC11[i] + dF_CFC12[i]+aerosols
        excess_planetary_energy=(dF_tot-dOLR)*surface
        dT=excess_planetary_energy/mc_tot
        anomaly+=dT
        dOLR=B*anomaly
        increase_temp+=dT
        temperature.append(increase_temp)
    return temperature

temperature_aero=temp_increase_corrected(len(time),co2,n2o,ch4,CFC11,CFC12,aero)
temp_max_aero=temp_increase_corrected(len(time),co2,n2o,ch4,CFC11,CFC12,aero_min)
temp_min_aero=temp_increase_corrected(len(time),co2,n2o,ch4,CFC11,CFC12,aero_max)
temperature_new=temp_increase_new(len(time),co2,n2o,ch4,CFC11,CFC12)
plt.plot(years,temp)
plt.fill_between(time,temp_min_aero,temp_max_aero,facecolor='pink')
plt.xlabel('Years')
plt.ylabel('Temperature Anomaly')
plt.plot(time,temperature_aero, color = 'red', label = 'model with aeorosols')
plt.plot(time,temperature_new, color = 'green', label = 'model without aeorsols')
plt.legend()

plt.show()

print(temperature_new[len(time)-1])
print(temperature_aero[len(time)-1])
array_aero=np.array(temperature_aero)
array_new=np.array(temperature_new)

#residuals=array_aero-array_new
#plt.plot(time,residuals)
#plt.show()
#%% add volcanoes to aerosols
#this is done in the cell with volcanoes, just modify the forcing for aerosols

#%%
Year,WMGHG,Ozone,Solar,Land_Use,SnowAlb_BC,Orbital,TropAerDir,TropAerInd,StratAer=np.loadtxt('Data/instant_forcings_1880_to_2020.csv',skiprows=1,delimiter=',',unpack=True)

def temp_increase_with_all(number_years,data_co2,data_n2o,data_ch4,data_cfc11,data_cfc12,ozone,solar,land_use,snow,orbital,aer_dir,aer_indir,strat_aer):
    equilibrium_temperature=287
    T=286.1
    anomaly=-0.09
    temperature=[]
    increase_temp=-0.09
    excess_planetary_energy=[]
    dOLR=B*anomaly 
    dF_CO2=forcing_CO2(data_co2)
    dF_N2O=forcing_n2o(data_ch4,data_n2o)
    dF_methane=forcing_methane(data_ch4,data_n2o)
    dF_CFC11 = forcing_CFC11(data_cfc11)
    dF_CFC12 = forcing_CFC12(data_cfc12)
    for i in range (0,number_years):
        dF_tot=dF_CO2[i]+dF_N2O[i]+dF_methane[i] + dF_CFC11[i] + dF_CFC12[i]+ozone[i]+solar[i]+land_use[i]+snow[i]+orbital[i]+aer_dir[i]+aer_indir[i]+strat_aer[i]
        excess_planetary_energy=(dF_tot-dOLR)*surface
        dT=excess_planetary_energy/mc_tot
        anomaly+=dT
        dOLR=B*anomaly
        increase_temp+=dT
        temperature.append(increase_temp)
    return temperature

temperature_old=temp_increase(len(time),co2,n2o,ch4)
temperature_new=temp_increase_with_all(len(time),co2,n2o,ch4,CFC11,CFC12,Ozone,Solar,Land_Use,SnowAlb_BC,Orbital,TropAerDir,TropAerInd,StratAer)
temp_min_all=temp_increase_with_all(len(time),co2,n2o,ch4,CFC11,CFC12,Ozone,Solar,Land_Use,SnowAlb_BC,Orbital,TropAerDir+0.5,TropAerInd+0.5,StratAer)
temp_max_all=temp_increase_with_all(len(time),co2,n2o,ch4,CFC11,CFC12,Ozone,Solar,Land_Use,SnowAlb_BC,Orbital,TropAerDir-0.5,TropAerInd-0.5,StratAer)
plt.fill_between(time,temp_min_all,temp_max_all,facecolor='pink')
plt.plot(years,temp)
plt.xlabel('Years')
plt.ylabel('Temperature Anomaly')
plt.plot(time,temperature_old, color = 'red', label = 'old model')
plt.plot(time,temperature_new, color = 'green', label = 'new model')
plt.legend()

plt.show()

print(temperature_new[len(time)-1])



