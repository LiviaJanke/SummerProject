# -*- coding: utf-8 -*-
"""
Created on Wed May 19 16:07:22 2021

@author: louis
"""

#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
import math
#%% plotting CO2 against time
#data=pd.read_csv('data/co2_mm_mlo.txt',delimiter=',')
#print(data.head())
data2=np.loadtxt('Data/co2_mm_mlo.txt') #loading CO2 against years data
print(data2)
plt.plot(data2[:,2],data2[:,3])
plt.show()

#sinusoidal part is probably due to seasonal fluctuations and can be ignored

#%% plotting temp increase against time
year1,temp1,temp2=sp.loadtxt('Data/graph.txt',skiprows=5,unpack=True)

plt.plot(year1[78:],temp2[78:])
print
#%% setting 1958 as temperature 0 and working out increase in temp from there
temp=temp2[78:]-temp2[78]
time=year1[79:]
plt.plot(year1[78:],temp)
print(np.max(temp))
print(year1[78:])
print(len(time))
#%% calculating beta
co2_data=np.loadtxt("Data/co2_mm_mlo.txt")
temp_data=np.loadtxt("Data/graph.txt", skiprows=5)
co2=co2_data[:,3]
beta=[]
x=0
avco2=[]
decimal_year=co2_data[:,2]
co24year=[]
temp=temp_data[78:,2]#slects only years from 1958 as we only have co2 data from then
#estimates average co2 for each year and appends it to avco2
year=1958
for i in decimal_year:
    integer=math.floor(decimal_year[x])
    if integer==year:
        co24year.append(co2[x])
    else:
        yearco2=np.array(co24year)
        mean=np.mean(yearco2)
        avco2.append(mean)
        co24year=[]
        year=year+1
    x=x+1
y=0
for i in range(1,63):
    fract=avco2[y+1]/avco2[y]
    beta.append(fract)
    y=y+1
print(beta)

print(len(beta))




#%% working out increase in radistive forcing and dq from year to year
#F=5.3ln(beta)
#tau=c/alpha 
#dT/dt=dq/mc
#dq=(F-Fg)*surface
#Fg=sigma*Tg^4
#radiative frocing only includes co2 and from 1958 as we are working on develloping the model while the others work on the data, so this is just the n=model without the data we will use
#we are hoping to add increase in rad forcing ude to increase in methane later
Tg=288
sigma=5.67e-8
Fg=sigma*Tg**4
print(Fg)
alpha=5.3
c=2

dF=alpha*np.log(beta)
print(dF)
surface=510e12

inc=np.sum(dF[0:39])

F0=1.46-inc #1.46 is rad forcing for CO" in 1998
print(F0)

dq=dF*surface
print(dq)
#dTt=
#plt.plot(time,temp)
#plt.plot(year[78:],dTt)
plt.show()
#%% calculating the increase in temeprature from year to year instead of since 1958
change=[]

for i in range(1,63):
    zz=temp[i]-temp[i-1]
    change.append(zz)
print(len(change))  

print(change)


#%% 
mc=(dq[0]-dq[61])/(change[0]-change[61])
print(mc)
#our estimate of mc is just taking last and first value for dq, since fitting the data did not work as trend is not strong enough
dTt=dq/mc
#print(len(dTt))
print(dTt)
Tt=[]
for i in range(0,len(dTt)): #working out the increase in temp since 1958 accordin to our model
    xyz=np.sum(dTt[0:i])
    Tt.append(xyz)
#print(len(Tt))
#print(Tt)
plt.plot(year1[78:],temp)
plt.plot(year1[79:],Tt)
plt.show()
print(Tt[61]) 
 #  trend is not that bad, could imporve on value of mc, though adding methane rad forcing might make up for the missing temp increase
#%%
#print(287+temp2[78])
#print(temp2[78])
#print(year1[78])
#%% 70% water, considering depths of 70m for well-mixed increases, 1000kgm^-3, heat cap 
c_water=4180 #Jkg^-1C-1
av_c_soil=1000 
depth_water=5
depth_soil=1.6
density_water=1000
density_soil=1300
mc_water=0.7*surface*depth_water*density_water*c_water
print(mc_water)
mc_soil=0.3*surface*depth_soil*density_soil*av_c_soil
print(mc_soil)
mc_tot=mc_water+mc_soil
print(mc_tot)
#%%

co2=co2_data[:,3]
time=year1[79:]
decimal_year=co2_data[:,2]
years=year1[78:]

#estimates average co2 for each year and appends it to avco2
#year=1958
def calc_beta(year):
    x=0
    co24year=[]
    beta=[]
    avco2=[]
    for i in decimal_year:
        integer=math.floor(decimal_year[x])
        if integer==year:
            co24year.append(co2[x])
        else:
            yearco2=np.array(co24year)
            mean=np.mean(yearco2)
            avco2.append(mean)
            co24year=[]
            year=year+1
        x=x+1
    y=0
    for i in range(1,63):
        fract=avco2[y+1]/avco2[y]
        beta.append(fract)
        y=y+1
    return beta

#print(beta)
#print(len(beta))
Beta=calc_beta(1958)

Tg=288
sigma=5.67e-8
Fg=sigma*Tg**4
surface=510e12
T_lw=0.2 #transmittance of the long wavelength of the atmosphere
#F_tot=[]
def forcing_CO2(alpha,beta):
    F=alpha*np.log(beta)
    return F

#olr of earth depends on sigma*T^4, so if T increases dT, olr increases by sigma*((T+dT)^4-dT^4)
## by taking base temperature as 287K, and thus 1958 temperature as 287K+anomaly


def yearly_temp_increase(number_years):
    T=287+temp2[78]
    inc_temp=[]
    excess_planetary_energy=[]
    dOLR=0
    for i in range (0,number_years):
        dF_CO2=forcing_CO2(5.3,Beta[i])
        dF_N2O=0
        dF_methane=0
        dF_tot=dF_CO2 +dF_N2O +dF_methane
        excess_planetary_energy=(dF_tot-dOLR)*surface
        dT=excess_planetary_energy/mc
        inc_temp.append(dT)
        dOLR=T_lw*sigma*((T+dT)**4-T**4)
        T=T+dT
    return inc_temp

def tot_temp_increase(timespan):
    array=yearly_temp_increase(timespan)
    tot_temp=[]
    for j in range (0,len(time)):
        xo=np.sum(array[0:j])
        tot_temp.append(xo)
    return tot_temp

#array_temp=yearly_temp_increase(len(time))
temperature=tot_temp_increase(len(time))
plt.plot(years,temp)
plt.plot(time,temperature)

plt.show()

print(temperature[len(time)-1])













  

