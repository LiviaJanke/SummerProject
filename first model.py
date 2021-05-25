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
year,temp1,temp2=sp.loadtxt('Data/graph.txt',skiprows=5,unpack=True)

plt.plot(year[78:],temp2[78:])

#%% setting 1958 as temperature 0 and working out increase in temp from there
temp=temp2[78:]-temp2[78]
time=year[78:]
plt.plot(time,temp)
print(np.max(temp))
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

Tt=[]
for i in range(0,len(dTt)): #working out the increase in temp since 1958 accordin to our model
    xyz=np.sum(dTt[0:i])
    Tt.append(xyz)
#print(len(Tt))
#print(Tt)
plt.plot(year[78:],temp)
plt.plot(year[79:],Tt)
plt.show()
    
 #  trend is not that bad, could imporve on value of mc, though adding methane rad forcing might make up for the missing temp increase
    
    

