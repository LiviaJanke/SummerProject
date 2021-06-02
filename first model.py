# -*- coding: utf-8 -*-
"""
Created on Wed May 19 16:07:22 2021

@author: louis
"""
# https://github.com/LiviaJanke/SummerProject.git
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
temper=temp2[78:]-temp2[78]
time=year1[79:]
plt.plot(year1[78:],temper)
print(np.max(temper))
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
print(avco2[1])
for i in range(1,63):
    fract=avco2[y+1]/avco2[y]
    beta.append(fract)
    y=y+1
#print(beta)

#print(len(beta))




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
#%% increase in temp due to co2 beteen 1958 and 2020
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
    for i in range(1,len(time)+1):
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
        dT=excess_planetary_energy/mc_tot
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

plt.plot(years,temper)
plt.plot(time,temperature)

plt.show()

print(temperature[len(time)-1])

#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op
import sympy as sym
import csv

#%%  predicting with only co2 if trend continues as it does now, temp increase from 2020 to 2070
#Section 1: Projecting Greenhouse gas levels

Co2_data = np.loadtxt('monthly_flask_co2_mlo.csv', skiprows = 59, delimiter = ',', unpack = 0)

#importing co2 data from Mauna Loa; relevant columns have index 3 and 9
#column 9 is seasonally adjusted data

co2_time = Co2_data[0:,3]      #extract columns of 3rd and 9th index
co2_conc = Co2_data[0:,9]


#function's job is to assume gas is on an exponential increase of form Ae^k(t^m), t=time
#produces all the same outputs as linear projection function
#requires good initial guesses to work
#initial guess array looks like int_guess = [A, k, m]
# {m=1 turns out to be ideal, but for some reason the fit fails if you remove int_guess[2] from exp_func
#and try to proceed without the m, even though it's pointless}


def exp_projection(time: list, conc:list, gas_name: str, int_guess: list):
    def exp_func(time, *int_guess):
        return int_guess[0]*np.e**(int_guess[1]*((time-1960)**int_guess[2]))
    x = np.linspace(2020, 2070, 51)
    params, pcov = op.curve_fit(exp_func, time, conc, int_guess)  
    plt.figure(figsize = (10, 7))
    plt.errorbar(x[::12], exp_func((x[::12]), *params), yerr = 0, fmt = 'y+', mew=2, ms=3, capsize = 2)
    plt.grid()
    plt.xlabel('Year')
    plt.ylabel('Concentration in ppm of: {}'.format(gas_name))
    plt.title('Concentration over time of: {}'.format(gas_name), fontsize = 15)
    plt.plot(x, exp_func(x, *params), 'r', linewidth = 2)

    plt.figure(figsize = (10, 7))
    plt.errorbar(time, conc, yerr = 0, fmt = 'r+', mew=1, ms=2, capsize = 2)
    plt.grid()
    plt.xlabel('Year')
    plt.ylabel('Concentration in ppm of: {}'.format(gas_name))
    plt.title('Concentration over time of: {}'.format(gas_name), fontsize = 15)
    plt.plot(time, exp_func(time, *params), 'b', linewidth = 3)
    #print("Exponential Projection Implies a CO2 increase in ppm, per year, by a factor of approx: {}".format(np.e**params[1])) 
    #print(params, ' = parameters' , pcov, ' = covariance')
    projected_data=exp_func(x,*params)
    return x,projected_data

#print(exp_projection(co2_time, co2_conc, 'CO2', int_guess = [2.465e-3, 5.99e-3, 1]))
w,z=exp_projection(co2_time, co2_conc, 'CO2', int_guess = [2.465e-3, 5.99e-3, 1])

co2=z
time=w
decimal_year=w

def calc_beta(year):
    x=0
    co24year=[]
    beta=[]
    avco2=z
    #for i in decimal_year:
        #integer=math.floor(decimal_year[x])
        #if integer==year:
            #co24year.append(co2[x])
        #else:
            #yearco2=np.array(co24year)
            #mean=np.mean(yearco2)
            #avco2.append(mean)
            #co24year=[]
            #year=year+1
        #x=x+1
    y=0
    for i in range(1,len(time)):
        fract=avco2[y+1]/avco2[y]
        beta.append(fract)
        y=y+1
    return beta

#print(beta)
#print(len(beta))
Beta=calc_beta(2020)

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
    T=288
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
temperature=tot_temp_increase(len(time)-1)
#plt.plot(years,temp)
plt.plot(time,temperature)

plt.show()

print(temperature[len(time)-1])
print(len(temperature))
#print(w[len(time)-1])
#%%

#%%  temp increase becasue of co2 between 1880 and 2020


t_gas,tt,co2,ch4,n2o=np.loadtxt('Data/cleandata.csv',skiprows=1,delimiter=',',unpack=True)
years,temp_no,temp=np.loadtxt('Data/graph.txt',skiprows=5,unpack=True)
decimal_year=t_gas
time=years[1:]
#estimates average co2 for each year and appends it to avco2
#year=1958
def calc_beta(year):
    x=0
    co24year=[]
    beta=[]
    avco2=co2
    y=0
    print(len(avco2))
    print(avco2[0])
    for i in range(1,len(time)+1):
        fract=avco2[y+1]/avco2[y]
        beta.append(fract)
        y=y+1
    return beta

#print(beta)
#print(len(beta))
Beta=calc_beta(1880)
print(len(Beta))
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
    T=286.1
    inc_temp=[]
    excess_planetary_energy=[]
    dOLR=0
    for i in range (0,number_years):
        dF_CO2=forcing_CO2(5.3,Beta[i])
        dF_N2O=0
        dF_methane=0
        dF_tot=dF_CO2 +dF_N2O +dF_methane
        excess_planetary_energy=(dF_tot-dOLR)*surface
        dT=excess_planetary_energy/mc_tot
        inc_temp.append(dT)
        dOLR=T_lw*sigma*((T+dT)**4-T**4)/10
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

#%%
plt.plot(t_gas,n2o)
plt.title('n2o')
plt.show()
plt.plot(t_gas,ch4)
plt.title('ch4')
plt.show()

#%% model with eq temp and correct olr but start olr 0
t_gas,tt,co2,ch4,n2o=np.loadtxt('Data/cleandata.csv',skiprows=1,delimiter=',',unpack=True)
years,temp_no,temp=np.loadtxt('Data/graph.txt',skiprows=5,unpack=True)
decimal_year=t_gas
time=years[1:]

surface=510e12
T_lw=0.2 #transmittance of the long wavelength of the atmosphere
B=2.218 #Wm^-2K^-1 constant for OLR, determined empirically
alpha_co2=5.35
ach4=0.036
an2o=0.12

#Beta=calc_beta(1880)
def forcing_CO2(data):
    avco2=data
    dfco2=[]
    for i in range(0,len(time)):
        beta=avco2[i+1]/avco2[i]
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
        dfch4 -= (meth_data[i])**0.5
        dfch4 *= ach4
        dfch4 -= f(meth_data[i+1],n2o_data[i+1])-f(meth_data[i],n2o_data[i])
        dfch4_list.append(dfch4)         
    return dfch4_list
      
def forcing_n2o(meth_data,n2o_data):
    dfn2o_list=[]
    for i in range(len(meth_data)-1):
        dfn2o=0
        dfn2o+=(n2o_data[i+1])**0.5-(n2o_data[i]**0.5)
        dfn2o*=0.12
        dfn2o-=f(meth_data[i],n2o_data[i+1])-f(meth_data[i],n2o_data[i])
        
        dfn2o_list.append(dfn2o)
    return dfn2o_list
print(forcing_methane(ch4,n2o))

## by taking base temperature as 287K, and thus 1958 temperature as 287K+anomaly
# according to penas report, the relationship between olr and temp is as follows: A+B*T where A=-339.647 and B=2.218 
mc_tot=calc_mc_tot(70, 1.6)

def temp_increase(number_years):
    T_eq=287
    T=286.1
    temperature=[]
    inc_temp=[]
    excess_planetary_energy=[]
    dOLR=0 #B*(T_eq-0.09)-B*T_eq 
    dF_CO2=forcing_CO2(co2)
    dF_N2O=forcing_n2o(ch4,n2o)
    dF_methane=forcing_methane(ch4,n2o)
    for i in range (0,number_years):
        dF_tot=+dF_CO2[i]+dF_N2O[i]+dF_methane[i]
        excess_planetary_energy=(dF_tot-dOLR)*surface
        dT=excess_planetary_energy/mc_tot
        inc_temp.append(dT)
        dOLR=B*dT
        T=T+dT
        temperature.append(T-286.1)
    return temperature





#array_temp=yearly_temp_increase(len(time))
temperature=temp_increase(len(time))
plt.plot(years,temp)
plt.plot(time,temperature)

plt.show()

print(temperature[len(time)-1])

#%% same as above but with start olr -0.09*B
t_gas,tt,co2,ch4,n2o=np.loadtxt('Data/cleandata.csv',skiprows=1,delimiter=',',unpack=True)
years,temp_no,temp=np.loadtxt('Data/graph.txt',skiprows=5,unpack=True)
decimal_year=t_gas
time=years[1:]

surface=510e12
T_lw=0.2 #transmittance of the long wavelength of the atmosphere
B=2.218 #Wm^-2K^-1 constant for OLR, determined empirically
alpha_co2=5.35
ach4=0.036
an2o=0.12

#Beta=calc_beta(1880)
def forcing_CO2(data):
    avco2=data
    dfco2=[]
    for i in range(0,len(time)):
        beta=avco2[i+1]/avco2[i]
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
        dfch4 -= (meth_data[i])**0.5
        dfch4 *= ach4
        dfch4 -= f(meth_data[i+1],n2o_data[i+1])-f(meth_data[i],n2o_data[i])
        dfch4_list.append(dfch4)         
    return dfch4_list
      
def forcing_n2o(meth_data,n2o_data):
    dfn2o_list=[]
    for i in range(len(meth_data)-1):
        dfn2o=0
        dfn2o+=(n2o_data[i+1])**0.5-(n2o_data[i]**0.5)
        dfn2o*=0.12
        dfn2o-=f(meth_data[i],n2o_data[i+1])-f(meth_data[i],n2o_data[i])
        
        dfn2o_list.append(dfn2o)
    return dfn2o_list


## by taking base temperature as 287K, and thus 1958 temperature as 287K+anomaly
# according to penas report, the relationship between olr and temp is as follows: A+B*T where A=-339.647 and B=2.218 
mc_tot=calc_mc_tot(70, 1.6)

def temp_increase(number_years):
    T_eq=287
    T=286.1
    temperature=[]
    inc_temp=[]
    excess_planetary_energy=[]
    dOLR=B*(-0.09) 
    dF_CO2=forcing_CO2(co2)
    dF_N2O=forcing_n2o(ch4,n2o)
    dF_methane=forcing_methane(ch4,n2o)
    
    for i in range (0,number_years):
        dF_tot=+dF_CO2[i]+dF_N2O[i]+dF_methane[i]
        excess_planetary_energy=(dF_tot-dOLR)*surface
        dT=excess_planetary_energy/mc_tot
        inc_temp.append(dT)
        dOLR=B*dT
        T=T+dT
        temperature.append(T-286.1)
    return temperature





#array_temp=yearly_temp_increase(len(time))
temperature=temp_increase(len(time))
plt.plot(years,temp)
plt.plot(time,temperature)

plt.show()

print(temperature[len(time)-1])


#%% model with C_0 as the first data point rather than the previous year
t_gas,tt,co2,ch4,n2o=np.loadtxt('Data/cleandata.csv',skiprows=1,delimiter=',',unpack=True)
years,temp_no,temp=np.loadtxt('Data/graph.txt',skiprows=5,unpack=True)
decimal_year=t_gas
time=years[1:]

surface=510e12
T_lw=0.2 #transmittance of the long wavelength of the atmosphere
B=2.218 #Wm^-2K^-1 constant for OLR, determined empirically
alpha_co2=5.35
ach4=0.036
an2o=0.12


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
mc_tot=calc_mc_tot(70, 1.6)

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


