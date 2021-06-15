# -*- coding: utf-8 -*-
"""
Created on Wed May 26 13:29:43 2021

@author: nyanraess
"""

# -*- coding: utf-8 -*-
"""
Created on Wed May 26 12:30:11 2021

@author: nrace
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op
import sympy as sym
import csv

#%%
#Section 1: Projecting Greenhouse gas levels

time, n, co2_conc, ch4_conc, n2o_conc = np.loadtxt('cleandata.csv', skiprows = 81, delimiter = ',', unpack = 1)

#importing co2 data from Mauna Loa; relevant columns have index 3 and 9
#column 9 is seasonally adjusted data


#%%
#Function's job is to assume gas is on a linear increase
#then plot a line of best fit over existing data
#then another plot extrapolating up to year 2070
#also returns an estimate of the yearly co2 increase based on the model

#to be able to use this function, your time and concentration values should be 1D arrays
#should work for data on any constituent of atmosphere
#if modelling for different gas, change the 'gas_name' input


def lin_projection(time: list, conc:list, gas_name_unit: str):
    fit, cov = np.polyfit(time, conc, 1, cov = 1)
    poly = np.poly1d(fit)
    x = np.linspace(2020, 2070, 1200)                      #600 divisions for 50 years i.e. approx monthly
    plt.figure(figsize = (10, 7))
    plt.errorbar(x[::24], poly(x[::24]) , yerr = 0, fmt = 'c+', mew=2, ms=3, capsize = 2) #marks levels at each year (every 12 months)
    plt.grid()
    plt.xlabel('Year')
    plt.ylabel('Concentration in ppm of: {}'.format(gas_name_unit))
    plt.title('Extrapolated Concentration over time of: {}'.format(gas_name_unit), fontsize = 15)
    plt.plot(x, poly(x))

    plt.figure(figsize = (10, 7))
    plt.errorbar(time, conc, yerr = 0, fmt = 'r+', mew=2, ms=3, capsize = 2)
    plt.grid()
    plt.xlabel('Year')
    plt.ylabel('Concentration in ppm of: {}'.format(gas_name_unit))
    plt.title('Concentration over time of: {}'.format(gas_name_unit), fontsize = 15)
    plt.plot(time, poly(time))
    print("Linear Projection Implies a CO2 increase in ppm, per year, of =%.3e +/- %.3e"%(fit[0], np.sqrt(cov[0,0])))
    lin_params = [fit[0], fit[1]]
    lin_proj_data = np.array([x, poly(x)])
    return lin_params, lin_proj_data
    

#returns parameters of linear fit as array
#returns projected data as array
#%%
#test the linear projection here; not great fit
lin_projection(time, co2_conc, gas_name_unit = 'CO2')

#%%

#function's job is to assume gas is on an exponential increase of form Ae^k(t^m), t=time
#produces all the same outputs as linear projection function
#requires good initial guesses to work
#initial guess array looks like int_guess = [A, k, m]
# note gas_name unit is of form 'CO2 (ppm)' or 'CH4 (ppb)' 

def exp_projection(time: list, conc:list, gas_name_unit: str, int_guess: list):
    x = np.linspace(2020, 2070, 1200)
    def exp_func(time, *int_guess):
        return int_guess[0]*np.e**(int_guess[1]*((time-1960)**int_guess[2]))                                           #plots twice a month
    params, pcov = op.curve_fit(exp_func, time, conc, int_guess)  
    plt.figure(figsize = (10, 7))
    plt.errorbar(x[::24], exp_func((x[::24]), *params), yerr = 0, fmt = 'y+', mew=2, ms=3, capsize = 2)
    plt.grid()
    plt.xlabel('Year')
    plt.ylabel('Concentration of: {}'.format(gas_name_unit))
    plt.title('Extrapolated Concentration over time of: {}'.format(gas_name_unit), fontsize = 15)
    plt.plot(x, exp_func(x, *params), 'r', linewidth = 2)

    plt.figure(figsize = (10, 7))
    plt.errorbar(time, conc, yerr = 0, fmt = 'r+', mew=2, ms=2, capsize = 2)
    plt.grid()
    plt.xlabel('Year')
    plt.ylabel('Concentration of: {}'.format(gas_name_unit))
    plt.title('Concentration over time of: {}'.format(gas_name_unit), fontsize = 15)
    plt.plot(time, exp_func(time, *params), 'b', linewidth = 3)
    print("Exponential Projection Implies an increase, per year, by a factor of approx: {}".format(np.e**params[1])) 
    print(params, ' = parameters' , pcov, ' = covariance')
    
    proj_data = [x, exp_func(x, *params)]
    data = np.array([time, exp_func(time, *params)])
    return proj_data, data, params
                        #returns projected data as an array
                        #returns fit to real data as 2d array
                        #also returns exponential parameters as array
    
#%%

#%%
#try the exponential projection here
#correspondence is really good

co2_exp_proj, co2_fit_data, co2_params = exp_projection(time, co2_conc, 'CO2 (ppm)', int_guess = [2.465e-3, 5.99e-3, 1])


#%%
#This cell for n20 (nitrous oxide) projections
    #convert to ppm for consistency

#n2o_lin_proj = lin_projection(time, n2o_conc, 'N2O')

n2o_exp_proj, n2o_fit_data, n2o_params = exp_projection(time, n2o_conc, 'N2O (ppb)', [0, 0.005, 1])

n2o_proj_data = n2o_exp_proj


#exponential fit still seems to win

#%%
#This cell for ch4 (methane) projections
#ch4 data follows no obvious pattern; might be worth looking into why this is
#neither fit particularly good
#exponential fit reduces to approximately linear in future; probably better


#ch4_lin_proj = lin_projection(time, ch4_conc, 'CH4')

ch4_exp_proj, ch4_fit_data, ch4_params = exp_projection(time, ch4_conc, 'CH4 (ppb)', [0, 0.005, 1])

print(co2_exp_proj)
print(ch4_exp_proj)
print(n2o_exp_proj)

#ch4_proj_data = n2o_exp_proj

#%%
#This cell for sulphate projections
#data only goes to 1996 so fit data from 1930 to be consistent with others
#i.e. projecting based on approx 60 years of data
#volcanoes make data too spiky

sulph_data = np.loadtxt('sulphates.txt', delimiter = ',', skiprows = 0, unpack = 0)

sulph_time = sulph_data[0:,0]
sulph_conc = sulph_data[0:,1]    

sulph_lin_proj = lin_projection(sulph_time, sulph_conc, 'Sulphates (ppb)')
sulph_exp_proj = exp_projection(sulph_time, sulph_conc, 'Sulphates (ppb)', [0.3, 0.0005, 1])


#%%

time, n, co2_conc, ch4_conc, n2o_conc = np.loadtxt('cleandata.csv', skiprows = 81, delimiter = ',', unpack = 1)




def exp_projection(time: list, conc:list, gas_name_unit: str, int_guess: list):
    x = np.linspace(2020, 2070, 1200)
    def exp_func(time, *int_guess):
        return int_guess[0]*np.e**(int_guess[1]*((time-1960)**int_guess[2]))                                           #plots twice a month
    params, pcov = op.curve_fit(exp_func, time, conc, int_guess)  
    plt.figure(figsize = (10, 7))
    plt.errorbar(x[::24], exp_func((x[::24]), *params), yerr = 0, fmt = 'y+', mew=2, ms=3, capsize = 2)
    plt.grid()
    plt.xlabel('Year')
    plt.ylabel('Concentration of: {}'.format(gas_name_unit))
    plt.title('Extrapolated Concentration over time of: {}'.format(gas_name_unit), fontsize = 15)
    plt.plot(x, exp_func(x, *params), 'r', linewidth = 2)

    plt.figure(figsize = (10, 7))
    plt.errorbar(time, conc, yerr = 0, fmt = 'r+', mew=2, ms=2, capsize = 3)
    plt.grid()
    plt.xlabel('Year')
    plt.ylabel('Concentration of: {}'.format(gas_name_unit))
    plt.title('Concentration over time of: {}'.format(gas_name_unit), fontsize = 15)
    plt.plot(time, exp_func(time, *params), 'b', linewidth = 3)
    print("Exponential Projection Implies an increase, per year, by a factor of approx: {}".format(np.e**params[1])) 
    print(params, ' = parameters' , pcov, ' = covariance')
    
    proj_data = [x, exp_func(x, *params)]
    data = np.array([time, exp_func(time, *params)])
    return proj_data, data, params

#%%
co2_exp_proj, co2_fit_data, co2_params = exp_projection(time, co2_conc, 'CO2 (ppm)', int_guess = [2.465e-3, 5.99e-3, 1])

n2o_exp_proj, n2o_fit_data, n2o_params = exp_projection(time, n2o_conc, 'N2O (ppb)', [0, 0.005, 1])

ch4_exp_proj, ch4_fit_data, ch4_params = exp_projection(time, ch4_conc, 'CH4 (ppb)', [0, 0.005, 1])
#%%
print(co2_exp_proj)
print(ch4_exp_proj)
print(n2o_exp_proj)

#%%
#Both CFCs have been on a linear decrease since the 90s; projected from then

cfc11_data = np.loadtxt('cfc11_data.csv', skiprows = 115, delimiter = ',', unpack = 0)
cfc11_conc = cfc11_data[0:len(cfc11_data)-1, 1]
cfc11_linparams, cfc11_proj_data = lin_projection(np.linspace(1994, 2020, 27), cfc11_conc, 'CFC-11 (ppb)')

cfc12_data = np.loadtxt('cfc12_data.csv', skiprows = 126, delimiter = ',', unpack = 0)
cfc12_conc = cfc12_data[0:len(cfc12_data)-1, 1]
cfc12_proj_data, cfc12_fit_data, cfc12_params = exp_projection(np.linspace(2002, 2020, 16), cfc12_conc, 'CFC-12 (ppb)',[0.6, -0.0001, 1])

#%%
from numpy import diff
from scipy.optimize import fsolve

#Target function takes in real data, projected data, target concnetration, target year
#i.e. by when the target concentration should be achieved, the parameters associated 
#with the exponential projection; works fine with ch4, co2, n2o


def exp(x, a, k, m):
    return a*np.e**(k*(x-1960)**(m)) 

def target_func(real_data:list, fit_data: list, target_conc: float, target_year: int,
                gas_params: list, gas_name_unit: str):
    T = target_conc
    Y = target_year
    fit_t = fit_data[0]
    fit_conc = fit_data[1]
    C_f = fit_conc[-1]
    dconc_dt = diff(fit_conc)/diff(fit_t)
    t = np.linspace(1990, target_year, 12*(target_year-1990))
    G = dconc_dt[-1]
    def func(x):
        A = x[0]
        B = x[1]
        C = x[2]
        D = x[3]
        return [3*A*Y**2 + 2*B*Y + C,
                3*A*1990**2 + 2*B*1990+C-G,
                ((A*Y**3)+(B*Y**2)+(C*Y+D))-T,
                (A*1990**3+B*1990**2+C*1990+D)-C_f]
    
    a, b, c, d = fsolve(func, [1e-6, 1e-6, 1e-6, 1e-6])
    y_t = a*t**3+b*t**2+c*t+d
    plt.figure(figsize = (10,7))
    plt.grid()
    plt.plot(fit_t, exp(fit_t, *gas_params), 'b', linewidth = 3)
    plt.plot(t, a*t**3+b*t**2+c*t+d, 'r', linewidth = 3)
    plt.errorbar(np.linspace(1960, 1990, 31), real_data, yerr = 0, fmt = 'm+', mew=2.5, ms=2.5, capsize = 2.5)
    plt.xlabel('Year')
    plt.ylabel('Concentration of: {}'.format(gas_name_unit))
    plt.suptitle('Concentration over time of: {}'.format(gas_name_unit), fontsize = 18, weight = 'bold', y=1)
    plt.legend(['Past Data fit', 'Future projection given constraints and current trends'], fontsize = 9)
    plt.title('Constraint: an atmospheric concentration of {} of {} by the year {}'
              .format(gas_name_unit, str(target_conc), str(target_year), fontsize = 14, y=1.07))
    return np.array([t, y_t])
    
    
    
    
#%%

#%%
   
co2_targ = target_func(co2_conc[:31:], co2_fit_data[:31:], 450, 2070, co2_params, 'CO2 (ppm)')

n2o_targ = target_func(n2o_conc[:31:], n2o_fit_data[:31:], 360, 2070, n2o_params, 'N2O (ppb)')

ch4_targ = target_func(ch4_conc[:31:], ch4_fit_data[:31:], 2000, 2070, ch4_params, 'CH4 (ppb)')


#%%

time_proj = co2_exp_proj[0]

import xlsxwriter

workbook = xlsxwriter.Workbook('Target Projections.xlsx')
worksheet = workbook.add_worksheet()

array = [time_proj[::24],
         co2_targ[1,::12],
         n2o_targ[1,::12],
         ch4_targ[1,::12] ]

row = 0

for col, data in enumerate(array):
    worksheet.write_column(row, col, data)

workbook.close()




#%%
#Saving all projections so far into excel file

co2_proj = co2_exp_proj[1]
ch4_proj = ch4_exp_proj[1]
n2o_proj = n2o_exp_proj[1]
cfc11_proj = cfc11_proj_data[1]
cfc12_proj = cfc12_proj_data[1]


workbook = xlsxwriter.Workbook('Projections.xlsx')
worksheet = workbook.add_worksheet()

array = [time_proj,
         co2_proj,
         ch4_proj,
         n2o_proj,
         cfc11_proj,
         cfc12_proj]

row = 0

for col, data in enumerate(array):
    worksheet.write_column(row, col, data)

workbook.close()

#%%
time, WMGHG, Ozone,	Solar, Land_Use, SnowAlb_BC, Orbital, TropAerDir, TropAerInd, StratAer = np.loadtxt('instant_forcings_1880_to_2020.csv', skiprows = 131, max_rows = 53, delimiter = ',', unpack = 1)

land_use_params, land_use_proj_data, = lin_projection(time, Land_Use, 'Land Use Radiative Forcing (W/m\u00b2)')

ozone_params, ozone_proj_data = lin_projection(time, Ozone, 'Ozone Radiative Forcing (W/m\u00b2)')

snow_params, snow_proj_data = lin_projection(time, SnowAlb_BC, 'Snow Albedo Radiative Forcing (W/m\u00b2)')

tropaerdir_params, tropaerdir_data = lin_projection(time, TropAerDir, 'Tropospheric Aerosols Direct Radiative Forcing (W/m\u00b2)')

tropaerind_params, tropaerind_proj_data = lin_projection(time, TropAerInd, 'Tropospheric Aerosols Indirect Radiative Forcing (W/m\u00b2)')
#%%
workbook = xlsxwriter.Workbook('Forcing Projections.xlsx')
worksheet = workbook.add_worksheet()

array = [np.linspace(2020, 2070, 51),
         land_use_proj_data[1],
         ozone_proj_data[1],
         snow_proj_data[1],
         tropaerdir_data[1],
         tropaerind_proj_data[1]]

row = 0

for col, data in enumerate(array):
    worksheet.write_column(row, col, data)

workbook.close()

#%%
