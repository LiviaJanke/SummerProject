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

Co2_data = np.loadtxt('monthly_flask_co2_mlo.csv', skiprows = 59, delimiter = ',', unpack = 0)

#importing co2 data from Mauna Loa; relevant columns have index 3 and 9
#column 9 is seasonally adjusted data

co2_time = Co2_data[0:,3]      #extract columns of 3rd and 9th index
co2_conc = Co2_data[0:,9]

#%%
#Function's job is to assume gas is on a linear increase
#then plot a line of best fit over existing data
#then another plot extrapolating up to year 2070
#also returns an estimate of the yearly co2 increase based on the model

#to be able to use this function, your time and concentration values should be 1D arrays
#should work for data on any constituent of atmosphere
#if modelling for different gas, change the 'gas_name' input


def lin_projection(time: list, conc:list, gas_name: str):
    fit, cov = np.polyfit(time, conc, 1, cov = 1)
    poly = np.poly1d(fit)
    x = np.linspace(2020, 2070, 600)                      #600 divisions for 50 years i.e. approx monthly
    plt.figure(figsize = (10, 7))
    plt.errorbar(x[::12], poly(x[::12]) , yerr = 0, fmt = 'c+', mew=2, ms=3, capsize = 2) #marks levels at each year (every 12 months)
    plt.grid()
    plt.xlabel('Year')
    plt.ylabel('Concentration in ppm of: {}'.format(gas_name))
    plt.title('Extrapolated Concentration over time of: {}'.format(gas_name), fontsize = 15)
    plt.plot(x, poly(x))

    plt.figure(figsize = (10, 7))
    plt.errorbar(time, conc, yerr = 0, fmt = 'r+', mew=2, ms=3, capsize = 2)
    plt.grid()
    plt.xlabel('Year')
    plt.ylabel('Concentration in ppm of: {}'.format(gas_name))
    plt.title('Concentration over time of: {}'.format(gas_name), fontsize = 15)
    plt.plot(time, poly(time))
    print("Linear Projection Implies a CO2 increase in ppm, per year, of =%.3e +/- %.3e"%(fit[0], np.sqrt(cov[0,0]))) 
#%%
#test the linear projection here; not great fit
lin_projection(co2_time, co2_conc, gas_name = 'CO2')

#%%

#function's job is to assume gas is on an exponential increase of form Ae^k(t^m), t=time
#produces all the same outputs as linear projection function
#requires good initial guesses to work
#initial guess array looks like int_guess = [A, k, m]
# {m=1 turns out to be ideal, but for some reason the fit fails if you remove int_guess[2] from exp_func
#and try to proceed without the m, even though it's pointless}


def exp_projection(time: list, conc:list, gas_name: str, int_guess: list):
    def exp_func(time, *int_guess):
        return int_guess[0]*np.e**(int_guess[1]*((time-1960)**int_guess[2]))
    x = np.linspace(2020, 2070, 600)
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
    print("Exponential Projection Implies a CO2 increase in ppm, per year, by a factor of approx: {}".format(np.e**params[1])) 
    print(params, ' = parameters' , pcov, ' = covariance')
    
    data = exp_func(x, *params)
    return data
                        #returns projected data as an array
    
#%%
   #Functions job is to automate the initial guesses 
   #produces set of 2 simultaenous equations with m=1 assumed; m/=1 produces usually unsolvable set 
   #of 3 equations
   #found my actual initial guesses (see cell below) by hand
   
#READ THIS BIT  --->   #because this function usually causes python to crash for me
                       #possibly because of how python solves equations of this form?
    
def exp_int_guess_solver(time:list, conc:list):
    A,k = sym.symbols('A,k')
    eq1 = sym.Eq(A*np.e**(time[0]*k), conc[0])
    eq2 = sym.Eq(A*np.e**(time[360]*k), conc[360])
    result = sym.solve([eq1,eq2],(A,k))
    print(result)

exp_int_guess_solver(co2_time, co2_conc)
#%%
#try the exponential projection here
#correspondence is really good

projected_data = exp_projection(co2_time, co2_conc, 'CO2', int_guess = [2.465e-3, 5.99e-3, 1])
#print(projected_data)

#%%




















