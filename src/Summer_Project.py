# -*- coding: utf-8 -*-
"""
Created on Wed May 26 13:29:43 2021

@author: nrace
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

co2_time = Co2_data[0:,3]      #extract columns of 3rd and 9th index
co2_conc = Co2_data[0:,9]


def lin_projection(time: list, conc:list, gas_name: str):
    fit, cov = np.polyfit(time, conc, 1, cov = 1)
    poly = np.poly1d(fit)
    x = np.linspace(2020, 2070, 600)
    plt.figure(figsize = (10, 7))
    plt.errorbar(x[::12], poly(x[::12]) , yerr = 0, fmt = 'c+', mew=2, ms=3, capsize = 2)
    plt.grid()
    plt.xlabel('Year')
    plt.ylabel('Concentration in ppm of: {}'.format(gas_name))
    plt.title('Concentration over time of: {}'.format(gas_name), fontsize = 15)
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
lin_projection(co2_time, co2_conc, gas_name = 'CO2')

#%%
def exp_projection(time: list, conc:list, gas_name: str, int_guess: list):
    def exp_func(time, *int_guess):
        return int_guess[0]*np.e**(int_guess[1]*(time-1960))
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
    print("Exponential Projection Implies a CO2 increase in ppm, per year, of =%.3e +/- %.3e"%(np.e**params[1], np.sqrt(pcov[1,1])))
    print(params, ' = parameters' , pcov, ' = covariance')
    
#%%
   #Functions job is to automate the initial guesses 
    
def exp_int_guess_solver(time:list, conc:list):
    A,k = sym.symbols('A,k')
    eq1 = sym.Eq(A*np.e**(time[0]*k), conc[0])
    eq2 = sym.Eq(A*np.e**(time[360]*k), conc[360])
    result = sym.solve([eq1,eq2],(A,k))
    print(result)

exp_int_guess_solver(co2_time, co2_conc)
#%%
exp_projection(co2_time, co2_conc, 'CO2', int_guess = [2.465e-3, 5.99e-3, 1])

#%%




















