# -*- coding: utf-8 -*-
"""
Created on Wed May 26 12:30:11 2021

@author: nrace
"""

import numpy as np
import matplotlib.pyplot as plt

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
    plt.errorbar(time, conc, yerr = 0, fmt = 'r+', mew=2, ms=3, capsize = 2)
    plt.grid()
    plt.xlabel('Year')
    plt.ylabel('Concentration of gas_name, ppm' )
    plt.title('Concentration of gas_name Over Time', fontsize = 15)
    plt.plot(time, poly(time))
    
    plt.figure(figsize = (10, 7))
    plt.errorbar(x, poly(x[::12] , yerr = 0, fmt = 'c+', mew=2, ms=3, capsize = 2)
    plt.grid()
    plt.xlabel('Year')
    plt.ylabel('Concentration of gas_name, ppm' )
    plt.title('Extrapolated Concentration of gas_name Over Time', fontsize = 15)
    plt.plot(x, poly(x))




