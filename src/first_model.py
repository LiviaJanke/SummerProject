#this short piece of code basically just generates the temperature of the earth
# as the transmittance of short and long wave radiation of the atmosphere
#changes. This just has each transmittance factors
#increasing and decreasing with time as co2 increases
import numpy as np
import matplotlib.pyplot as plt
#sw and lw represent the transmittance of short and long wave radiation for CO2
#respectively every year for the next 50 years
#co2 increase assumed to double over next 50 years
co2=sw=np.linspace(4.12e-4, 8.24e-4,num=50)
#using numbers from the atmo book and assuming co2 is a perfect greenhouse gas
#i estimated the transmitance values for co2 and rest of the atmosphere without co2 and
#came up with this formula for the transmittance of the atmosphere as a whole
#as a function of co2 for shortwave(sw) and longwave(lw) radiation
sw=co2+(1-co2)*0.89935
lw=(1-co2)*0.2008
time=np.linspace(1,50,num=50)
#takes value for transmission coefficents and calculates equilibrium temperature, 240 represents irradiance of the sun and other number is the boltzmann constant
def temp(s,l):
    value=((1+s)*240/(5.67e-8*(1+l)))**0.25-273
    return value
print(temp(sw,lw))
plt.plot(time,temp(sw,lw))
#is this working