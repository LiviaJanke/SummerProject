import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
import math
#%%

clean_data=np.loadtxt("Data/cleandata.csv", delimiter=',', skiprows=1)
ch4=clean_data[1:,3]
ach4=0.036
an2o=0.12
n2o=clean_data[1:,4]
def f(M,N):
    f_function=0.47*np.log(1+2.01e-5*(M*N)**0.75+5.31e-15*M*(M*N)**1.52)
    return f_function
def forcing_methane(meth_data,n2o_data):
    dfch4_list=[]
    for i in range(len(meth_data)-1):
        # dfch4=ach4(
        #     (meth_data[i+1])**0.5
        #     -(meth_data[i])**0.5)
        # -(f(meth_data[i+1],n2o_data[i+1]-f(meth_data[i],n2o_data[i]))
        dfch4 = 0
        dfch4 += (meth_data[i+1])**0.5
        dfch4 -= (meth_data[i])**0.5
        dfch4 -= f(meth_data[i+1],n2o_data[i+1]-f(meth_data[i],n2o_data[i]))
        dfch4 *= ach4
        dfch4_list.append(dfch4)                                                       
    return dfch4_list
def forcing_n2o(meth_data,n2o_data):
    dfn2o_list=[]
    for i in range(len(meth_data)-1):
        dfn2o=0
        dfn2o+=(n2o_data[i+1])**0.5-(n2o_data[i]**0.5)
        dfn2o-=f(meth_data[i],n2o_data[i+1])-f(meth_data[i],n2o_data[i])
        dfn2o*=0.12
        dfn2o_list.append(dfn2o)
    return dfn2o_list
print(forcing_n2o(ch4,n2o))
print(len(ch4))
#