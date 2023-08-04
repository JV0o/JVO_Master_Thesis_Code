# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 15:57:53 2023

@author: s210212
"""
#%% Import models
#from mmodels.monod_bb.im_monod_bb import im_monod_bb
#from mmodels.xu_bb.im_xu_bb import im_xu_bb
#from mmodels.anane_bb.im_anane_bb import im_anane_bb # If this doesn't work run the two lines below for anane model.
import sys
sys.path.append('C:/Users/s210212/Documents/DTU/Thesis/Python_code/Xu_IF_allvol/mmodels')
from anane_bb.im_anane_bb import im_anane_bb
from xu_bb.im_xu_bb import im_xu_bb

#from mmodels.anane_bb.im_anane_bb import im_anane_bb
#%% 
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import os
import pandas as pd


#%%
def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)
################################## XU MODEL ##################################

############# Batch mode #######################
#%%
method=0
im_mmodel = im_xu_bb(
    strain_id='DDB35',
    strain_description= 'WT E. coli'
    )

## WT strain
'''
im_mmodel.define_strain_params(
    qS_max=np.mean([0.74036, 0.72265]),#np.mean([0.81528,0.7752]),## #1.3 np.mean([0.74036, 0.72265])
    qm_max=0, #0.04
    qA_c_max=np.mean([0.03444,0.036252]), #0.15
    qO_max=np.mean([6.7453,6.9238]), #15 np.mean([5.9305,5.8531]) np.mea(n[5.4322, 5.5012])
    Ysx_ox=np.mean([0.4592,0.4714]), # np.mean([0.4592,0.4714])    np.mean([0.43349, 44498])
    Ysx_of=np.mean([0.4592,0.4714]), # np.mean([0.4592,0.4714])      0.15 np.mean([0.46473,0.50662])-np.mean([0.41704,0.43812]) np.mean([0.41704,0.43812])
    Ysa=0.667, #0.667 np.mean([0.25974,0.23213])
    Yax=0, #0.4
    Ki_s=100, #5
    Ks=0.05, #0.05
    Ka=0.05, #0.05
    Ki_o=100 #4
)
'''

#New parameters with HPLC better calibration
im_mmodel.define_strain_params(
    qS_max=np.mean([0.78433, 0.76553]),#np.mean([0.81528,0.7752]),## #1.3 np.mean([0.74036, 0.72265])
    qm_max=0, #0.04
    #qA_c_max=np.mean([0.03444,0.036252]), #0.15
    qA_c_max=np.mean([0.03444,0.03726]), #0.15
    qO_max=np.mean([6.7453,6.9238]), #15 np.mean([5.9305,5.8531]) np.mea(n[5.4322, 5.5012])
    Ysx_ox=np.mean([0.43349, 0.44498]), # np.mean([0.4592,0.4714])    np.mean([0.43349, 44498])
    Ysx_of=np.mean([0.43349, 0.44498]), # np.mean([0.4592,0.4714])      0.15 np.mean([0.46473,0.50662])-np.mean([0.41704,0.43812]) np.mean([0.41704,0.43812])
    Ysa=0.667, #0.667 np.mean([0.25974,0.23213])
    Yax=0, #0.4
    Ki_s=100, #5
    Ks=0.05, #0.05
    Ka=0.05, #0.05
    Ki_o=100 #4
)


# ODE functions
func = im_mmodel.im_bmodel_odes
# initial values
y0 = [
    0.3*0.072,    # X [g/L] (OD WT_ini=0.072 HMP_ini=0.036 CDW OD factor=0.3)
    20,     # S [g/L]
    0.2095, # O [mol/L]
    0      # A [g/L]
]
# timesteps to evaluate
t_eval = np.linspace(0,35,100001)

# process parameter function arguments
kla = 200 # kla [1/h] # WT:200 HMP:150
#pabs = 1 # pabs [bar]
#yO = 0.2095 # yO [-]
#Sf = 0 # Sf [g/L]
#mu_set = 0 # mu_set [1/h]
#V_fixed = True
#X_fixed = False

pabs = 1 # pabs [bar]
yO = 0.2095 # yO [-]
Sf = 0 # Sf [g/L]
mu_set = 0 # mu_set [1/h] #0.032
V_fixed = True
X_fixed = False
pulse_cycle_time = 10*60 # [s] 180
pulse_on_ratio = 0.38 # [-] 0.23
#D_set = mu_set/im_mmodel.Ysx_ox*X/Sf
##Values from compartment model
D_set = 0 #0.01
perc_starv_comp=  0.5926  # old value0.5556
anae_ratio_comp= 0.3084# old value 0.4188

# im_bmodel_odes(self,t,y,kla,pabs,yO,Sf:float,mu_set:float,V_fixed:bool=False,X_fixed:bool=False,pulse_cycle_time:float=None,pulse_on_ratio:float=None,returns:str='dydt'):
#args_batch=(kla,pabs,yO,Sf,mu_set,V_fixed,X_fixed)
args_batch=(kla,pabs,yO,method,pulse_cycle_time,pulse_on_ratio,anae_ratio_comp,Sf,mu_set,D_set,V_fixed,X_fixed,'dydt') 

sol_batch = odeint(func=func,y0=y0,tfirst=True,t=t_eval,args=args_batch)

Tmpnts=pd.read_csv('C:/Users/s210212/Documents/DTU/Thesis/Lab/E5/tmpnts.csv')
glc_conc=pd.read_csv('C:/Users/s210212/Documents/DTU/Thesis/Lab/E5/glucose_conc.csv')
X_conc=pd.read_csv('C:/Users/s210212/Documents/DTU/Thesis/Lab/E5/X_conc.csv')
Ac_conc=pd.read_csv('C:/Users/s210212/Documents/DTU/Thesis/Lab/E5/Ac_conc.csv')
DO_conc=pd.read_csv('C:/Users/s210212/Documents/DTU/Thesis/Lab/E5/DO_conc.csv')
# Create the figure and two axes objects
#fig, ax1 = plt.subplots(figsize=(7.4/2.54,6/2.54))
fig, ax1 = plt.subplots(figsize=cm2inch(7.4, 6.5))
#fig = plt.figure()
plt.rcParams['font.size']=7 

# Make a second axes that shares the same x-axis
ax2 = ax1.twinx()
Reactor=0
# Plot the first data series on the first axis
ax1.plot(t_eval, sol_batch[:,0], 'g-', label='X sim')
ax1.plot(t_eval, sol_batch[:,1], 'y-', label='Glc sim')
ax1.plot(t_eval, sol_batch[:,3], 'r-', label='Ac sim')
ax1.plot(Tmpnts.iloc[1:36,Reactor],X_conc.iloc[:,Reactor],'go', label='X exp',markersize=3)
ax1.plot(Tmpnts.iloc[1:36,Reactor],glc_conc.iloc[:,Reactor],'yo', label='Glc exp',markersize=3)
ax1.plot(Tmpnts.iloc[1:36,Reactor],Ac_conc.iloc[:,Reactor],'ro', label='Ac exp',markersize=3)
ax1.set_xlabel('Time (h)')
ax1.set_ylabel('c (g/L)')


# Plot the second data series on the second axis
ax2.plot(t_eval, sol_batch[:,2], 'b-', label='O2 sim')
ax2.plot(Tmpnts.iloc[1:36,Reactor],DO_conc.iloc[:,Reactor],'bo', label='O2 exp',markersize=3)
ax2.set_ylabel('c (mmol/L)')


# Add a legend
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
lines = lines1 + lines2
labels = labels1 + labels2
#ax1.legend(lines, labels, loc='center left')
#plt.savefig('C:/Users/s210212/Documents/DTU/Thesis/Report/batch_Exp_simu_WT_Xu.jpeg',dpi=300)
# Show the plot
#plt.figure(figsize=(10/2.54,6/2.54))
plt.savefig('C:/Users/s210212/Documents/DTU/Thesis/Report/Figures/batch_Exp_simu_WT_Xu.jpeg',dpi=600,bbox_inches='tight')
plt.show()

#% Rates
rates_batch = im_mmodel.im_bmodel_odes(
    t=t_eval,
    y=sol_batch.T,
    kla=kla,
    #kla=kla,
    pabs=pabs,
    yO=yO,
    method=method,
    pulse_cycle_time=pulse_cycle_time,
    pulse_on_ratio=pulse_on_ratio,
    anae_ratio_comp=anae_ratio_comp,
    Sf=Sf,
    mu_set=mu_set,
    D_set=D_set,
    V_fixed=V_fixed,
    X_fixed=X_fixed,
    returns='rates'
    )
#rates_batch = im_mmodel.im_bmodel_odes(
#    t=t_eval,
#    y=sol_batch.T,
#    kla=kla,
#    pabs=pabs,
#    yO=yO,
#    Sf=Sf,
#    mu_set=mu_set,
#    V_fixed=V_fixed,
#    X_fixed=X_fixed,
#    returns='rates'
#)
#%%
## HMP strain
method=0 # this part is for the pulse design
im_mmodel = im_xu_bb(
    strain_id='HMP3071',
    strain_description= 'E. coli TRP. prod. strain'
    )
'''
im_mmodel.define_strain_params(
    qS_max=np.mean([1.0346, 1.0549]), #np.mean([1.0536,1.0631])  1.3 ###real np.mean([1.0346, 1.0549])
    qm_max=0, #0.04
    qA_c_max=np.mean([0.05848,0.061354]), #0.15 np.mean([0.05848,0.061354])
    qO_max=np.mean([9.0033,8.8746]), #15np.mean([7.4722,7.292]) #### real np.mean([9.0033,8.8746])
    Ysx_ox=np.mean([0.35331,0.34492]), #0.5np.mean([0.34694,0.34227])     ### real np.mean([0.35331,0.34492])
    Ysx_of=np.mean([0.2745, 0.27167]), #0.15np.mean([0.27253,0.27026])   ### real np.mean([0.2745, 0.27167])
    Ysa=0.667, #np.mean([0.15666,0.16875]), #0.667
    Yax=0, #0.4
    Ki_s=100, #5
    Ks=0.05, #0.05
    Ka=0.05, #0.05
    Ki_o=100 #4
)
'''
#New parameters with calibration correction HPLC
im_mmodel.define_strain_params(
    qS_max=np.mean([1.0956, 1.1169]), #np.mean([1.0536,1.0631])  1.3 ###real np.mean([1.0346, 1.0549])
    qm_max=0, #0.04
    #qA_c_max=np.mean([0.05848,0.061354]), #0.15 np.mean([0.05848,0.061354])
    qA_c_max=np.mean([0.05636,0.06379]),
    qO_max=np.mean([9.0033,8.8746]), #15np.mean([7.4722,7.292]) #### real np.mean([9.0033,8.8746])
    Ysx_ox=np.mean([0.33364,0.3258]), #0.5np.mean([0.34694,0.34227])     ### real np.mean([0.35331,0.34492])
    Ysx_of=np.mean([0.25907, 0.25642]), #0.15np.mean([0.27253,0.27026])   ### real np.mean([0.2745, 0.27167])
    Ysa=0.667, #np.mean([0.15666,0.16875]), #0.667
    Yax=0, #0.4
    Ki_s=100, #5
    Ks=0.05, #0.05
    Ka=0.05, #0.05
    Ki_o=100 #4
)

# ODE functions
func = im_mmodel.im_bmodel_odes
# initial values
y0 = [
    0.3*0.036,    # X [g/L] (OD WT_ini=0.072 HMP_ini=0.036 CDW OD factor=0.3)
    20,     # S [g/L]
    0.2095, # O [mol/L]
    0      # A [g/L]
]
# timesteps to evaluate
t_eval = np.linspace(0,35,10001)

# process parameter function arguments
kla = 150 # kla [1/h] # WT:200 HMP:150
#pabs = 1 # pabs [bar]
#yO = 0.2095 # yO [-]
#Sf = 0 # Sf [g/L]
#mu_set = 0 # mu_set [1/h]
#V_fixed = True
#X_fixed = False

pabs = 1 # pabs [bar]
yO = 0.2095 # yO [-]
Sf = 0 # Sf [g/L]
mu_set = 0 # mu_set [1/h] #0.032
V_fixed = True
X_fixed = False
pulse_cycle_time = 10*60 # [s] 180
pulse_on_ratio = 0.38 # [-] 0.23
#D_set = mu_set/im_mmodel.Ysx_ox*X/Sf
##Values from compartment model
D_set = 0 #0.01
perc_starv_comp=  0.5926  # old value0.5556
anae_ratio_comp= 0.3084# old value 0.4188


# im_bmodel_odes(self,t,y,kla,pabs,yO,Sf:float,mu_set:float,V_fixed:bool=False,X_fixed:bool=False,pulse_cycle_time:float=None,pulse_on_ratio:float=None,returns:str='dydt'):
#args_batch=(kla,pabs,yO,Sf,mu_set,V_fixed,X_fixed)
args_batch=(kla,pabs,yO,method,pulse_cycle_time,pulse_on_ratio,anae_ratio_comp,Sf,mu_set,D_set,V_fixed,X_fixed,'dydt') 

sol_batch = odeint(func=func,y0=y0,tfirst=True,t=t_eval,args=args_batch)

Tmpnts=pd.read_csv('C:/Users/s210212/Documents/DTU/Thesis/Lab/E5/tmpnts.csv')
glc_conc=pd.read_csv('C:/Users/s210212/Documents/DTU/Thesis/Lab/E5/glucose_conc.csv')
X_conc=pd.read_csv('C:/Users/s210212/Documents/DTU/Thesis/Lab/E5/X_conc.csv')
Ac_conc=pd.read_csv('C:/Users/s210212/Documents/DTU/Thesis/Lab/E5/Ac_conc.csv')
DO_conc=pd.read_csv('C:/Users/s210212/Documents/DTU/Thesis/Lab/E5/DO_conc.csv')
# Create the figure and two axes objects
fig, ax1 = plt.subplots(figsize=cm2inch(7.4, 6.5))
#fig = plt.figure()
plt.rcParams['font.size']=7 

# Make a second axes that shares the same x-axis
ax2 = ax1.twinx()
Reactor=2
# Plot the first data series on the first axis
ax1.plot(t_eval, sol_batch[:,0], 'g-', label='X sim')
ax1.plot(t_eval, sol_batch[:,1], 'y-', label='Glc sim')
ax1.plot(t_eval, sol_batch[:,3], 'r-', label='Ac sim')
ax1.plot(Tmpnts.iloc[1:36,Reactor],X_conc.iloc[:,Reactor],'go', label='X exp',markersize=3)
ax1.plot(Tmpnts.iloc[1:36,Reactor],glc_conc.iloc[:,Reactor],'yo', label='Glc exp',markersize=3)
ax1.plot(Tmpnts.iloc[1:36,Reactor],Ac_conc.iloc[:,Reactor],'ro', label='Ac exp',markersize=3)
ax1.set_xlabel('Time (h)')
ax1.set_ylabel('c (g/L)')


# Plot the second data series on the second axis
ax2.plot(t_eval, sol_batch[:,2], 'b-', label='O2 sim')
ax2.plot(Tmpnts.iloc[1:36,Reactor],DO_conc.iloc[:,Reactor],'bo', label='O2 exp',markersize=3)
ax2.set_ylabel('c (mmol/L)')


# Add a legend
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
lines = lines1 + lines2
labels = labels1 + labels2
#ax1.legend(lines, labels, loc='center left')
plt.savefig('C:/Users/s210212/Documents/DTU/Thesis/Report/Figures/batch_Exp_simu_HMP_Xu.jpeg',dpi=600,bbox_inches='tight')
# Show the plt
plt.show()
#% Rates
rates_batch = im_mmodel.im_bmodel_odes(
    t=t_eval,
    y=sol_batch.T,
    kla=kla,
    #kla=kla,
    pabs=pabs,
    yO=yO,
    method=method,
    pulse_cycle_time=pulse_cycle_time,
    pulse_on_ratio=pulse_on_ratio,
    anae_ratio_comp=anae_ratio_comp,
    Sf=Sf,
    mu_set=mu_set,
    D_set=D_set,
    V_fixed=V_fixed,
    X_fixed=X_fixed,
    returns='rates'
    )
#rates_batch = im_mmodel.im_bmodel_odes(
#    t=t_eval,
#    y=sol_batch.T,
#    kla=kla,
#    pabs=pabs,
#    yO=yO,
#    Sf=Sf,
#    mu_set=mu_set,
#    V_fixed=V_fixed,
#    X_fixed=X_fixed,
#    returns='rates'
#)
#%%Batch mode only experiments!!!!
method=0
im_mmodel = im_xu_bb(
    strain_id='DDB35',
    strain_description= 'WT E. coli'
    )

## WT strain
'''
im_mmodel.define_strain_params(
    qS_max=np.mean([0.74036, 0.72265]),#np.mean([0.81528,0.7752]),## #1.3 np.mean([0.74036, 0.72265])
    qm_max=0, #0.04
    qA_c_max=np.mean([0.03444,0.036252]), #0.15
    qO_max=np.mean([6.7453,6.9238]), #15 np.mean([5.9305,5.8531]) np.mea(n[5.4322, 5.5012])
    Ysx_ox=np.mean([0.4592,0.4714]), # np.mean([0.4592,0.4714])    np.mean([0.43349, 44498])
    Ysx_of=np.mean([0.4592,0.4714]), # np.mean([0.4592,0.4714])      0.15 np.mean([0.46473,0.50662])-np.mean([0.41704,0.43812]) np.mean([0.41704,0.43812])
    Ysa=0.667, #0.667 np.mean([0.25974,0.23213])
    Yax=0, #0.4
    Ki_s=100, #5
    Ks=0.05, #0.05
    Ka=0.05, #0.05
    Ki_o=100 #4
)
'''

#New parameters with HPLC better calibration
im_mmodel.define_strain_params(
    qS_max=np.mean([0.78433, 0.76553]),#np.mean([0.81528,0.7752]),## #1.3 np.mean([0.74036, 0.72265])
    qm_max=0, #0.04
    #qA_c_max=np.mean([0.03444,0.036252]), #0.15
    qA_c_max=np.mean([0.03444,0.03726]), #0.15
    qO_max=np.mean([6.7453,6.9238]), #15 np.mean([5.9305,5.8531]) np.mea(n[5.4322, 5.5012])
    Ysx_ox=np.mean([0.43349, 0.44498]), # np.mean([0.4592,0.4714])    np.mean([0.43349, 44498])
    Ysx_of=np.mean([0.43349, 0.44498]), # np.mean([0.4592,0.4714])      0.15 np.mean([0.46473,0.50662])-np.mean([0.41704,0.43812]) np.mean([0.41704,0.43812])
    Ysa=0.667, #0.667 np.mean([0.25974,0.23213])
    Yax=0, #0.4
    Ki_s=100, #5
    Ks=0.05, #0.05
    Ka=0.05, #0.05
    Ki_o=100 #4
)


# ODE functions
func = im_mmodel.im_bmodel_odes
# initial values
y0 = [
    0.3*0.072,    # X [g/L] (OD WT_ini=0.072 HMP_ini=0.036 CDW OD factor=0.3)
    20,     # S [g/L]
    0.2095, # O [mol/L]
    0      # A [g/L]
]
# timesteps to evaluate
t_eval = np.linspace(0,35,100001)

# process parameter function arguments
kla = 200 # kla [1/h] # WT:200 HMP:150
#pabs = 1 # pabs [bar]
#yO = 0.2095 # yO [-]
#Sf = 0 # Sf [g/L]
#mu_set = 0 # mu_set [1/h]
#V_fixed = True
#X_fixed = False

pabs = 1 # pabs [bar]
yO = 0.2095 # yO [-]
Sf = 0 # Sf [g/L]
mu_set = 0 # mu_set [1/h] #0.032
V_fixed = True
X_fixed = False
pulse_cycle_time = 10*60 # [s] 180
pulse_on_ratio = 0.38 # [-] 0.23
#D_set = mu_set/im_mmodel.Ysx_ox*X/Sf
##Values from compartment model
D_set = 0 #0.01
perc_starv_comp=  0.5926  # old value0.5556
anae_ratio_comp= 0.3084# old value 0.4188

# im_bmodel_odes(self,t,y,kla,pabs,yO,Sf:float,mu_set:float,V_fixed:bool=False,X_fixed:bool=False,pulse_cycle_time:float=None,pulse_on_ratio:float=None,returns:str='dydt'):
#args_batch=(kla,pabs,yO,Sf,mu_set,V_fixed,X_fixed)
args_batch=(kla,pabs,yO,method,pulse_cycle_time,pulse_on_ratio,anae_ratio_comp,Sf,mu_set,D_set,V_fixed,X_fixed,'dydt') 

sol_batch = odeint(func=func,y0=y0,tfirst=True,t=t_eval,args=args_batch)

Tmpnts=pd.read_csv('C:/Users/s210212/Documents/DTU/Thesis/Lab/E5/tmpnts.csv')
glc_conc=pd.read_csv('C:/Users/s210212/Documents/DTU/Thesis/Lab/E5/glucose_conc.csv')
X_conc=pd.read_csv('C:/Users/s210212/Documents/DTU/Thesis/Lab/E5/X_conc.csv')
Ac_conc=pd.read_csv('C:/Users/s210212/Documents/DTU/Thesis/Lab/E5/Ac_conc.csv')
DO_conc=pd.read_csv('C:/Users/s210212/Documents/DTU/Thesis/Lab/E5/DO_conc.csv')
# Create the figure and two axes objects
#fig, ax1 = plt.subplots(figsize=(7.4/2.54,6/2.54))
fig, ax1 = plt.subplots(figsize=cm2inch(7.4, 6.5))
#fig = plt.figure()
plt.rcParams['font.size']=7 

# Make a second axes that shares the same x-axis
ax2 = ax1.twinx()
Reactor=0
# Plot the first data series on the first axis
#ax1.plot(t_eval, sol_batch[:,0], 'g-', label='X sim')
#ax1.plot(t_eval, sol_batch[:,1], 'y-', label='Glc sim')
#ax1.plot(t_eval, sol_batch[:,3], 'r-', label='Ac sim')
ax1.plot(Tmpnts.iloc[1:36,Reactor],X_conc.iloc[:,Reactor],'go', label='X exp',markersize=3)
ax1.plot(Tmpnts.iloc[1:36,Reactor],glc_conc.iloc[:,Reactor],'yo', label='Glc exp',markersize=3)
ax1.plot(Tmpnts.iloc[1:36,Reactor],Ac_conc.iloc[:,Reactor],'ro', label='Ac exp',markersize=3)
ax1.set_xlabel('t (h)')
ax1.set_ylabel('c (g/L)')

# Plot the second data series on the second axis
#ax2.plot(t_eval, sol_batch[:,2], 'b-', label='O2 sim')
ax2.plot(Tmpnts.iloc[1:36,Reactor],DO_conc.iloc[:,Reactor],'bo', label='O2 exp',markersize=3)
ax2.set_ylabel('c (mmol/L)')

# Add vertical line
plt.axvline(x = 15.53, color = 'k', linewidth=1, linestyle=':')
plt.axvline(x = 17.55, color = 'k', linewidth=1, linestyle=':')

plt.text(12, 0.2, '(I)', fontsize = 8)
plt.text(15.7, 0.2, '(II)', fontsize = 8)
plt.text(18, 0.2, '(III)', fontsize = 8)

# Add a legend
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
lines = lines1 + lines2
labels = labels1 + labels2
#ax1.legend(lines, labels, loc='center left')
#plt.savefig('C:/Users/s210212/Documents/DTU/Thesis/Report/batch_Exp_simu_WT_Xu.jpeg',dpi=300)
# Show the plot
#plt.figure(figsize=(10/2.54,6/2.54))

ax1.set_ylim(0,20)
ax2.set_ylim(0,0.24)
plt.xlim(9,30)
plt.savefig('C:/Users/s210212/Documents/DTU/Thesis/Report/Figures/batch_Exp_WT_Xu.jpeg',dpi=600,bbox_inches='tight')
plt.show()

#% Rates
rates_batch = im_mmodel.im_bmodel_odes(
    t=t_eval,
    y=sol_batch.T,
    kla=kla,
    #kla=kla,
    pabs=pabs,
    yO=yO,
    method=method,
    pulse_cycle_time=pulse_cycle_time,
    pulse_on_ratio=pulse_on_ratio,
    anae_ratio_comp=anae_ratio_comp,
    Sf=Sf,
    mu_set=mu_set,
    D_set=D_set,
    V_fixed=V_fixed,
    X_fixed=X_fixed,
    returns='rates'
    )
#rates_batch = im_mmodel.im_bmodel_odes(
#    t=t_eval,
#    y=sol_batch.T,
#    kla=kla,
#    pabs=pabs,
#    yO=yO,
#    Sf=Sf,
#    mu_set=mu_set,
#    V_fixed=V_fixed,
#    X_fixed=X_fixed,
#    returns='rates'
#)
#%% ONly experimetn!!!!!
## HMP strain
method=0 # this part is for the pulse design
im_mmodel = im_xu_bb(
    strain_id='HMP3071',
    strain_description= 'E. coli TRP. prod. strain'
    )
'''
im_mmodel.define_strain_params(
    qS_max=np.mean([1.0346, 1.0549]), #np.mean([1.0536,1.0631])  1.3 ###real np.mean([1.0346, 1.0549])
    qm_max=0, #0.04
    qA_c_max=np.mean([0.05848,0.061354]), #0.15 np.mean([0.05848,0.061354])
    qO_max=np.mean([9.0033,8.8746]), #15np.mean([7.4722,7.292]) #### real np.mean([9.0033,8.8746])
    Ysx_ox=np.mean([0.35331,0.34492]), #0.5np.mean([0.34694,0.34227])     ### real np.mean([0.35331,0.34492])
    Ysx_of=np.mean([0.2745, 0.27167]), #0.15np.mean([0.27253,0.27026])   ### real np.mean([0.2745, 0.27167])
    Ysa=0.667, #np.mean([0.15666,0.16875]), #0.667
    Yax=0, #0.4
    Ki_s=100, #5
    Ks=0.05, #0.05
    Ka=0.05, #0.05
    Ki_o=100 #4
)
'''
#New parameters with calibration correction HPLC
im_mmodel.define_strain_params(
    qS_max=np.mean([1.0956, 1.1169]), #np.mean([1.0536,1.0631])  1.3 ###real np.mean([1.0346, 1.0549])
    qm_max=0, #0.04
    #qA_c_max=np.mean([0.05848,0.061354]), #0.15 np.mean([0.05848,0.061354])
    qA_c_max=np.mean([0.05636,0.06379]),
    qO_max=np.mean([9.0033,8.8746]), #15np.mean([7.4722,7.292]) #### real np.mean([9.0033,8.8746])
    Ysx_ox=np.mean([0.33364,0.3258]), #0.5np.mean([0.34694,0.34227])     ### real np.mean([0.35331,0.34492])
    Ysx_of=np.mean([0.25907, 0.25642]), #0.15np.mean([0.27253,0.27026])   ### real np.mean([0.2745, 0.27167])
    Ysa=0.667, #np.mean([0.15666,0.16875]), #0.667
    Yax=0, #0.4
    Ki_s=100, #5
    Ks=0.05, #0.05
    Ka=0.05, #0.05
    Ki_o=100 #4
)

# ODE functions
func = im_mmodel.im_bmodel_odes
# initial values
y0 = [
    0.3*0.036,    # X [g/L] (OD WT_ini=0.072 HMP_ini=0.036 CDW OD factor=0.3)
    20,     # S [g/L]
    0.2095, # O [mol/L]
    0      # A [g/L]
]
# timesteps to evaluate
t_eval = np.linspace(0,35,10001)

# process parameter function arguments
kla = 150 # kla [1/h] # WT:200 HMP:150
#pabs = 1 # pabs [bar]
#yO = 0.2095 # yO [-]
#Sf = 0 # Sf [g/L]
#mu_set = 0 # mu_set [1/h]
#V_fixed = True
#X_fixed = False

pabs = 1 # pabs [bar]
yO = 0.2095 # yO [-]
Sf = 0 # Sf [g/L]
mu_set = 0 # mu_set [1/h] #0.032
V_fixed = True
X_fixed = False
pulse_cycle_time = 10*60 # [s] 180
pulse_on_ratio = 0.38 # [-] 0.23
#D_set = mu_set/im_mmodel.Ysx_ox*X/Sf
##Values from compartment model
D_set = 0 #0.01
perc_starv_comp=  0.5926  # old value0.5556
anae_ratio_comp= 0.3084# old value 0.4188


# im_bmodel_odes(self,t,y,kla,pabs,yO,Sf:float,mu_set:float,V_fixed:bool=False,X_fixed:bool=False,pulse_cycle_time:float=None,pulse_on_ratio:float=None,returns:str='dydt'):
#args_batch=(kla,pabs,yO,Sf,mu_set,V_fixed,X_fixed)
args_batch=(kla,pabs,yO,method,pulse_cycle_time,pulse_on_ratio,anae_ratio_comp,Sf,mu_set,D_set,V_fixed,X_fixed,'dydt') 

sol_batch = odeint(func=func,y0=y0,tfirst=True,t=t_eval,args=args_batch)

Tmpnts=pd.read_csv('C:/Users/s210212/Documents/DTU/Thesis/Lab/E5/tmpnts.csv')
glc_conc=pd.read_csv('C:/Users/s210212/Documents/DTU/Thesis/Lab/E5/glucose_conc.csv')
X_conc=pd.read_csv('C:/Users/s210212/Documents/DTU/Thesis/Lab/E5/X_conc.csv')
Ac_conc=pd.read_csv('C:/Users/s210212/Documents/DTU/Thesis/Lab/E5/Ac_conc.csv')
DO_conc=pd.read_csv('C:/Users/s210212/Documents/DTU/Thesis/Lab/E5/DO_conc.csv')
# Create the figure and two axes objects
fig, ax1 = plt.subplots(figsize=cm2inch(7.4, 6.5))
#fig = plt.figure()
plt.rcParams['font.size']=7 

# Make a second axes that shares the same x-axis
ax2 = ax1.twinx()
Reactor=2
# Plot the first data series on the first axis
#ax1.plot(t_eval, sol_batch[:,0], 'g-', label='X sim')
#ax1.plot(t_eval, sol_batch[:,1], 'y-', label='Glc sim')
#ax1.plot(t_eval, sol_batch[:,3], 'r-', label='Ac sim')
ax1.plot(Tmpnts.iloc[1:36,Reactor],X_conc.iloc[:,Reactor],'go', label='X exp',markersize=3)
ax1.plot(Tmpnts.iloc[1:36,Reactor],glc_conc.iloc[:,Reactor],'yo', label='Glc exp',markersize=3)
ax1.plot(Tmpnts.iloc[1:36,Reactor],Ac_conc.iloc[:,Reactor],'ro', label='Ac exp',markersize=3)
ax1.set_xlabel('t (h)')
ax1.set_ylabel('c (g/L)')


# Plot the second data series on the second axis
#ax2.plot(t_eval, sol_batch[:,2], 'b-', label='O2 sim')
ax2.plot(Tmpnts.iloc[1:36,Reactor],DO_conc.iloc[:,Reactor],'bo', label='O2 exp',markersize=3)
ax2.set_ylabel('c (mmol/L)')

# Add vertical line
plt.axvline(x = 16.78, color = 'k', linewidth=1, linestyle=':')
plt.axvline(x = 20.28, color = 'k', linewidth=1, linestyle=':')

plt.text(15, 0.2, '(I)', fontsize = 8)
plt.text(17.7, 0.2, '(II)', fontsize = 8)
plt.text(21, 0.2, '(III)', fontsize = 8)

# Add a legend
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
lines = lines1 + lines2
labels = labels1 + labels2
#ax1.legend(lines, labels, loc='center left')
ax1.set_ylim(0,20)
ax2.set_ylim(0,0.24)
plt.savefig('C:/Users/s210212/Documents/DTU/Thesis/Report/Figures/batch_Exp_HMP_Xu.jpeg',dpi=600,bbox_inches='tight')
# Show the plt
plt.show()
#% Rates
rates_batch = im_mmodel.im_bmodel_odes(
    t=t_eval,
    y=sol_batch.T,
    kla=kla,
    #kla=kla,
    pabs=pabs,
    yO=yO,
    method=method,
    pulse_cycle_time=pulse_cycle_time,
    pulse_on_ratio=pulse_on_ratio,
    anae_ratio_comp=anae_ratio_comp,
    Sf=Sf,
    mu_set=mu_set,
    D_set=D_set,
    V_fixed=V_fixed,
    X_fixed=X_fixed,
    returns='rates'
    )
#rates_batch = im_mmodel.im_bmodel_odes(
#    t=t_eval,
#    y=sol_batch.T,
#    kla=kla,
#    pabs=pabs,
#    yO=yO,
#    Sf=Sf,
#    mu_set=mu_set,
#    V_fixed=V_fixed,
#    X_fixed=X_fixed,
#    returns='rates'
#)
#%%

###############################CHemostat ##################################

# ODE functions
func = im_mmodel.im_bmodel_odes
# initial values
y0 = [
    40,    # X [g/L]
    0.05,     # S [g/L]
    0.2095, # O [mol/L]
    0      # A [g/L]
]
# timesteps to evaluate
t_eval = np.linspace(0,500,10001)

# process parameter function arguments
kla = 2 # kla [1/h]
pabs = 1 # pabs [bar]
yO = 0.2095 # yO [-]
Sf = 20 # Sf [g/L]
mu_set = 0.13 # mu_set [1/h]
V_fixed = True
X_fixed = False

# im_bmodel_odes(self,t,y,kla,pabs,yO,Sf:float,mu_set:float,V_fixed:bool=False,X_fixed:bool=False,pulse_cycle_time:float=None,pulse_on_ratio:float=None,returns:str='dydt'):
args_chemostat=(kla,pabs,yO,Sf,mu_set,V_fixed,X_fixed)

sol_chemostat = odeint(func=func,y0=y0,tfirst=True,t=t_eval,args=args_chemostat)

#%% find critical qS
func = im_mmodel.im_bmodel_odes
# initial values
y0 = [
    40,    # X [g/L]
    0.05,     # S [g/L]
    0.2095, # O [mol/L]
    0      # A [g/L]
]
# timesteps to evaluate
t_eval = np.linspace(0,500,100001)
mu_test=np.linspace(0.01,0.2,50)
qS=np.zeros(len(mu_test))
qS_of=np.zeros(len(mu_test))
for i in range(len(mu_test)):
    # process parameter function arguments
    kla = 200 # kla [1/h]
    pabs = 1 # pabs [bar]
    yO = 0.2095 # yO [-]
    Sf = 20 # Sf [g/L]
    mu_set = mu_test[i] # mu_set [1/h]
    V_fixed = True
    X_fixed = False
    
    
    # im_bmodel_odes(self,t,y,kla,pabs,yO,Sf:float,mu_set:float,V_fixed:bool=False,X_fixed:bool=False,pulse_cycle_time:float=None,pulse_on_ratio:float=None,returns:str='dydt'):
    args_chemostat=(kla,pabs,yO,Sf,mu_set,V_fixed,X_fixed)
    
    sol_chemostat = odeint(func=func,y0=y0,tfirst=True,t=t_eval,args=args_chemostat)
    rates_chemostat = im_mmodel.im_bmodel_odes(
        t=t_eval,
        y=sol_chemostat.T,
        kla=kla,
        pabs=pabs,
        yO=yO,
        Sf=Sf,
        mu_set=mu_set,
        V_fixed=V_fixed,
        X_fixed=X_fixed,
        returns='rates'
    )
    
    qS[i]=rates_chemostat[2][-1]
    qS_of[i]=rates_chemostat[11][-1]





#plt.plot(qS,qS_of,'o')
plt.plot(mu_test,qS_of,'o')
#plt.plot(qS_of,qS,'o')
#%% Rates
rates_chemostat = im_mmodel.im_bmodel_odes(
    t=t_eval,
    y=sol_chemostat.T,
    kla=kla,
    pabs=pabs,
    yO=yO,
    Sf=Sf,
    mu_set=mu_set,
    V_fixed=V_fixed,
    X_fixed=X_fixed,
    returns='rates'
)

#%%
################################ Pulse sequence ###############################
#### FOR AMBR ####
#%%%             WT
method=1 #method meaning if new SD method is used 1 or old one 0
im_mmodel = im_xu_bb(
    strain_id='DDB35',
    strain_description= 'WT E. coli'
    )

## WT strain
'''
im_mmodel.define_strain_params(
    qS_max=np.mean([0.74036, 0.72265]),#np.mean([0.81528,0.7752]),## #1.3 np.mean([0.74036, 0.72265])
    qm_max=0, #0.04
    qA_c_max=np.mean([0.03444,0.036252]), #0.15
    qO_max=np.mean([6.7453,6.9238]), #15 np.mean([5.9305,5.8531]) np.mea(n[5.4322, 5.5012])
    Ysx_ox=np.mean([0.4592,0.4714]), # np.mean([0.4592,0.4714])    np.mean([0.43349, 44498])
    Ysx_of=np.mean([0.4592,0.4714]), # np.mean([0.4592,0.4714])      0.15 np.mean([0.46473,0.50662])-np.mean([0.41704,0.43812]) np.mean([0.41704,0.43812])
    Ysa=0.667, #0.667 np.mean([0.25974,0.23213])
    Yax=0, #0.4
    Ki_s=100, #5
    Ks=0.05, #0.05
    Ka=0.05, #0.05
    Ki_o=100 #4
)
'''

#New parameters with HPLC better calibration
im_mmodel.define_strain_params(
    qS_max=np.mean([0.78433, 0.76553]),#np.mean([0.81528,0.7752]),## #1.3 np.mean([0.74036, 0.72265])
    qm_max=0, #0.04
    #qA_c_max=np.mean([0.03444,0.036252]), #0.15 The one used the first time 
    qA_c_max=np.mean([0.03444,0.03726]), #0.15
    qO_max=np.mean([6.7453,6.9238]), #15 np.mean([5.9305,5.8531]) np.mea(n[5.4322, 5.5012])
    Ysx_ox=np.mean([0.43349, 0.44498]), # np.mean([0.4592,0.4714])    np.mean([0.43349, 44498])
    Ysx_of=np.mean([0.43349, 0.44498]), # np.mean([0.4592,0.4714])      0.15 np.mean([0.46473,0.50662])-np.mean([0.41704,0.43812]) np.mean([0.41704,0.43812])
    Ysa=0.667, #0.667 np.mean([0.25974,0.23213])
    Yax=0, #0.4
    Ki_s=100, #5
    Ks=0.05, #0.05
    Ka=0.05, #0.05
    Ki_o=100 #4
)

#ODE functions
func = im_mmodel.im_bmodel_odes
# initial values
X = 45#45
y0 = [
    X,    # X [g/L]
    0,     # S [g/L]
    0.2095, # O [mol/L]
    0      # A [g/L]
]

# process parameter function arguments
if method==0:
    kla = 565 # kla [1/h] 565
else:
    kla=1000
pabs = 1 # pabs [bar]
yO = 0.2095 # yO [-]
Sf = 500 # Sf [g/L]
mu_set = 0.05 # mu_set [1/h] #0.032
V_fixed = True
X_fixed = True
pulse_cycle_time = 10*60 # [s] 180
pulse_on_ratio = 0.38 # [-] 0.23
#D_set = mu_set/im_mmodel.Ysx_ox*X/Sf
##Values from compartment model
D_set = 0.0086 #0.01
perc_starv_comp=  0.5926  # old value0.5556
anae_ratio_comp= 0.3084# old value 0.4188

t_eval = np.linspace(0,pulse_cycle_time/3600,100001)
#args_pulse=(kla,pabs,yO,Sf,mu_set,V_fixed,X_fixed,D_set,pulse_cycle_time,pulse_on_ratio,anae_ratio_comp,'dydt')
#args_pulse=(kla,pabs,yO,pulse_cycle_time,pulse_on_ratio,anae_ratio_comp,Sf,mu_set,V_fixed,X_fixed,D_set,'dydt')
args_pulse=(kla,pabs,yO,method,pulse_cycle_time,pulse_on_ratio,anae_ratio_comp,Sf,mu_set,D_set,V_fixed,X_fixed,'dydt')
sol_pulse = odeint(func=func,y0=y0,tfirst=True,t=t_eval,args=args_pulse)

#% Rates
if method==0:
    #args_pulse=(kla,pabs,yO,method,pulse_cycle_time,pulse_on_ratio,anae_ratio_comp,Sf,mu_set,D_set,V_fixed,X_fixed,'dydt')   
    rates_pulse = im_mmodel.im_bmodel_odes(
        t=t_eval,
        y=sol_pulse.T,
        kla=kla,
        #kla=kla,
        pabs=pabs,
        yO=yO,
        method=method,
        pulse_cycle_time=pulse_cycle_time,
        pulse_on_ratio=pulse_on_ratio,
        anae_ratio_comp=anae_ratio_comp,
        Sf=Sf,
        mu_set=mu_set,
        D_set=D_set,
        V_fixed=V_fixed,
        X_fixed=X_fixed,
        returns='rates'
        )
elif method==1:
    #args_pulse=(kla,pabs,yO,method,pulse_cycle_time,pulse_on_ratio,anae_ratio_comp,Sf,mu_set,D_set,V_fixed,X_fixed,'dydt')   
    rates_pulse = im_mmodel.im_bmodel_odes(
        t=t_eval,
        y=sol_pulse.T,
        kla=0,
        #kla=kla,
        pabs=pabs,
        yO=yO,
        method=method,
        pulse_cycle_time=pulse_cycle_time,
        pulse_on_ratio=pulse_on_ratio,
        anae_ratio_comp=anae_ratio_comp,
        Sf=Sf,
        mu_set=mu_set,
        D_set=D_set,
        V_fixed=V_fixed,
        X_fixed=X_fixed,
        returns='rates'
        ) 


starvation_limit=0.05
count=0
counto=0
for i in range(len(t_eval)):
    if rates_pulse[2][i]/im_mmodel.qS_max<starvation_limit:
        count=count+1
perc_starv=round(count/len(t_eval),3)

anae_ratio_avg=np.mean(rates_pulse[11])/np.mean(rates_pulse[2])
#% Figure plotted

fig = plt.figure(figsize=cm2inch(7.4, 6.5))
plt.rcParams['font.size']=7 
plt.plot(t_eval*3600/pulse_cycle_time, rates_pulse[2]/im_mmodel.qS_max, label='qS/qS_max')
plt.plot(t_eval*3600/pulse_cycle_time, np.zeros(len(t_eval))+starvation_limit, label='qS/qS_max=0.05')
plt.plot(t_eval*3600/pulse_cycle_time, rates_pulse[11]/im_mmodel.qS_max, label='qS_of/qS_max')
plt.plot(t_eval*3600/pulse_cycle_time, sol_pulse[:, 2], label='O2')

plt.fill_between([0,pulse_on_ratio], 0, 1, facecolor='gray', alpha=0.2)
#plt.text(pulse_on_ratio/2, 0.5, 'Feed pump on', fontsize=12, ha='center', va='center')
if method==1:
    plt.fill_between([pulse_on_ratio-anae_ratio_avg*pulse_on_ratio,pulse_on_ratio], 0, 1, facecolor='red', alpha=0.1)
#plt.text((pulse_on_ratio-anae_ratio_avg*pulse_on_ratio), 0.4, 'O2 off', fontsize=12, ha='left', va='center')

plt.xlabel('t/t_pulse (s/s)')
#plt.legend(loc='center right')
#plt.legend(bbox_to_anchor=(0.42,0.26))
plt.ylabel('rel. value (-)')

# Create the bar plot and add values to bars
x1 = np.array(['A', 'B',])
y1 = np.array([round(perc_starv_comp,2), round(perc_starv,2)])
x2 = np.array(['C', 'D'])
y2 = np.array([round(anae_ratio_comp,2), round(anae_ratio_avg,2)])

ax = fig.add_subplot(3, 2, 2)
ax.bar(x1, y1)
ax.bar(x2, y2)

# Add the values to the bars
for i, v in enumerate(y1):
    ax.text(i-0.45, v-0.13, str(v), color='black', fontsize=7)
for i, v in enumerate(y2):
    ax.text(i+1.55, v+0.01, str(v), color='black', fontsize=7)
if method==0:
    plt.savefig('C:/Users/s210212/Documents/DTU/Thesis/Report/Figures/SD_design/SD_desgin_WT.jpeg',dpi=600,bbox_inches='tight')
elif method==1:
    plt.savefig('C:/Users/s210212/Documents/DTU/Thesis/Report/Figures/SD_design/SD_desgin_WT_new.jpeg',dpi=600,bbox_inches='tight')
#plt.savefig('C:/Users/s210212/Documents/DTU/Thesis/Report/Figures/SD_design/SD_desgin_WT.jpeg',dpi=600,bbox_inches='tight')
plt.show()
#%%%             HMP
## HMP strain
method=1 #method meaning if new SD method is used 1 or old one 0
im_mmodel = im_xu_bb(
    strain_id='HMP3071',
    strain_description= 'E. coli TRP. prod. strain'
    )
'''
im_mmodel.define_strain_params(
    qS_max=np.mean([1.0346, 1.0549]), #np.mean([1.0536,1.0631])  1.3 ###real np.mean([1.0346, 1.0549])
    qm_max=0, #0.04
    qA_c_max=np.mean([0.05848,0.061354]), #0.15 np.mean([0.05848,0.061354])
    qO_max=np.mean([9.0033,8.8746]), #15np.mean([7.4722,7.292]) #### real np.mean([9.0033,8.8746])
    Ysx_ox=np.mean([0.35331,0.34492]), #0.5np.mean([0.34694,0.34227])     ### real np.mean([0.35331,0.34492])
    Ysx_of=np.mean([0.2745, 0.27167]), #0.15np.mean([0.27253,0.27026])   ### real np.mean([0.2745, 0.27167])
    Ysa=0.667, #np.mean([0.15666,0.16875]), #0.667
    Yax=0, #0.4
    Ki_s=100, #5
    Ks=0.05, #0.05
    Ka=0.05, #0.05
    Ki_o=100 #4
)
'''
#New parameters with calibration correction HPLC
im_mmodel.define_strain_params(
    qS_max=np.mean([1.0956, 1.1169]), #np.mean([1.0536,1.0631])  1.3 ###real np.mean([1.0346, 1.0549])
    qm_max=0, #0.04
    #qA_c_max=np.mean([0.05848,0.061354]), #0.15 np.mean([0.05848,0.061354])
    qA_c_max=np.mean([0.05636,0.06379]),
    qO_max=np.mean([9.0033,8.8746]), #15np.mean([7.4722,7.292]) #### real np.mean([9.0033,8.8746])
    Ysx_ox=np.mean([0.33364,0.3258]), #0.5np.mean([0.34694,0.34227])     ### real np.mean([0.35331,0.34492])
    Ysx_of=np.mean([0.25907, 0.25642]), #0.15np.mean([0.27253,0.27026])   ### real np.mean([0.2745, 0.27167])
    Ysa=0.667, #np.mean([0.15666,0.16875]), #0.667
    Yax=0, #0.4
    Ki_s=100, #5
    Ks=0.05, #0.05
    Ka=0.05, #0.05
    Ki_o=100, #4
)
#ODE functions
func = im_mmodel.im_bmodel_odes
# initial values
X = 34#45
y0 = [
    X,    # X [g/L]
    0,     # S [g/L]
    0.2095, # O [mol/L]
    0      # A [g/L]
]

# process parameter function arguments
if method==0:
    kla = 530 # kla [1/h] 530
elif method==1:
    kla=1000
pabs = 1 # pabs [bar]
yO = 0.2095 # yO [-]
Sf = 500 # Sf [g/L]
#mu_set = 0.032 # mu_set [1/h]
V_fixed = True
X_fixed = True
pulse_cycle_time = 10*60 # [s] 180
pulse_on_ratio = 0.36 # [-] 0.23
# D_set = mu_set/im_mmodel.Ysx_ox*X/Sf
##Values from compartment model
D_set = 0.0066 #0.01
perc_starv_comp=0.6296 # 0.6296 olde value as well
anae_ratio_comp=0.3710 # 0.4294 old value

t_eval = np.linspace(0,pulse_cycle_time/3600,100001)

#
args_pulse=(kla,pabs,yO,method,pulse_cycle_time,pulse_on_ratio,anae_ratio_comp,Sf,mu_set,D_set,V_fixed,X_fixed,'dydt') 
sol_pulse = odeint(func=func,y0=y0,tfirst=True,t=t_eval,args=args_pulse)

#% Rates
if method==0:
    #args_pulse=(kla,pabs,yO,method,pulse_cycle_time,pulse_on_ratio,anae_ratio_comp,Sf,mu_set,D_set,V_fixed,X_fixed,'dydt')   
    rates_pulse = im_mmodel.im_bmodel_odes(
        t=t_eval,
        y=sol_pulse.T,
        kla=kla,
        #kla=kla,
        pabs=pabs,
        yO=yO,
        method=method,
        pulse_cycle_time=pulse_cycle_time,
        pulse_on_ratio=pulse_on_ratio,
        anae_ratio_comp=anae_ratio_comp,
        Sf=Sf,
        mu_set=mu_set,
        D_set=D_set,
        V_fixed=V_fixed,
        X_fixed=X_fixed,
        returns='rates'
        )
elif method==1:
    #args_pulse=(kla,pabs,yO,method,pulse_cycle_time,pulse_on_ratio,anae_ratio_comp,Sf,mu_set,D_set,V_fixed,X_fixed,'dydt')   
    rates_pulse = im_mmodel.im_bmodel_odes(
        t=t_eval,
        y=sol_pulse.T,
        kla=0,
        #kla=kla,
        pabs=pabs,
        yO=yO,
        method=method,
        pulse_cycle_time=pulse_cycle_time,
        pulse_on_ratio=pulse_on_ratio,
        anae_ratio_comp=anae_ratio_comp,
        Sf=Sf,
        mu_set=mu_set,
        D_set=D_set,
        V_fixed=V_fixed,
        X_fixed=X_fixed,
        returns='rates'
        )    

starvation_limit=0.05
count=0
counto=0
for i in range(len(t_eval)):
    if rates_pulse[2][i]/im_mmodel.qS_max<starvation_limit:
        count=count+1
perc_starv=round(count/len(t_eval),3)

anae_ratio_avg=np.mean(rates_pulse[11])/np.mean(rates_pulse[2])

#% Figure plotted
fig = plt.figure(figsize=cm2inch(7.4, 6.5))
plt.rcParams['font.size']=7 
plt.plot(t_eval*3600/pulse_cycle_time, rates_pulse[2]/im_mmodel.qS_max, label='qS/qS_max')
plt.plot(t_eval*3600/pulse_cycle_time, np.zeros(len(t_eval))+starvation_limit, label='qS/qS_max=0.05')
plt.plot(t_eval*3600/pulse_cycle_time, rates_pulse[11]/im_mmodel.qS_max, label='qS_of/qS_max')
plt.plot(t_eval*3600/pulse_cycle_time, sol_pulse[:, 2], label='O2')

plt.fill_between([0,pulse_on_ratio], 0, 1, facecolor='gray', alpha=0.2)
if method==1:
    plt.fill_between([pulse_on_ratio-anae_ratio_avg*pulse_on_ratio,pulse_on_ratio], 0, 1, facecolor='red', alpha=0.1)

plt.xlabel('t/t_pulse (s/s)')
#plt.legend(loc='center right')
#plt.legend(bbox_to_anchor=(0.42,0.26))
plt.ylabel('rel. value (-)')

# Create the bar plot and add values to bars
x1 = np.array(['A', 'B',])
y1 = np.array([round(perc_starv_comp,2), round(perc_starv,2)])
x2 = np.array(['C', 'D'])
y2 = np.array([round(anae_ratio_comp,2), round(anae_ratio_avg,2)])

ax = fig.add_subplot(3, 2, 2)
ax.bar(x1, y1)
ax.bar(x2, y2)

# Add the values to the bars
for i, v in enumerate(y1):
    ax.text(i-0.45, v-0.13, str(v), color='black', fontsize=7)
for i, v in enumerate(y2):
    ax.text(i+1.55, v+0.01, str(v), color='black', fontsize=7)
if method==0:
    plt.savefig('C:/Users/s210212/Documents/DTU/Thesis/Report/Figures/SD_design/SD_desgin_TPS.jpeg',dpi=600,bbox_inches='tight')
elif method==1:
    plt.savefig('C:/Users/s210212/Documents/DTU/Thesis/Report/Figures/SD_design/SD_desgin_TPS_new.jpeg',dpi=600,bbox_inches='tight')
plt.show()

#%% 
###To get the dependence of biomass to starvation ratio and anaerobic ratio####

#%%%             WT

im_mmodel = im_xu_bb(
    strain_id='DDB35',
    strain_description= 'WT E. coli'
    )

## WT strain
'''
im_mmodel.define_strain_params(
    qS_max=np.mean([0.74036, 0.72265]),#np.mean([0.81528,0.7752]),## #1.3 np.mean([0.74036, 0.72265])
    qm_max=0, #0.04
    qA_c_max=np.mean([0.03444,0.036252]), #0.15
    qO_max=np.mean([6.7453,6.9238]), #15 np.mean([5.9305,5.8531]) np.mea(n[5.4322, 5.5012])
    Ysx_ox=np.mean([0.4592,0.4714]), # np.mean([0.4592,0.4714])    np.mean([0.43349, 44498])
    Ysx_of=np.mean([0.4592,0.4714]), # np.mean([0.4592,0.4714])      0.15 np.mean([0.46473,0.50662])-np.mean([0.41704,0.43812]) np.mean([0.41704,0.43812])
    Ysa=0.667, #0.667 np.mean([0.25974,0.23213])
    Yax=0, #0.4
    Ki_s=100, #5
    Ks=0.05, #0.05
    Ka=0.05, #0.05
    Ki_o=100 #4
)
'''

#New parameters with HPLC better calibration
im_mmodel.define_strain_params(
    qS_max=np.mean([0.78433, 0.76553]),#np.mean([0.81528,0.7752]),## #1.3 np.mean([0.74036, 0.72265])
    qm_max=0, #0.04
    qA_c_max=np.mean([0.03444,0.036252]), #0.15
    qO_max=np.mean([6.7453,6.9238]), #15 np.mean([5.9305,5.8531]) np.mea(n[5.4322, 5.5012])
    Ysx_ox=np.mean([0.43349, 0.44498]), # np.mean([0.4592,0.4714])    np.mean([0.43349, 44498])
    Ysx_of=np.mean([0.43349, 0.44498]), # np.mean([0.4592,0.4714])      0.15 np.mean([0.46473,0.50662])-np.mean([0.41704,0.43812]) np.mean([0.41704,0.43812])
    Ysa=0.667, #0.667 np.mean([0.25974,0.23213])
    Yax=0, #0.4
    Ki_s=100, #5
    Ks=0.05, #0.05
    Ka=0.05, #0.05
    Ki_o=100 #4
)



#ODE functions

func = im_mmodel.im_bmodel_odes
X_vector=list(range(20,71,5))
X_vector.remove(50)
ratios_starv_anae=pd.DataFrame(columns=['X (g/L)','perc. starv', 'perc. of'])
for j in range(len(X_vector)):
    
    # initial values
    X = X_vector[j]#45
    y0 = [
        X,    # X [g/L]
        0,     # S [g/L]
        0.2095, # O [mol/L]
        0      # A [g/L]
    ]
    
    # process parameter function arguments
    kla = 1000 # kla [1/h] 570
    pabs = 1 # pabs [bar]
    yO = 0.2095 # yO [-]
    Sf = 500 # Sf [g/L]
    mu_set = 0.05 # mu_set [1/h]
    V_fixed = True
    X_fixed = True
    pulse_cycle_time = 10*60 # [s] 180
    pulse_on_ratio = 0.38 # [-] 0.23
    D_set = mu_set/im_mmodel.Ysx_ox*X/Sf
    ##Values from compartment model
    #D_set = 0.0086 #0.01
    perc_starv_comp=0.5556
    anae_ratio_comp=0.4188
    
    t_eval = np.linspace(0,pulse_cycle_time/3600,10001)
    args_pulse=(kla,pabs,yO,Sf,mu_set,V_fixed,X_fixed,D_set,pulse_cycle_time,pulse_on_ratio,'dydt')
    
    sol_pulse = odeint(func=func,y0=y0,tfirst=True,t=t_eval,args=args_pulse)
    
    #% Rates
    rates_pulse = im_mmodel.im_bmodel_odes(
        t=t_eval,
        y=sol_pulse.T,
        kla=kla,
        pabs=pabs,
        yO=yO,
        Sf=Sf,
        mu_set=mu_set,
        V_fixed=V_fixed,
        X_fixed=X_fixed,
        returns='rates'
    )
    
    starvation_limit=0.05
    count=0
    counto=0
    for i in range(len(t_eval)):
        if rates_pulse[2][i]/im_mmodel.qS_max<starvation_limit:
            count=count+1
    perc_starv=round(count/len(t_eval),3)
    
    anae_ratio_avg=np.mean(rates_pulse[11])/np.mean(rates_pulse[2])
    
    ratios_starv_anae.loc[j,'X (g/L)']=X
    ratios_starv_anae.loc[j,'perc. starv']=perc_starv
    ratios_starv_anae.loc[j,'perc. of']=anae_ratio_avg

plt.plot(ratios_starv_anae.loc[:,'X (g/L)'],ratios_starv_anae.loc[:,'perc. starv'],label='perc starv')
plt.plot(ratios_starv_anae.loc[:,'X (g/L)'],ratios_starv_anae.loc[:,'perc. of'],label='perc overflow')

plt.xlabel('X (g/L)')
plt.ylabel('%/100 ')
plt.legend(loc='lower right')

#%%%             HMP
## HMP strain
im_mmodel = im_xu_bb(
    strain_id='HMP3071',
    strain_description= 'E. coli TRP. prod. strain'
    )
im_mmodel.define_strain_params(
    qS_max=np.mean([1.0346, 1.0549]), #np.mean([1.0536,1.0631])  1.3 ###real np.mean([1.0346, 1.0549])
    qm_max=0, #0.04
    qA_c_max=np.mean([0.05848,0.061354]), #0.15 np.mean([0.05848,0.061354])
    qO_max=np.mean([9.0033,8.8746]), #15np.mean([7.4722,7.292]) #### real np.mean([9.0033,8.8746])
    Ysx_ox=np.mean([0.35331,0.34492]), #0.5np.mean([0.34694,0.34227])     ### real np.mean([0.35331,0.34492])
    Ysx_of=np.mean([0.2745, 0.27167]), #0.15np.mean([0.27253,0.27026])   ### real np.mean([0.2745, 0.27167])
    Ysa=0.667, #np.mean([0.15666,0.16875]), #0.667
    Yax=0, #0.4
    Ki_s=100, #5
    Ks=0.05, #0.05
    Ka=0.05, #0.05
    Ki_o=100 #4
)
#ODE functions
func = im_mmodel.im_bmodel_odes
X_vector=list(range(20,71,5))
X_vector.remove(35)
ratios_starv_anae=pd.DataFrame(columns=['X (g/L)','perc. starv', 'perc. of'])
for j in range(len(X_vector)):
        
    
    # initial values
    X = X_vector[j]#45
    y0 = [
        X,    # X [g/L]
        0,     # S [g/L]
        0.2095, # O [mol/L]
        0      # A [g/L]
    ]
    
    # process parameter function arguments
    kla = 460 # kla [1/h] 570
    pabs = 1 # pabs [bar]
    yO = 0.2095 # yO [-]
    Sf = 500 # Sf [g/L]
    #mu_set = 0.032 # mu_set [1/h]
    V_fixed = True
    X_fixed = True
    pulse_cycle_time = 180 # [s] 180
    pulse_on_ratio = 0.3 # [-] 0.23
    # D_set = mu_set/im_mmodel.Ysx_ox*X/Sf
    ##Values from compartment model
    D_set = 0.0066 #0.01
    perc_starv_comp=0.6296
    anae_ratio_comp=0.4294
    
    
    t_eval = np.linspace(0,pulse_cycle_time/3600,180001)
    args_pulse=(kla,pabs,yO,Sf,mu_set,V_fixed,X_fixed,D_set,pulse_cycle_time,pulse_on_ratio,'dydt')
    
    sol_pulse = odeint(func=func,y0=y0,tfirst=True,t=t_eval,args=args_pulse)
    
    #% Rates
    rates_pulse = im_mmodel.im_bmodel_odes(
        t=t_eval,
        y=sol_pulse.T,
        kla=kla,
        pabs=pabs,
        yO=yO,
        Sf=Sf,
        mu_set=mu_set,
        V_fixed=V_fixed,
        X_fixed=X_fixed,
        returns='rates'
    )
    starvation_limit=0.05
    count=0
    counto=0
    for i in range(len(t_eval)):
        if rates_pulse[2][i]/im_mmodel.qS_max<starvation_limit:
            count=count+1
    perc_starv=round(count/len(t_eval),3)
    
    anae_ratio_avg=np.mean(rates_pulse[11])/np.mean(rates_pulse[2])
    
    ratios_starv_anae.loc[j,'X (g/L)']=X
    ratios_starv_anae.loc[j,'perc. starv']=perc_starv
    ratios_starv_anae.loc[j,'perc. of']=anae_ratio_avg
   
plt.plot(ratios_starv_anae.loc[:,'X (g/L)'],ratios_starv_anae.loc[:,'perc. starv'],label='perc starv')
plt.plot(ratios_starv_anae.loc[:,'X (g/L)'],ratios_starv_anae.loc[:,'perc. of'],label='perc overflow')

plt.xlabel('X (g/L)')
plt.ylabel('%/100 ')
plt.legend(loc='lower right')

#%%
################################ Pulse sequence ###############################
#### FOR Chemostat ####
#%%%             WT
im_mmodel = im_xu_bb(
    strain_id='DDB35',
    strain_description= 'WT E. coli'
    )

## WT strain
im_mmodel.define_strain_params(
    qS_max=np.mean([0.74036, 0.72265]),#np.mean([0.81528,0.7752]),## #1.3 np.mean([0.74036, 0.72265])
    qm_max=0, #0.04
    qA_c_max=np.mean([0.03444,0.036252]), #0.15
    qO_max=np.mean([6.7453,6.9238]), #15 np.mean([5.9305,5.8531]) np.mea(n[5.4322, 5.5012])
    Ysx_ox=np.mean([0.4592,0.4714]), # np.mean([0.4592,0.4714])    0.5 np.mean([0.41704,0.43812]) 
    Ysx_of=np.mean([0.4592,0.4714]), # np.mean([0.4592,0.4714])      0.15 np.mean([0.46473,0.50662])-np.mean([0.41704,0.43812]) np.mean([0.41704,0.43812])
    Ysa=0.667, #0.667 np.mean([0.25974,0.23213])
    Yax=0, #0.4
    Ki_s=100, #5
    Ks=0.05, #0.05
    Ka=0.05, #0.05
    Ki_o=100 #4
)
#ODE functions
func = im_mmodel.im_bmodel_odes
# initial values
X = 10#45
y0 = [
    X,    # X [g/L]
    0,     # S [g/L]
    0.2095, # O [mol/L] 0.035
    0      # A [g/L]
]

# process parameter function arguments
kla = 180 # kla [1/h] 570
pabs = 1 # pabs [bar]
yO = 0.2095 # yO [-]  0.035
Sf = 20 # Sf [g/L]
mu_set = 0.032 # mu_set [1/h]
V_fixed = True
X_fixed = True
pulse_cycle_time = 600 # [s] 180
pulse_on_ratio = 0.1 # [-] 0.23
# D_set = mu_set/im_mmodel.Ysx_ox*X/Sf
##Values from compartment model
D_set = 0.1 #0.01
perc_starv_comp=0.5556
anae_ratio_comp=0.4188

t_eval = np.linspace(0,pulse_cycle_time/3600,200001)
args_pulse=(kla,pabs,yO,Sf,mu_set,V_fixed,X_fixed,D_set,pulse_cycle_time,pulse_on_ratio,'dydt')

sol_pulse = odeint(func=func,y0=y0,tfirst=True,t=t_eval,args=args_pulse)

#% Rates
rates_pulse = im_mmodel.im_bmodel_odes(
    t=t_eval,
    y=sol_pulse.T,
    kla=kla,
    pabs=pabs,
    yO=yO,
    Sf=Sf,
    mu_set=mu_set,
    V_fixed=V_fixed,
    X_fixed=X_fixed,
    returns='rates'
)
starvation_limit=0.05
count=0
counto=0
for i in range(len(t_eval)):
    if rates_pulse[2][i]/im_mmodel.qS_max<starvation_limit:
        count=count+1
perc_starv=round(count/len(t_eval),3)

anae_ratio_avg=np.mean(rates_pulse[11])/np.mean(rates_pulse[2])
'''
for i in range(len(t_eval)):
    if rates_pulse[2][i]/im_mmodel.qS_max>anae_ratio_avg:
        counto=counto+1
perc_O2_lim=round(counto/len(t_eval),3)
'''
#% Figure plotted

fig = plt.figure(figsize=(8, 6))

plt.plot(t_eval*3600/pulse_cycle_time, rates_pulse[2]/im_mmodel.qS_max, label='qS/qS_max')
plt.plot(t_eval*3600/pulse_cycle_time, np.zeros(len(t_eval))+starvation_limit, label='qS/qS_max=0.05')
plt.plot(t_eval*3600/pulse_cycle_time, rates_pulse[11]/im_mmodel.qS_max, label='qS_of/qS_max')
plt.plot(t_eval*3600/pulse_cycle_time, sol_pulse[:, 2], label='O2')

plt.xlabel('t/t_pulse (s/s)')
plt.legend(loc='center right')
plt.ylabel('%/100')

# Create the bar plot and add values to bars
x1 = np.array(['A', 'B',])
y1 = np.array([round(perc_starv_comp,3), round(perc_starv,3)])
x2 = np.array(['C', 'D'])
y2 = np.array([round(anae_ratio_comp,3), round(anae_ratio_avg,3)])

ax = fig.add_subplot(3, 2, 2)
ax.bar(x1, y1)
ax.bar(x2, y2)

# Add the values to the bars
for i, v in enumerate(y1):
    ax.text(i-0.32, v-0.13, str(v), color='black', fontsize=12)
for i, v in enumerate(y2):
    ax.text(i+1.7, v+0.01, str(v), color='black', fontsize=12)

plt.show()

#%%%             HMP
## HMP strain
im_mmodel = im_xu_bb(
    strain_id='HMP3071',
    strain_description= 'E. coli TRP. prod. strain'
    )
im_mmodel.define_strain_params(
    qS_max=np.mean([1.0346, 1.0549]), #np.mean([1.0536,1.0631])  1.3 ###real np.mean([1.0346, 1.0549])
    qm_max=0, #0.04
    qA_c_max=np.mean([0.05848,0.061354]), #0.15 np.mean([0.05848,0.061354])
    qO_max=np.mean([9.0033,8.8746]), #15np.mean([7.4722,7.292]) #### real np.mean([9.0033,8.8746])
    Ysx_ox=np.mean([0.35331,0.34492]), #0.5np.mean([0.34694,0.34227])     ### real np.mean([0.35331,0.34492])
    Ysx_of=np.mean([0.2745, 0.27167]), #0.15np.mean([0.27253,0.27026])   ### real np.mean([0.2745, 0.27167])
    Ysa=0.667, #np.mean([0.15666,0.16875]), #0.667
    Yax=0, #0.4
    Ki_s=100, #5
    Ks=0.05, #0.05
    Ka=0.05, #0.05
    Ki_o=100 #4
)
#ODE functions
func = im_mmodel.im_bmodel_odes
# initial values
X = 10#45
y0 = [
    X,    # X [g/L]
    0,     # S [g/L]
    0.032, # O [mol/L]
    0      # A [g/L]
]

# process parameter function arguments
kla = 180 # kla [1/h] 570
pabs = 1 # pabs [bar]
yO = 0.032 # yO [-]
Sf = 50 # Sf [g/L]
mu_set = 0.032 # mu_set [1/h]
V_fixed = True
X_fixed = True
pulse_cycle_time = 180 # [s] 180
pulse_on_ratio = 0.5 # [-] 0.23
# D_set = mu_set/im_mmodel.Ysx_ox*X/Sf
##Values from compartment model
D_set = 0.0066 #0.01
perc_starv_comp=0.6296
anae_ratio_comp=0.4294

t_eval = np.linspace(0,pulse_cycle_time/3600,10001)
args_pulse=(kla,pabs,yO,Sf,mu_set,V_fixed,X_fixed,D_set,pulse_cycle_time,pulse_on_ratio,'dydt')

sol_pulse = odeint(func=func,y0=y0,tfirst=True,t=t_eval,args=args_pulse)

#% Rates
rates_pulse = im_mmodel.im_bmodel_odes(
    t=t_eval,
    y=sol_pulse.T,
    kla=kla,
    pabs=pabs,
    yO=yO,
    Sf=Sf,
    mu_set=mu_set,
    V_fixed=V_fixed,
    X_fixed=X_fixed,
    returns='rates'
)
starvation_limit=0.05
count=0
counto=0
for i in range(len(t_eval)):
    if rates_pulse[2][i]/im_mmodel.qS_max<starvation_limit:
        count=count+1
perc_starv=round(count/len(t_eval),3)

anae_ratio_avg=np.mean(rates_pulse[11])/np.mean(rates_pulse[2])

#% Figure plotted
fig = plt.figure(figsize=(8, 6))

plt.plot(t_eval*3600/pulse_cycle_time, rates_pulse[2]/im_mmodel.qS_max, label='qS/qS_max')
plt.plot(t_eval*3600/pulse_cycle_time, np.zeros(len(t_eval))+starvation_limit, label='qS/qS_max=0.05')
plt.plot(t_eval*3600/pulse_cycle_time, rates_pulse[11]/im_mmodel.qS_max, label='qS_of/qS_max')
plt.plot(t_eval*3600/pulse_cycle_time, sol_pulse[:, 2], label='O2')

plt.xlabel('t/t_pulse (s/s)')
plt.legend(loc='center right')
plt.ylabel('%/100')

# Create the bar plot and add values to bars
x1 = np.array(['A', 'B',])
y1 = np.array([round(perc_starv_comp,3), round(perc_starv,3)])
x2 = np.array(['C', 'D'])
y2 = np.array([round(anae_ratio_comp,3), round(anae_ratio_avg,3)])

ax = fig.add_subplot(3, 2, 2)
ax.bar(x1, y1)
ax.bar(x2, y2)

# Add the values to the bars
for i, v in enumerate(y1):
    ax.text(i-0.32, v-0.13, str(v), color='black', fontsize=12)
for i, v in enumerate(y2):
    ax.text(i+1.7, v+0.01, str(v), color='black', fontsize=12)

plt.show()




#%% with qS_crit tried stuff

'''        
qS_crit=0.4 #gS/gX/h
YSO_FB=14.285 #mmol/gS from DDB_PD_117
qO_crit=qS_crit*YSO_FB
OTR_max=qO_crit*X
rs_Omax=OTR_max/YSO_FB
rs_max=im_mmodel.qS_max*X
O2_critline=rs_Omax/rs_max
#%% Calculate OTR max from simulation 350000 X=45 g/L
anae_ratio_comp=0.42
#OTR=kla*(O*-O)

#rs_omax=anae_ratio_comp*rs_max
#otr_max=rs_omax*YSO_FB
#print(otr_max)
#plt.plot(t_eval*3600,sol_pulse[:,2])
#%%
'''

#%%
################################# Anane MODEL #################################

#%%    WT
method=0
im_mmodel = im_anane_bb(
    strain_id='DDB35',
    strain_description= 'WT E. coli'
    )
## WT strain
im_mmodel.define_strain_params(
    #qS_max=np.mean([0.74036, 0.72265]), #1.37 (0.8)
    qS_max=np.mean([0.78433, 0.76553]),
    qm_max=0, #0.04 (0.0129)
    #qA_max=np.mean([0.03444,0.036252]), #0.15
    qA_max=np.mean([0.03444,0.03726]),
    qO_max=np.mean([6.7453,6.9238]), #15 (0.0068)
    pA_max=0.8, #0.25 0.8
    #Ysx_ox=np.mean([0.4592,0.4714]) ,
    Ysx_ox=np.mean([0.43349, 0.44498]),
    #Ysx_of=np.mean([0.4592,0.4714]) ,
    Ysx_of=np.mean([0.43349, 0.44498]),
    Ysa=0.667, #np.mean([0.25974,0.23213])
    Yax=0,
    Ki_s=2.1231, #(2.1231)  100 inhibition of glucose on acetate uptake.
    Ki_a=100, #(1.2399)     =Ki_s in Xu    inhibition of acetate on glucose uptake
    Ks=0.05, #(0.037)  0.05 
    Ksa=0.05, #(0.0134)  =Ka in Xu
    Kap=0.5052, #0.5052 (0.1) 
    Ko=0.0001, #10 (0.0001)
    Ki_o=100 #none
)

#% Batch Mode
# ODE functions
func = im_mmodel.im_bmodel_odes
# initial values
y0 = [
    0.3*0.072,    # X [g/L] (OD WT_ini=0.072 HMP_ini=0.036 CDW OD factor=0.3)
    20,     # S [g/L]
    0.2095, # O [mol/L]
    0      # A [g/L]
]
# timesteps to evaluate
t_eval = np.linspace(0,35,10001)

# process parameter function arguments
kla = 200 # kla [1/h] # WT:120 HMP:100
pabs = 1 # pabs [bar]
yO = 0.2095 # yO [-]
Sf = 0 # Sf [g/L]
mu_set = 0 # mu_set [1/h]
V_fixed = True
X_fixed = False
###########Allibi
pulse_cycle_time = 10*60 # [s] 180
pulse_on_ratio = 0.38 # [-] 0.23
#D_set = mu_set/im_mmodel.Ysx_ox*X/Sf
##Values from compartment model
D_set = 0 #0.01
perc_starv_comp=  0.5926  # old value0.5556
anae_ratio_comp= 0.3084# old value 0.4188
############
# im_bmodel_odes(self,t,y,kla,pabs,yO,Sf:float,mu_set:float,V_fixed:bool=False,X_fixed:bool=False,pulse_cycle_time:float=None,pulse_on_ratio:float=None,returns:str='dydt'):
#args_batch=(kla,pabs,yO,Sf,mu_set,V_fixed,X_fixed)
args_batch=(kla,pabs,yO,method,pulse_cycle_time,pulse_on_ratio,anae_ratio_comp,Sf,mu_set,D_set,V_fixed,X_fixed,'dydt')

sol_batch = odeint(func=func,y0=y0,tfirst=True,t=t_eval,args=args_batch)

# Create the figure and two axes objects
#fig, ax1 = plt.subplots()
fig, ax1 = plt.subplots(figsize=cm2inch(7.4, 6.5))
#fig = plt.figure()
plt.rcParams['font.size']=7 

# Make a second axes that shares the same x-axis
ax2 = ax1.twinx()
Reactor=0
# Plot the first data series on the first axis
ax1.plot(t_eval, sol_batch[:,0], 'g-', label='X sim')
ax1.plot(t_eval, sol_batch[:,1], 'y-', label='Glc sim')
ax1.plot(t_eval, sol_batch[:,3], 'r-', label='Ac sim')
ax1.plot(Tmpnts.iloc[1:36,Reactor],X_conc.iloc[:,Reactor],'go', label='X exp',markersize=3)
ax1.plot(Tmpnts.iloc[1:36,Reactor],glc_conc.iloc[:,Reactor],'yo', label='Glc exp',markersize=3)
ax1.plot(Tmpnts.iloc[1:36,Reactor],Ac_conc.iloc[:,Reactor],'ro', label='Ac exp',markersize=3)
ax1.set_xlabel('Time (h)')
ax1.set_ylabel('c (g/L)')


# Plot the second data series on the second axis
ax2.plot(t_eval, sol_batch[:,2], 'b-', label='O2 sim')
ax2.plot(Tmpnts.iloc[1:36,Reactor],DO_conc.iloc[:,Reactor],'bo', label='O2 exp',markersize=3)
ax2.set_ylabel('c (mmol/L)')


# Add a legend
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
lines = lines1 + lines2
labels = labels1 + labels2
#ax1.legend(lines, labels, loc='center left')
plt.savefig('C:/Users/s210212/Documents/DTU/Thesis/Report/Figures/batch_Exp_simu_WT_anane.jpeg',dpi=600,bbox_inches='tight')
# Show the plot
plt.show()
#%%  HMP strain

im_mmodel = im_anane_bb(
    strain_id='HMP3071',
    strain_description= 'E. coli TRP. prod. strain'
    )
im_mmodel.define_strain_params(
    qS_max=np.mean([1.0346, 1.0549]), #1.37 (0.8)
    qm_max=0, #0.04 (0.0129)
    qA_max=np.mean([0.05848,0.061354]), #0.15
    qO_max=np.mean([9.0033,8.8746]), #15 (0.0068)
    pA_max=0.7, #0.25
    Ysx_ox=np.mean([0.35331,0.34492]) ,
    Ysx_of=np.mean([0.2745, 0.27167]) ,
    Ysa=0.667, #np.mean([0.25974,0.23213])
    Yax=0,
    Ki_s=2.1231, #(2.1231)  100 inhibition of glucose on acetate uptake.
    Ki_a=100, #(1.2399)     =Ki_s in Xu    inhibition of acetate on glucose uptake
    Ks=0.05, #(0.037)  0.05 
    Ksa=0.05, #(0.0134)  =Ka in Xu
    Kap=0.5052, #0.5052 (0.1) 
    Ko=0.0001, #10 (0.0001)
    Ki_o=100 #none
)


#% Batch Mode
# ODE functions
func = im_mmodel.im_bmodel_odes
# initial values
y0 = [
    0.3*0.036,    # X [g/L] (OD WT_ini=0.072 HMP_ini=0.036 CDW OD factor=0.3)
    20,     # S [g/L]
    0.2095, # O [mol/L]
    0      # A [g/L]
]
# timesteps to evaluate
t_eval = np.linspace(0,35,10001)

# process parameter function arguments
kla = 150 # kla [1/h] # WT:120 HMP:100
pabs = 1 # pabs [bar]
yO = 0.2095 # yO [-]
Sf = 0 # Sf [g/L]
mu_set = 0 # mu_set [1/h]
V_fixed = True
X_fixed = False
###########Allibi
pulse_cycle_time = 10*60 # [s] 180
pulse_on_ratio = 0.38 # [-] 0.23
#D_set = mu_set/im_mmodel.Ysx_ox*X/Sf
##Values from compartment model
D_set = 0 #0.01
perc_starv_comp=  0.5926  # old value0.5556
anae_ratio_comp= 0.3084# old value 0.4188
############

# im_bmodel_odes(self,t,y,kla,pabs,yO,Sf:float,mu_set:float,V_fixed:bool=False,X_fixed:bool=False,pulse_cycle_time:float=None,pulse_on_ratio:float=None,returns:str='dydt'):
args_batch=(kla,pabs,yO,method,pulse_cycle_time,pulse_on_ratio,anae_ratio_comp,Sf,mu_set,D_set,V_fixed,X_fixed,'dydt')

sol_batch = odeint(func=func,y0=y0,tfirst=True,t=t_eval,args=args_batch)

# Create the figure and two axes objects
#fig, ax1 = plt.subplots()
fig, ax1 = plt.subplots(figsize=cm2inch(7.4, 6.5))
#fig = plt.figure()
plt.rcParams['font.size']=7 
# Make a second axes that shares the same x-axis
ax2 = ax1.twinx()
Reactor=2
# Plot the first data series on the first axis
ax1.plot(t_eval, sol_batch[:,0], 'g-', label='Biomass')
ax1.plot(t_eval, sol_batch[:,1], 'y-', label='Glucose')
ax1.plot(t_eval, sol_batch[:,3], 'r-', label='Acetate')
ax1.plot(Tmpnts.iloc[1:36,Reactor],X_conc.iloc[:,Reactor],'go', label='X exp',markersize=3)
ax1.plot(Tmpnts.iloc[1:36,Reactor],glc_conc.iloc[:,Reactor],'yo', label='Glc exp',markersize=3)
ax1.plot(Tmpnts.iloc[1:36,Reactor],Ac_conc.iloc[:,Reactor],'ro', label='Ac exp',markersize=3)
ax1.set_xlabel('Time (h)')
ax1.set_ylabel('c (g/L)')


# Plot the second data series on the second axis
ax2.plot(t_eval, sol_batch[:,2], 'b-', label='Oxygen')
ax2.plot(Tmpnts.iloc[1:36,Reactor],DO_conc.iloc[:,Reactor],'bo', label='DO exp',markersize=3)
ax2.set_ylabel('c (mmol/L)')


# Add a legend
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
lines = lines1 + lines2
labels = labels1 + labels2
#ax1.legend(lines, labels, loc='center left')
plt.savefig('C:/Users/s210212/Documents/DTU/Thesis/Report/Figures/batch_Exp_simu_HMP_anane.jpeg',dpi=600,bbox_inches='tight')
# Show the plot
plt.show()