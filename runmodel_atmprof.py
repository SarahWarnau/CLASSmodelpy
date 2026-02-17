# -*- coding: utf-8 -*-
"""
Created on Tue Feb  3 15:55:01 2026

@author: swar
"""

from model import *
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

"""
Get pandas dataframe with atmosheric profiles on z levels
"""
dz = 0.1
z = np.arange(0, 750, dz, dtype=np.float64())
df = pd.DataFrame({'z':z,
                   'theta':293+z*6e-3,
                   'q':9e-3-z*1e-6,
                   'u':4 + z*1e-2,
                   'v':z*0}
                  )
# df.loc[:10, 'theta'] = 288
dx = 100
x = np.arange(0, 25e3, dx)
X = np.zeros(len(x))
X[int(15e3/dx):] = 1
X[int(5e3/dx):int(10e3/dx)] = 2

x_spray = 8e3


""" 
Create empty model_input and set up case
"""
run1input = model_input()

run1input.dt         = 60.       # time step [s]
run1input.runtime    = 3*3600    # total run time [s]

# atmospheric profile input
# When active, the provided dataframe will be used to calculate the initial, gradient, 
# and jump values of the variables present. 
# The dataframe should contain a column labeled 'z' and optionally 'theta', 'q', 'u', 'v'.
run1input.sw_ap      = True      # atmospheric profile switch
run1input.ap = df

run1input.sw_x       = True     # switch distance instead of time
run1input.col_vel    = 4.        # column velocity [m s-1]
run1input.x          = x         # numpy array with distances [m], evenly spaced
run1input.X          = X         # numpy array with surface codes (0=sea, 1=land, 2=solar evaporator)

run1input.sw_sp      = True     # spray switch
run1input.tspray     = x_spray/run1input.col_vel      # time of spraying [s]
run1input.zspray     = 100       # height of spraying [m]

run1input.sw_so      = False      # solar evaporation technology switch
run1input.rstech     = 0.        # surface resistance to evaporation of solar evaporator
run1input.z0mtech    = 0.01      # roughness length for momentum solar evaporator [m]
run1input.z0htech    = 0.01      # roughness length for scalars solar evaporator [m]

# mixed-layer input
run1input.sw_ml      = True     # mixed-layer model switch
run1input.sw_shearwe = True     # shear growth mixed-layer switch
run1input.sw_fixft   = True      # Fix the free-troposphere switch
run1input.h          = 200.      # initial ABL height [m]
run1input.Ps         = 101300.   # surface pressure [Pa]
run1input.divU       = 0.        # horizontal large-scale divergence of wind [s-1]
run1input.fc         = 0#1.e-4     # Coriolis parameter [m s-1]

run1input.theta      = 288.      # initial mixed-layer potential temperature [K]
run1input.dtheta     = 1.        # initial temperature jump at h [K]
run1input.gammatheta = 0.006     # free atmosphere potential temperature lapse rate [K m-1]
run1input.advtheta   = 0.        # advection of heat [K s-1]
run1input.beta       = 0.2       # entrainment ratio for virtual heat [-]
run1input.wtheta     = 0.1       # surface kinematic heat flux [K m s-1]

run1input.q          = 0.008     # initial mixed-layer specific humidity [kg kg-1]
run1input.dq         = -0.001    # initial specific humidity jump at h [kg kg-1]
run1input.gammaq     = 0.        # free atmosphere specific humidity lapse rate [kg kg-1 m-1]
run1input.advq       = 0.        # advection of moisture [kg kg-1 s-1]
run1input.wq         = 0.1e-3    # surface kinematic moisture flux [kg kg-1 m s-1]

run1input.CO2        = 422.      # initial mixed-layer CO2 [ppm]
run1input.dCO2       = -44.      # initial CO2 jump at h [ppm]
run1input.gammaCO2   = 0.        # free atmosphere CO2 lapse rate [ppm m-1]
run1input.advCO2     = 0.        # advection of CO2 [ppm s-1]
run1input.wCO2       = 0.        # surface kinematic CO2 flux [ppm m s-1]

run1input.sw_wind    = True     # prognostic wind switch
run1input.u          = 6.        # initial mixed-layer u-wind speed [m s-1]
run1input.du         = 4.        # initial u-wind jump at h [m s-1]
run1input.gammau     = 0.        # free atmosphere u-wind speed lapse rate [s-1]
run1input.advu       = 0.        # advection of u-wind [m s-2]

run1input.v          = -4.0      # initial mixed-layer u-wind speed [m s-1]
run1input.dv         = 4.0       # initial u-wind jump at h [m s-1]
run1input.gammav     = 0.        # free atmosphere v-wind speed lapse rate [s-1]
run1input.advv       = 0.        # advection of v-wind [m s-2]

run1input.sw_sl      = True     # surface layer switch
run1input.ustar      = 0.3       # surface friction velocity [m s-1]
run1input.z0m        = 0.02      # roughness length for momentum [m]
run1input.z0h        = 0.002     # roughness length for scalars [m]

run1input.sw_rad     = False     # radiation switch
run1input.lat        = 30.     # latitude [deg]
run1input.lon        = 0.     # longitude [deg]
run1input.doy        = 180.      # day of the year [-]
run1input.tstart     = 9.        # time of the day [h UTC]
run1input.cc         = 0.0       # cloud cover fraction [-]
run1input.Q          = 600.      # net radiation [W m-2] 
run1input.Qin        = 1000.     # net incoming radiation (Swin + Lwin)
run1input.dFz        = 0.        # cloud top radiative divergence [W m-2] 

run1input.sw_ss      = False     # sea surface switch
run1input.SST        = df['theta'][0]+2 # sea surface temperature
run1input.z0msea    = 1e-4      # roughness length for momentum solar evaporator [m]
run1input.z0hsea    = 1e-4      # roughness length for scalars solar evaporator [m]

run1input.sw_ls      = True     # land surface switch
run1input.ls_type    = 'js'      # land-surface parameterization ('js' for Jarvis-Stewart or 'ags' for A-Gs)
run1input.wg         = 0.21      # volumetric water content top soil layer [m3 m-3]
run1input.w2         = 0.21      # volumetric water content deeper soil layer [m3 m-3]
run1input.cveg       = 0.85      # vegetation fraction [-]
run1input.Tsoil      = 285.      # temperature top soil layer [K]
run1input.T2         = 286.      # temperature deeper soil layer [K]
run1input.a          = 0.219     # Clapp and Hornberger retention curve parameter a
run1input.b          = 4.90      # Clapp and Hornberger retention curve parameter b
run1input.p          = 4.        # Clapp and Hornberger retention curve parameter c
run1input.CGsat      = 3.56e-6   # saturated soil conductivity for heat

run1input.wsat       = 0.472     # saturated volumetric water content ECMWF config [-]
run1input.wfc        = 0.323     # volumetric water content field capacity [-]
run1input.wwilt      = 0.171     # volumetric water content wilting point [-]

run1input.C1sat      = 0.132     
run1input.C2ref      = 1.8

run1input.LAI        = 2.        # leaf area index [-]
run1input.gD         = 0.0       # correction factor transpiration for VPD [-]
run1input.rsmin      = 110.      # minimum resistance transpiration [s m-1]
run1input.rssoilmin  = 50.       # minimun resistance soil evaporation [s m-1]
run1input.alpha      = 0.25      # surface albedo [-]

run1input.Ts         = 300.      # initial surface temperature [K]

run1input.Wmax       = 0.0002    # thickness of water layer on wet vegetation [m]
run1input.Wl         = 0.0000    # equivalent water layer depth for wet vegetation [m]

run1input.Lambda     = 5.9       # thermal diffusivity skin layer [-]

run1input.c3c4       = 'c3'      # Plant type ('c3' or 'c4')

run1input.sw_cu      = False     # Cumulus parameterization switch
run1input.dz_h       = 150.      # Transition layer thickness [m]

"""
Init and run the model
"""
r1 = model(run1input)
r1.run()

"""
Plot output
"""

# plt.plot(r1.out.h)
fig, axs = plt.subplots(2,2)
ax = axs.flatten()
for i, var in enumerate(['theta', 'q', 'u', 'v']):
    mpt = 1
    if var == 'q':
        mpt = 1000
    if run1input.sw_x:
        img = ax[i].contourf(x[:-1]/1000, r1.out_NetCDF['z']/1000, 
                             r1.out_NetCDF[var].T*mpt, levels=20)
    else:
        img = ax[i].contourf(r1.out_NetCDF['time']/3600, r1.out_NetCDF['z']/1000, 
                             r1.out_NetCDF[var].T*mpt, levels=20)
    plt.colorbar(img)
    ax[i].set_title(var)

fig.tight_layout()
plt.show()
# plt.ylim(0,len(df))

variables = {'theta [K]' : 'theta', 
            'gamma theta [K km-1]' : 'gammatheta', 
            'd thetav [K]' : 'dthetav',
            'q [g kg-1]' : 'q',  
            'gamma q [g kg-1 km-1]' : 'gammaq', 
            'd q [g kg-1]': 'dq', 
            'theta surf. [K]':'thetasurf',
            'LE [W m-2]':'LE',
            'H [W m-2]':'H',
            'surface code':'X',
            'ABL (orange=LCL) [m]':'h',
            'Net Radiation (orange=in) [W m-2]':'Q'}

fig, axs = plt.subplots(4,3, figsize=(10,6))
ax = axs.flatten()
for i, key in enumerate(variables):
    var = variables[key]
    mpt = 1
    mpt2 = 1
    if 'q' in var:
        mpt = 1000
    if 'gamma' in var:
        mpt2 = 1000
    if run1input.sw_x:
        img = ax[i].plot(x[:-1]/1000, r1.out.__dict__[var]*mpt*mpt2)
    else:
        img = ax[i].plot(r1.out.__dict__[var]*mpt*mpt2)
    ax[i].set_title(key)

if run1input.sw_x:
    img = ax[-1].plot(x[:-1]/1000, r1.out.__dict__['Qin']*mpt*mpt2)
    img = ax[-2].plot(x[:-1]/1000, r1.out.__dict__['zlcl']*mpt*mpt2)
else:
    img = ax[-1].plot(r1.out.__dict__['Qin'])
    img = ax[-2].plot(r1.out.__dict__['zlcl'])
fig.tight_layout()
plt.show()

#%% Test energy balance with atmospheric profile cases
if run1input.sw_ap:
    thetachange = (r1.out_NetCDF['theta'].sum('z').values[1:] - r1.out_NetCDF['theta'].sum('z').values[:-1])*run1input.col_vel*dz/dx
    qchange = 1e3*(r1.out_NetCDF['q'].sum('z').values[1:] - r1.out_NetCDF['q'].sum('z').values[:-1])*run1input.col_vel*dz/dx
    
    fig, ax = plt.subplots(2, 2, figsize=(9,6))
    ax[0,0].set_title('potential temperature change [K s-1]')
    ax[0,0].plot(r1.out.x[1:], thetachange, label='change from profile output')
    ax[0,0].plot(r1.out.x[1:], r1.out.wtheta[1:], ls='--', label='surface flux')
    ax[0,1].set_title('specific humidity change [g kg-1 s-1]')
    ax[0,1].plot(r1.out.x[1:], qchange)
    ax[0,1].plot(r1.out.x[1:], 1e3*r1.out.wq[1:], ls='--')
    
    if run1input.sw_sp:
        thetachange[np.where(thetachange==min(thetachange))] = np.nan
        qchange[np.where(qchange==max(qchange))] = np.nan
        for axs in ax.flatten():
            axs.axvline(x_spray, zorder=0, linestyle=':', label='moment of spray')
    
    # remove timestep of spray from data
    ax[1,0].set_title('cumulative potential temperature change [K]')
    ax[1,0].plot(r1.out.x[1:], np.nancumsum(thetachange), label='change from profile output')
    ax[1,0].plot(r1.out.x[1:], np.cumsum(r1.out.wtheta[1:]), label='surface flux', ls='--')
    ax[1,1].set_title('cumulative specific humidity change [g kg-1]')
    ax[1,1].plot(r1.out.x[1:], np.nancumsum(qchange))
    ax[1,1].plot(r1.out.x[1:], np.cumsum(1e3*r1.out.wq[1:]), ls='--')
    
    ax[0,0].set_ylim(-0.1, 0.5)
    ax[0,1].set_ylim(-0.1, 0.5)
    ax[0,0].legend()
    fig.suptitle("Energy and water conservation check")
    fig.tight_layout()