# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 16:40:12 2023

@author: sarah
"""

#
# Example of how to run the Python code, and access the output
# This case is identical to the default setup of CLASS (the version with interface) 
#

import copy as cp
from pylab import *
from modelMC import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='whitegrid')

""" 
Create empty model_input and set up case
"""
run1input = model_input()

run1input.dt         = 60.       # time step [s]
run1input.runtime    = 12*3600    # total run time [s]

# mixed-layer input
run1input.sw_ml      = True      # mixed-layer model switch
run1input.sw_shearwe = False     # shear growth mixed-layer switch
run1input.sw_fixft   = False     # Fix the free-troposphere switch
run1input.h          = 200.      # initial ABL height [m]
run1input.Ps         = 101300.   # surface pressure [Pa]
run1input.divU       = 0.        # horizontal large-scale divergence of wind [s-1]
run1input.fc         = 1.e-4     # Coriolis parameter [m s-1]

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

run1input.sw_wind    = False     # prognostic wind switch
run1input.u          = 6.        # initial mixed-layer u-wind speed [m s-1]
run1input.du         = 4.        # initial u-wind jump at h [m s-1]
run1input.gammau     = 0.        # free atmosphere u-wind speed lapse rate [s-1]
run1input.advu       = 0.        # advection of u-wind [m s-2]

run1input.v          = -4.0      # initial mixed-layer u-wind speed [m s-1]
run1input.dv         = 4.0       # initial u-wind jump at h [m s-1]
run1input.gammav     = 0.        # free atmosphere v-wind speed lapse rate [s-1]
run1input.advv       = 0.        # advection of v-wind [m s-2]

run1input.sw_sl      = False     # surface layer switch
run1input.ustar      = 0.3       # surface friction velocity [m s-1]
run1input.z0m        = 0.02      # roughness length for momentum [m]
run1input.z0h        = 0.002     # roughness length for scalars [m]

run1input.sw_rad     = False     # radiation switch
run1input.lat        = 51.97     # latitude [deg]
run1input.lon        = -4.93     # longitude [deg]
run1input.doy        = 268.      # day of the year [-]
run1input.tstart     = 6.8       # time of the day [h UTC]
run1input.cc         = 0.0       # cloud cover fraction [-]
run1input.Q          = 400.      # net radiation [W m-2] 
run1input.dFz        = 0.        # cloud top radiative divergence [W m-2] 

run1input.sw_ls      = False     # land surface switch
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

run1input.Ts         = 290.      # initial surface temperature [K]

run1input.Wmax       = 0.0002    # thickness of water layer on wet vegetation [m]
run1input.Wl         = 0.0000    # equivalent water layer depth for wet vegetation [m]

run1input.Lambda     = 5.9       # thermal diffusivity skin layer [-]

run1input.c3c4       = 'c3'      # Plant type ('c3' or 'c4')

run1input.sw_cu      = False     # Cumulus parameterization switch
run1input.dz_h       = 150.      # Transition layer thickness [m]

# initialize moving column scheme
run1input.slope      = 2.15      # Mountain slope angle (degrees)
run1input.U          = 1.5       # Velocity of the moving column in x direction [m/s]



"""
Init and run the model
"""

# # run1input.wg         = 0.005      # volumetric water content top soil layer [m3 m-3]
# # run1input.w2         = 0.005      # volumetric water content deeper soil layer [m3 m-3]
# run1input.ls_type    = 'js'      # land-surface parameterization ('js' for Jarvis-Stewart or 'ags' for A-Gs)
# run1input.sw_sl      = False     # surface layer switch
# run1input.z0m        = 0.0      # roughness length for momentum [m]
# run1input.z0h        = 0.0     # roughness length for scalars [m]
# run1input.sw_cu = True

run1input.sw_ml      = True      # mixed-layer model switch
run1input.sw_shearwe = True     # shear growth mixed-layer switch
run1input.sw_fixft   = False     # Fix the free-troposphere switch
run1input.sw_wind    = False
run1input.sw_sl      = True       # surface layer switch
run1input.sw_rad     = True     # radiation switch
run1input.sw_ls      = True     # land surface switch --> makes everything super sensitive to initial theta.
run1input.sw_cu      = False    # Cumulus parameterization switch
# run1input.sw_mc      = False    # Moving column switch
run1input.sw_tech    = False    # Evaporation technology switch. Part of the land surface module.
<<<<<<< Updated upstream
run1input.sw_test = False 
=======
<<<<<<< Updated upstream
=======
run1input.sw_test = True 
>>>>>>> Stashed changes
>>>>>>> Stashed changes

run1input.lat        = 40
run1input.lon        = 0
run1input.doy        = 106 # mid April
run1input.tstart     = 6

run1input.alpha = 0.33
run1input.wg = 0.255
run1input.w2 = run1input.wg



run1input.theta = 294
run1input.gammatheta = 4.5e-3 
run1input.gammaq = -1.5e-6
run1input.q = 8.5e-3
run1input.h = 600
run1input.slope = np.rad2deg(np.arctan(1.3/55))
run1input.U          = 3.3
# run1input.runtime    = 55e3 /run1input.U
run1input.runtime = -17741.9 * np.log(3548.4/58548.4)


run1input.sw_mc      = True

r1 = model(run1input)
r1.run()


<<<<<<< Updated upstream

run2input = cp.deepcopy(run1input)

run2input.sw_tech    = True
run2input.sw_ls      = False  # needs to be false for technology to work.
run2input.tech_cov   = 1. 
run2input.rstech     = 0. # [s m-1]
run2input.alpha      = 0.1
run2input.dt         = 60.
run2input.tech_cutoff = 10e3



=======
<<<<<<< Updated upstream
run2input.sw_mc      = True
=======
# run2input.sw_tech    = True
# run2input.sw_ls      = False  # needs to be false for technology to work.
# run2input.tech_cov   = 1. 
# run2input.rstech     = 0. # [s m-1]
# run2input.alpha      = 0.1
# run2input.dt         = 60.
# run2input.tech_cutoff = 10e3



>>>>>>> Stashed changes
>>>>>>> Stashed changes


r2 = model(run2input)
# r2.run()
        
# run3input = cp.deepcopy(run2input)
# run3input.sw_test = True
# r3 = model(run3input)
# r3.run()



if __name__ == "__main__":
    """
    Plot output
    """
    def find_nearest(array, value):
        """
        Function to find the index of the value in an array nearest to given value
        By Gijs van Leeuwen https://github.com/gjvanleeuwen

        Parameters
        ----------
        array : 1d.array
            1d np array wih random values, floats or ints
        value : float or int
            value to look for in array

        Returns
        -------
        idx : int
            index for nearest value

        """
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx

    
<<<<<<< Updated upstream
=======
<<<<<<< Updated upstream
    ax[0,0].plot(r2.out.t, r2.out.LE)
    ax[0,1].plot(r2.out.t, r2.out.thetasurf)
    ax[1,0].plot(r2.out.t, r2.out.q*1000)
    ax[1,1].plot(r2.out.t, r2.out.h)
    ax[2,0].plot(r2.out.t, r2.out.RH_h)
    ax[2,1].plot(r2.out.t, r2.out.H)
=======
>>>>>>> Stashed changes
    fig, ax = plt.subplots(4,3, figsize=(10,10), dpi=300)
    # ax[0,0].plot(r1.out.t[1:], r1.out.LE[1:])
    # ax[0,1].plot(r1.out.t[2:], r1.out.H[2:])
    # ax[1,0].plot(r1.out.t, r1.out.q*1000)
    # ax[1,1].plot(r1.out.t, r1.out.theta)
    # ax[2,0].plot(r1.out.t, r1.out.RH_h)
    # ax[2,1].plot(r1.out.t[3:], r1.out.T2m[3:])
    # ax[0,2].plot(r1.out.t, r1.out.altitude)
    # ax[1,2].plot(r1.out.t, r1.out.h)
    # ax[1,2].plot(r1.out.t, r1.out.zlcl, linestyle='--', c='tab:blue')
    # ax[2,2].plot(r1.out.t, r1.out.dtheta)
    
    n_runs = 5
    colors=plt.cm.nipy_spectral(np.linspace(0,1,n_runs+2))
    
    df = pd.DataFrame(columns=['H2Otech', 'H2Onotech', 'diff'])
    
    times = np.linspace(5,9,n_runs)
    
    for i, t in enumerate(times):
        run2input.tstart = t
        r2 = model(run2input)
        r2.run()
        
        H2Otech = np.sum(r2.out.LE*run2input.dt/(2257*1000*1000))*1000
        
        run1input.tstart = t
        r1 = model(run1input)
        r1.run()
        
        H2Onotech = np.sum(r1.out.LE*run1input.dt/(2257*1000*1000))*1000
        
        df.loc[i] = [H2Otech, H2Onotech, H2Otech - H2Onotech]
        
        
        i += 1

        # ax[0,0].plot(r2.out.t[1:], r2.out.LE[1:])
        # ax[0,1].plot(r2.out.t[2:], r2.out.H[2:])
        # ax[1,0].plot(r2.out.t, r2.out.q*1000)
        # ax[1,1].plot(r2.out.t, r2.out.theta)
        # ax[2,0].plot(r2.out.t, r2.out.RH_h)
        # ax[2,1].plot(r2.out.t[3:], r2.out.T2m[3:])
        # ax[0,2].plot(r2.out.t, r2.out.altitude)
        # ax[1,2].plot(r2.out.t, r2.out.h)
        # ax[1,2].plot(r2.out.t, r2.out.zlcl, linestyle='--', c='tab:orange')
        # ax[2,2].plot(r2.out.t, r2.out.dtheta)
        
        ax[0,0].plot(r2.out.t[1:], r2.out.LE[1:], c=colors[i])
        ax[0,1].plot(r2.out.t[2:], r2.out.H[2:], c=colors[i])
        ax[1,0].plot(r2.out.x/1000, r2.out.q*1000, c=colors[i])
        ax[1,1].plot(r2.out.x/1000, r2.out.theta, c=colors[i])
        ax[2,0].plot(r2.out.x/1000, r2.out.RH_h, c=colors[i])
        # ax[2,0].plot(r2.out.t, r2.out.RH_h, c=colors[i])
        ax[2,1].plot(r2.out.x[3:]/1000, r2.out.T2m[3:], c=colors[i])
        ax[0,2].plot(r2.out.t, r2.out.Q, c=colors[i])
<<<<<<< Updated upstream
        ax[1,2].plot(r2.out.x/1000, r2.out.h, c=colors[i], label=np.round(t,2))
=======
        ax[1,2].plot(r2.out.x/1000, r2.out.h - r2.out.altitude, c=colors[i], label=np.round(t,2))
>>>>>>> Stashed changes
        # ax[1,2].plot(r2.out.t-t, r2.out.zlcl, linestyle='--', c='tab:orange')
        ax[2,2].plot(r2.out.x/1000, r2.out.dtheta, c=colors[i])



        # ax[0,0].plot(r1.out.t[1:], r1.out.LE[1:], c=colors[i], linestyle='--')
        # ax[0,1].plot(r1.out.t[2:], r1.out.H[2:], c=colors[i], linestyle='--')
        # ax[1,0].plot(r1.out.x/1000, r1.out.q*1000, c=colors[i], linestyle='--')
        # ax[1,1].plot(r1.out.x/1000, r1.out.theta, c=colors[i], linestyle='--')
        # ax[2,0].plot(r1.out.x/1000, r1.out.RH_h, c=colors[i], linestyle='--')
        # # ax[2,0].plot(r1.out.t, r1.out.RH_h, c=colors[i], linestyle='--')
        # ax[2,1].plot(r1.out.x[3:]/1000, r1.out.T2m[3:], c=colors[i], linestyle='--')
        # ax[0,2].plot(r1.out.t, r1.out.Q, c=colors[i], linestyle='--')
        # ax[1,2].plot(r1.out.x/1000, r1.out.h, c=colors[i], linestyle='--')#, label=np.round(t,2))
        # # ax[1,2].plot(r1.out.t-t, r1.out.zlcl, linestyle='--', c='tab:orange')
        # ax[2,2].plot(r1.out.x/1000, r1.out.dtheta, c=colors[i], linestyle='--')  
        
        # t12 = find_nearest(r1.out.t, 12)
        # ax[3,2].scatter(r1.out.x[t12]/1000, r1.out.H[t12], color='r')
        # ax[3,2].scatter(r1.out.x[t12]/1000, r1.out.LE[t12], color='b')
    # ax[0,0].plot(r3.out.t[1:], r3.out.LE[1:])
    # ax[0,1].plot(r3.out.t[2:], r3.out.H[2:])
    # ax[1,0].plot(r3.out.t, r3.out.q*1000)
    # ax[1,1].plot(r3.out.t, r3.out.theta)
    # ax[2,0].plot(r3.out.t, r3.out.RH_h)
    # ax[2,1].plot(r3.out.t[3:], r3.out.T2m[3:])
    # ax[0,2].plot(r3.out.t, r3.out.altitude)
    # ax[1,2].plot(r3.out.t, r3.out.h)
    # ax[1,2].plot(r3.out.t, r3.out.zlcl, linestyle='--', c='tab:green')
    # ax[2,2].plot(r3.out.t, r3.out.dtheta)

<<<<<<< Updated upstream
=======
>>>>>>> Stashed changes
>>>>>>> Stashed changes
    
    ax[0,0].set_ylabel('LE')
    ax[0,1].set_ylabel('H')
    ax[1,0].set_ylabel('q*1000')
    ax[1,1].set_ylabel(r'$\theta$')
    ax[2,0].set_ylabel('RH(z=h)')
    ax[2,1].set_ylabel(r'T$_{2m}$')
    ax[0,2].set_ylabel('Q')
    ax[1,2].set_ylabel('ABL top')
    ax[2,2].set_ylabel(r'$\Delta \theta$')
    
    for axis in ax:
        for i in range(len(axis)):
            axis[i].set_xlabel('x (km)')
    ax[0,0].set_xlabel('time (hour)')
    ax[0,1].set_xlabel('time (hour)')
    ax[0,2].set_xlabel('time (hour)')
    
    
    ax[3,0].scatter(times, df['diff'].values)
    ax[3,0].set_xlabel('column start time (hour)')
    ax[3,0].set_ylabel('Extra evap. [mm]')
    dt = np.mean(times[1:] - times[:-1])
    ax[3,0].text(x=times[0], y=df['diff'].values[-2], 
                      s=f"{np.round(np.sum(df['diff'].values)*dt, 2)} mm per day")
    
    # ax[3,1].plot(r1.out.x, r1.out.Ps)
    # ax[3,1].plot(r1.out.x, r1.out.P_h)
    
    ax[3,2].plot(r1.out.x, r1.out.altitude)
    ax[3,2].set_ylabel('altitude')

    
    ax[1,2].legend()
    
    fig.tight_layout()
    # plt.show()
    
    
    
    
    # print(r1.out.theta.max() - r1.out.theta.min())
    # print(r1.out.theta[-1] - r1.out.theta[0])
    # print(r1.out.q.max()*1000)
    # print(r1.out.q[-1]*1000)
    # print(np.sum(r1.out.wq)*run1input.dt)
    

    # subplot(234)
    # plot(r1.out.t, )
    
    