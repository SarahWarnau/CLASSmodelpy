# -*- coding: utf-8 -*-
"""
Created on Thu Feb  5 09:51:01 2026

@author: swar
"""

import numpy as np
import metpy.calc as mpcalc
from metpy.units import units

class SprayModule:
    def __init__(self, theta0, q0, p):
        
        print(f"before spray theta={np.round(theta0,2)} K, q={np.round(q0*1e3,2)} g/kg")
        
        # constants as in MicroHH paper (van Heerwaarden et al. 2017)
        self.kappa = 0.4 # Von Karman constant [-]
        self.g = 9.81 # gravitational acceleration [m s-2]
        self.cp = 1005. # specific heat of dry air at constant pressure [J k-1 K-1]
        self.p0 = 1.e5 # reference pressure [Pa]
        self.Rd = 287.04 # gas contant for dry air [J K-1 kg-1]
        self.Rv = 461.5 # gas constant for water vapor [J K-1 kg-1]
        self.Lv = 2.5e6 # latent heat of vaporization [J kg-1]
        self.epsilon = self.Rd/self.Rv
        
        self.theta=theta0
        self.q=q0
        self.p=p
        
        self.theta_sprayed, self.q_sprayed = self.calc_sprayed_parcel()
        
        print(f"after spray theta={np.round(self.theta_sprayed,2)} K, q={np.round(self.q_sprayed*1e3,2)} g/kg")
        
    # atmospheric thermodynamics in mass units
    def calc_esat(self, T):
        return 611.2*np.exp(17.67*(T-273.15)/(T - 273.15 + 243.5))

    def calc_qsat(self, T):
        """output kg/kg"""
        es = self.calc_esat(T)
        rs = self.epsilon*es/(self.p - es)
        return rs/(1+rs)

    def exner(self):
        return ((self.p0/self.p)**(self.Rd/self.cp))

    def calc_T_from_theta(self):
        return self.theta/self.exner()

    def calc_Td(self):
        """input q in kg/kg"""
        r = self.q/(1-self.q)
        e = self.p*r/(r + self.epsilon)
        return mpcalc.dewpoint(vapor_pressure=e*units.Pa).to('K').magnitude

    def calc_Tw(self):
        Td = self.calc_Td()
        T = self.calc_T_from_theta()
        return mpcalc.wet_bulb_temperature(pressure=self.p*units.Pa, 
                                                    temperature=T*units.K, 
                                                    dewpoint=Td*units.K).to('K').magnitude

    def calc_theta_from_T(self, T):
        return T*self.exner()
    
    # calculate new specific humidity after spray evaporation
    def calc_sprayed_parcel(self):
        """output in kg/kg"""
        T_sprayed = self.calc_Tw()
        theta_sprayed = self.calc_theta_from_T(T_sprayed).item()
        q_sprayed = self.calc_qsat(T=T_sprayed).item()
        return theta_sprayed, q_sprayed