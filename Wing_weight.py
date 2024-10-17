#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 14:44:08 2024

@author: shiyu
"""
import numpy as np
import math
from unit_conversion import convert_units

def Weight_wing(S_w,W_fw,Aspect_ratio,sweep,q,taper,tc,N_z,W_dg):
    return 0.036*(S_w**0.758)*(W_fw**0.0035)*(((Aspect_ratio/(np.cos(sweep)**2)))**0.6)*(q**0.006)*(taper**0.04)*(((100*tc)/np.cos(sweep))**-0.3)*((N_z*W_dg)**0.49)
    

def Weight_fuselage(S_f, N_z, W_dg, L_t, Lift_drag, q, W_press):
    return 0.052*(S_f**1.086)*((N_z*W_dg)**0.177)*(L_t**-0.051)*(Lift_drag**-0.0072)*(q**0.241)+W_press


# Horizontal tail weight calculation
def horizontal_tail_weight(N_z, W_dg, q, S_ht, tc, sweep, Aspect_ratio, sweep_ht, taper_h):
    W_ht = 0.016*((N_z*W_dg)**0.414)*(q**0.168)*(S_ht**0.896)*(((100*tc)/math.cos(sweep))**(-0.12))*((Aspect_ratio/math.cos(sweep_ht)** 2)**0.043)*(taper_h**(-0.02))
    return W_ht

def vertical_tail_weight(N_z, W_dg, q, S_vt, tc, sweep_vt, Aspect_ratio, H_t_H_v, taper_vt):
    W_vt = 0.073*(1+0.2*(H_t_H_v))*((N_z*W_dg)**0.376)*(q**0.122)*(S_vt**0.873)*((100*tc)/math.cos(sweep_vt))**(-0.49)*((Aspect_ratio/math.cos(sweep_vt)**2)**0.357)*(taper_vt**0.039)
    return W_vt

def Weight_LandingGear(Nl,Wl,Lm,Ln):
    
    Wmlg=0.095*((Nl*Wl)**0.768)*((Lm/12)**0.409)
    Wnlg=0.125*((Nl*Wl)**0.566)*((Ln/12)**0.845)
    
    return Wmlg+Wnlg

def class_II_weight(S_w, W_fw, Aspect_ratio, sweep, q, taper, tc, N_z, W_dg, S_f, L_t, Lift_drag, W_press, S_ht, sweep_ht, taper_h, S_vt, sweep_vt, H_t_H_v, taper_vt, Nl, Wl, Lm, Ln):
    #unit convert
    S_w = convert_units(S_w, 'm^2', False)
    W_fw = convert_units(W_fw, 'newtons', False)
    q = convert_units(q, 'pascals', False)
    W_dg = convert_units(W_dg, 'newtons', False)
    S_f = convert_units(S_f, 'm^2', False)
    L_t = convert_units(L_t, 'm', False)
    W_press = convert_units(W_press, 'pascal_psi', False)
    S_ht = convert_units(S_ht, 'm^2', False)
    S_vt = convert_units(S_vt, 'm^2', False)
    Nl = convert_units(S_vt, 'm^2', False)
    Wl = convert_units(S_vt, 'm^2', False)
    Lm = convert_units(S_vt, 'meters_inches', False)
    Ln = convert_units(S_vt, 'meters_inches', False)

    
    #weight estimate
    W_wing = Weight_wing(S_w,W_fw,Aspect_ratio,sweep,q,taper,tc,N_z,W_dg)
    W_fus = Weight_fuselage(S_f, N_z, W_dg, L_t, Lift_drag, q, W_press)
    W_hor_tail = horizontal_tail_weight(N_z, W_dg, q, S_ht, tc, sweep, Aspect_ratio, sweep_ht, taper_h)
    W_vert_tail = vertical_tail_weight(N_z, W_dg, q, S_vt, tc, sweep_vt, Aspect_ratio, H_t_H_v, taper_vt)
    W_land = Weight_LandingGear(Nl,Wl,Lm,Ln)
    
    powerplant_mass = 2 * 6078
    W_power = powerplant_mass*32.1740 #ft/s^2
    
    W_tot = W_wing + W_fus + W_hor_tail + W_vert_tail + W_land + W_power
    
    #convert back
    return W_tot



#nr_passenger = 11
#nr_pilot = 3
#misc_weight = (32*nr_passenger) + (nr_pilot*60)



