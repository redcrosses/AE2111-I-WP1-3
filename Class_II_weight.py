#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 14:44:08 2024

@author: shiyu
"""
import numpy as np
import math
#from unit_conversion import *
#---------------------------------------------------------------------------------------------------------
def convert_units(value, unit, to_metric=True):
    metric_to_imperial = {
        'cm': value / 2.54,
        'meters_inches': value / 0.0254,  # meters to inches
        'm': value / 0.3048,
        'km': value / 1.60934,
        'g': value / 0.453592,
        'kg': value / 0.45359237,
        'L': value / 3.78541,
        'mL': value / 29.5735,
        'cm^2': value / 6.4516,
        'm^2': value / 0.092903,
        'cm^3': value / 16.3871,
        'm^3': value / 0.0283168,
        'newtons': value / 4.44822,  # N to lbf
        'pascals': value / 47.8803,  # Pa to psf
        'pascals_psi': value / 6894.76,  # Pa to psi
        'newton_meters': value / 1.35582,  # Nm to lb-ft
    }

    imperial_to_metric = {
        'in': value * 2.54,
        'inches_meters': value * 0.0254,  # inches to meters
        'ft': value * 0.3048,
        'yds': value * 0.9144,
        'miles': value * 1.60934,
        'pounds': value * 0.453592,
        'ounces': value * 28.3495,
        'gallons': value * 3.78541,
        'fluid_ounces': value * 29.5735,
        'square_inches': value * 6.4516,
        'ft^2': value * 0.092903,
        'square_yards': value * 0.836127,
        'cubic_inches': value * 16.3871,
        'ft^3': value * 0.0283168,
        'pounds_force': value * 4.44822,  # lbf to N
        'psf': value * 47.8803,  # psf to Pa
        'psi': value * 6894.76,  # psi to Pa
        'pound_feet': value * 1.35582,  # lb-ft to Nm
    }

    if to_metric:
        return imperial_to_metric.get(unit, "Unknown unit")
    else:
        return metric_to_imperial.get(unit, "Unknown unit")

#-----------------------------------------------------------------------------------------------------------

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
    W_press = convert_units(W_press, 'pascals_psi', False)
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
    W_tot = convert_units(W_tot, 'pounds_force', True)
    W_wing = convert_units(W_wing, 'pounds_force', True)
    W_fus = convert_units(W_fus, 'pounds_force', True)
    return W_tot, W_wing, W_fus



#nr_passenger = 11
#nr_pilot = 3
#misc_weight = (32*nr_passenger) + (nr_pilot*60)



