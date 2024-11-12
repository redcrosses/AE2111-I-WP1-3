#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 14:05:52 2024

@author: shiyu
"""


from scipy.optimize import root
from scipy.optimize import fsolve
import numpy as np

def fuselage(req_vol_fuel):
        
    seat = 0.46 #estimate value based off reader
    aisle = 0.6 #estimate value based off reader
    aisle_nr = 1 #estimate value based off reader
    armrest = 0.1
    passenger_nr = 86
    clearance = 0.05
    seats_aisle = np.round(0.45 * np.sqrt(passenger_nr))

    total_width =seats_aisle*seat+(seats_aisle+aisle_nr+1)*armrest+aisle*aisle+2*clearance

    legroom = total_width - 2*clearance - 2*armrest

    headroom = legroom - seat

    print("\n\033[1m\033[4m Cross section dimentions [m] \033[0m")

    print("{:24} {:.5f} {:16}\n{:24} {:.5f} {:16}\n{:24} {:.5f} {:16}\n{:24} {:.5f} {:16}\n{:24} {:.5f} {:16}\n{:24} {:.5f} {:16}\n{:24} {:.5f} {:16}\n{:24} {:.5f} {:16}".format("# of Passenger:", passenger_nr, "[-]", "Seat width:", seat, "[m]", "Aisle width:", aisle, "[m]", "Armrest width:", armrest, "[m]", "Clearance:", clearance, "[m]", "Total width:", total_width, "[m]", "Leg room:", legroom, "[m]", "Head room:", headroom, "[m]"))
    #chatgpt goated here^^
    inner_D = 3.12653
    outer_D = 1.045*inner_D + 0.084

    print("Inner diameter:", inner_D, "[m]")
    print("outer diameter:", outer_D, "[m]")
    print()

    #--------------------------------------------------------------

    print("\033[1m\033[4m Outer dimentions [m] \033[0m")
    double_bubble = False
    if not double_bubble:
        w_fus = outer_D #[m]
        h_fus = outer_D #[m]
    else:
        w_fus = 3.35123 #[m]
        h_fus = 4.05625 #[m]


    k_cabin = 1.17 #for long range airplanes

    w = 2.87375 #[m]

    d_fus = (w_fus+h_fus)/2

    l_cabin = k_cabin * passenger_nr/seats_aisle
    l_n = 4.5 #np.mean([3.0,4.5]) #average of range in reader
    l_ncratio = 1.5 #np.mean([1.5,2.5]) #average of range in reader
    l_nc = l_ncratio*d_fus
    l_tcratio = 2 #np.mean([2.0,5.0]) #average of range in reader
    l_tc = l_tcratio*d_fus
    l_tratio = 1 #np.mean([0.5,1.0]) #average of range in reader
    l_t = l_tratio*d_fus
    l_fus = l_cabin + l_n + l_t
    l_cyl = l_fus - l_tc - l_nc

    # req_vol_fuel = 83.1 #m^3

    not_optimal = True

    def f(x):
        value = (2*(x**2) - w**2) / (2*(x**2))
        # Ensure the value is within the valid range for arccos
        value = np.clip(value, -1, 1)
        return np.pi * (x**2) - (1/2) * (x**2) * (np.arccos(value) - np.sin(np.arccos(value))) - area_req

    while not_optimal and double_bubble:
        #1:start with initial radius dimentions
        #2:get cylynder length
        #3:calculate area required
        #4:calculate new dimentions of cross section
        #5:calculate new cylinder length
        #6:check if it meets volume requirement
            #if yes: return value
            #if no: repeat from 3
        area_req = req_vol_fuel/l_cyl
        
        initial_guesses = [1,2]
        solution = fsolve(f,initial_guesses)
        
        r_chin = abs(solution[1])
        
        inner_r = inner_D/2

        in_h_1 = r_chin + r_chin*np.cos(np.arccos(((2*r_chin**2)-w**2)/(2*r_chin**2))/2)
        in_h_2 = inner_r +inner_r*np.cos(np.arccos(((2*inner_r**2)-w**2)/(2*inner_r**2))/2)
        
        h_1 = in_h_1 + ((r_chin*2)*1.045+0.084)/2-r_chin
        h_2 = in_h_2 + ((inner_r*2)*1.045+0.084)/2 - inner_r
        
        h_fus = h_1 + h_2
        
        if r_chin >= inner_r:
            w_fus = 1.045*r_chin*2 + 0.084
        else:
            w_fus = 1.045*inner_r*2 + 0.084

        d_fus = (w_fus+h_fus)/2

        l_cabin = k_cabin * passenger_nr/seats_aisle
        l_n = 4.5 #np.mean([3.0,4.5]) #average of range in reader
        l_ncratio = 1.5 #np.mean([1.5,2.5]) #average of range in reader
        l_nc = l_ncratio*d_fus
        l_tcratio = 2 #np.mean([2.0,5.0]) #average of range in reader
        l_tc = l_tcratio*d_fus
        l_tratio = 1 #np.mean([0.5,1.0]) #average of range in reader
        l_t = l_tratio*d_fus
        l_fus = l_cabin + l_n + l_t
        l_cyl = l_fus - l_tc - l_nc
        
        new_fuel_vol = l_cyl * area_req
        
        if new_fuel_vol == req_vol_fuel:
            not_optimal = False

    S_wfus = (np.pi*w_fus/4)*((1/(3*l_nc**2))*((((4*l_nc**2)+(w_fus**2/4))**1.5)-((l_tc**3)/8))-w_fus+4*l_cyl+2*(np.sqrt((l_tc**2)+((w_fus**2)/4))))

    print("{:24} {:.5f} {:16}\n{:24} {:.5f} {:16}\n{:24} {:.5f} {:16}\n{:24} {:.5f} {:16}\n{:24} {:.5f} {:16}\n{:24} {:.5f} {:16}\n{:24} {:.5f} {:16}\n{:24} {:.5f} {:16}\n{:24} {:.5f} {:16}\n{:24} {:.5f} {:16}".format("h_fus:", h_fus, "[m]", "w_fus:", w_fus, "[m]", "l_cabin:", l_cabin, "[m]", "l_n:", l_n, "[m]", "l_nc:", l_nc, "[m]", "l_tc:", l_tc, "[m]", "l_t:", l_t, "[m]", "l_cyl:", l_cyl, "[m]", "l_fus:", l_fus, "[m]", "S_wfus:", S_wfus, "[m^2]"))
    #^^chatgpt goated for this
    return S_wfus, l_fus, l_cabin, l_nc,w_fus,w






