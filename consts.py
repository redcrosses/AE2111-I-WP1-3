import numpy as np
inputs = 91

max_to_mass = 265625 #kg
aspect_ratio = 8.02

approach_speed = 78 
landing_massfraction = 0.79
landing_temp_diff = 15
landing_fieldlengthreq = 1981.2
cruise_massfraction = 0.95
cruise_altitude = 9449.8
ref_altitude = 100
cruise_temp = 288.15-0.0065*(cruise_altitude-ref_altitude)
cruise_minmach = 0.85
cruise_speed = 0.85 * np.sqrt(1.4*287*cruise_temp)

to_massfraction = 0.85
to_altitude = 7400
to_cL = 2.036181377 #to be changed
to_field_length = 3048
to_oswald_efficiency = 0.8286
to_temp = 298
to_pressure = 95457.84253

bypass_ratio = 10
wetted_ratio = 6
friction_coefficient = 0.0028
parasite_drag = 0.0075
span_efficiency = 0.97

climbgradient_massfraction = [1, 1, 1, 1, 0.79]
climbgradient_zerodrag = [0.0758, 0.0563, 0.0363, 0.0168, 0.0558] # taken from drag polar D:38-
climbgradient_gradient = [3.2, 0 , 0, 1.2, 2.1]
climbgradient_oswaldfactor = [0.867622569, 0.828622569, 0.828622569, 0.789622569, 0.867622569]

## HLDs and Control Surfaces
sweep_quarter = np.radians(28.39) #actually leading edge lmao
taper_ratio = 0.2*(2-sweep_quarter)
hld_margin = 1

Clratio =  0.8  #1.04
Clmax = 1.797
CLmax_wingclean = Clratio * Clmax
clmax_landing = 2.6 #C_L design

delta_c_to_cf = 0.64
c_ratio_TE = 1+ 0.35 * delta_c_to_cf #single-slotted fowler flap trailing edge
delta_clmax = 1.3 * c_ratio_TE #single-slotted fowler 
C_lalpha = 6.7614
c_ratio_LE = 1.1 #slat at the leading edge
cl_leadingedge = 0.4*c_ratio_LE
C_d0 = 0.0008 #wing 

P = np.radians(20)
stall_speed = 69.44
dalpha = 0.4581

tau = 0.4

specific_fuel_energy = 44e6
jet_efficiency = 0.75

# tsfc = 22*np.power(bypass_ratio,-0.19)*1E-6 #from adsee estimation for turbofans?? or how else do you do it lol
# print(tsfc)
c_L_cruise = 0.52924 #from the cruise angle of attack section
C_D0_cruise = 0.0168 #from Fred's excel Drag polar section

