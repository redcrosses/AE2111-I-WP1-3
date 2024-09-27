import numpy as np
inputs = 91

max_to_mass = 265625 #kg

approach_speed = 78
landing_massfraction = 0.79
landing_temp_diff = 15
landing_fieldlengthreq = 1981.2
cruise_massfraction = 0.95
cruise_altitude = 9449.8
cruise_minmach = 0.85

aspect_ratio = 9

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
ref_altitude = 100
cruise_temp = 288.15-0.0065*(cruise_altitude-ref_altitude)
parasite_drag = 0.0075
span_efficiency = 0.97

climbgradient_massfraction = [1, 1, 1, 1, 0.79]
climbgradient_zerodrag = [0.0758, 0.0563, 0.0363, 0.0168, 0.0558] # taken from drag polar D:38-
climbgradient_gradient = [3.2, 0 , 0, 1.2, 2.1]
climbgradient_oswaldfactor = [0.867622569, 0.828622569, 0.828622569, 0.789622569, 0.867622569]

## HLDs and Control Surfaces
sweep_sixc = 0.445791998
sweep_quarter = np.radians(31.35)
taper_ratio = 0.2*(2-sweep_quarter)
hld_margin = 1

CLratio = 0.96
Clmax = 2.097
clmax_landing = 2.6 #CL design
CLmax_wingclean = 2.013

c_ratio = 1.21 #Fowler flap ratio
delta_clmax = 1.3 * c_ratio
C_lalpha = 6.7614
C_d0 = 0.0008

P = np.radians(20)
stall_speed = 69.44
dalpha = 0.4581

tau = 0.4

TSFC = 1

