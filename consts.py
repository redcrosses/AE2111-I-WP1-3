inputs = 91
clmax_landing = 2.1
approach_speed = 78
landing_massfraction = 0.79
landing_temp_diff = 15
landing_fieldlengthreq = 1786
cruise_massfraction = 0.95
cruise_altitude = 9449.8
cruise_minmach = 0.85

aspect_ratio = 10

to_massfraction = 0.85
to_altitude = 7400
to_cL = 0.633 #to be changed
to_field_length = 3048
to_oswald_efficiency = 0.8286
to_temp = 288.15-0.0065*(to_altitude)

bypass_ratio = 10
wetted_ratio = 6
friction_coefficient = 0.0027
ref_altitude = 100
cruise_temp = 288.15-0.0065*(cruise_altitude-ref_altitude)
parasite_drag = 0.0075
span_efficiency = 0.97

climbgradient_massfraction = [1, 1, 1, 1, 0.79]
climbgradient_zerodrag = [0.0162, 0.0362, 0.0357, 0.0557, 0.0552, 0.0752] # taken from drag polar D:38-
climbgradient_gradient = [3.2, 0 , 2.4, 1.2, 2.1]
climbgradient_oswaldfactor = [0.867622569, 0.828622569, 0.828622569, 0.789622569, 0.867622569]