import numpy as np
from consts import *
from fuselage import *
import math
import matplotlib.pyplot as plt
from intersect import intersection
import inspect

###FUNCTIONS AND CLASSES###
wing_loading = np.arange(0.1,9100,100) #<- 0.1 avoids the division by zero warning
plt.figure(figsize=(15,5))
def Class_1_est(Liftoverdrag,h_CR,V_CR,jet_eff,energy_fuel,R_nom, R_div,t_E, f_con, m_OE, M_pl):
    # energy fuel is like the weird 41sth/kw or idk
	#t_E is the Loiter time in emergencies
    g=9.81
    R_lost = 1/0.7 * (Liftoverdrag) * (h_CR + (V_CR**2)/(2*g)) /1000 # lost range via : LIft/Drag , height of cruise, velocity cruise
    R_eq = (R_nom + R_lost)*(1+f_con)+1.2 
    R_div + t_E * V_CR # nominal and lost range plus fraction trip fuel for contingency
    m_f = 1- np.exp((-R_eq * g * 1000)/(jet_eff * energy_fuel * liftoverdrag))
    M_MTO = M_pl /(1-(m_OE)-(m_f))  # m_OE taken from reference or hallucinated
    M_f = m_f * M_MTO
    M_OE= m_OE * M_MTO
    print(R_lost, R_eq, R_div, m_f)
    print("Operating empty: {:f}, \nFuel: {:f}, \nMax TO {:f}, \nFuel mass fraction: {:f}".format(M_OE, M_f, M_MTO, m_f))
    return(M_OE, M_f, M_MTO, m_f) #in kilos small m is mass fraction, big M is acutal mass

def min_speed_list(clmax_landing):
	landing_airdensity = 101325/(287*( landing_temp_diff+288.15))
	return  clmax_landing*(( approach_speed/1.23)**2)*landing_airdensity/(2* landing_massfraction)


def field_length_list(clmax_landing):
	landing_airdensity = 101325/(287*( landing_temp_diff+288.15))
	return 1/ landing_massfraction* landing_fieldlengthreq/0.45*landing_airdensity* clmax_landing/2


def cruise_speed_list():
	cruise_pressure = 101325*(1+(-0.0065* cruise_altitude/288.15))**(-9.81/(-0.0065*287))
	cruise_totalpressure = cruise_pressure*(1+(1.4-1)/2* cruise_minmach**2)**3.5
	cruise_deltapressure = cruise_totalpressure/101325
	cruise_thrustlapse = cruise_deltapressure*(1-(0.43+0.014* bypass_ratio)*math.sqrt( cruise_minmach))
	zerolift_drag =  friction_coefficient* wetted_ratio
	oswald_efficiency = 1/(math.pi* aspect_ratio* parasite_drag+(1/ span_efficiency))
	cruise_density = cruise_pressure/(287* cruise_temp)
	cruise_speed = math.sqrt(1.4*287* cruise_temp)* cruise_minmach
	return  cruise_massfraction/cruise_thrustlapse*((zerolift_drag*0.5*cruise_density*cruise_speed**2)/( cruise_massfraction*wing_loading)+( cruise_massfraction*wing_loading)/(math.pi* aspect_ratio*oswald_efficiency*0.5*cruise_speed**2*cruise_density))

class climb_gradient():
	def __init__(self, num):
		self.num=num

	def climb_gradient_list(self):
		climbgradient_maxcl = math.sqrt( climbgradient_zerodrag[self.num]*math.pi* bypass_ratio* climbgradient_oswaldfactor[self.num])
		climbgradient_best_speed = np.sqrt((wing_loading*2)/(1.225*climbgradient_maxcl))
		climbgradient_mach_number = climbgradient_best_speed/np.sqrt(1.4*287*288.15)
		climbgradient_total_pressure = (1+(1.4-1)/2*np.power(climbgradient_mach_number,2))*101325
		climbgradient_delta_pressure = climbgradient_total_pressure/101325
		climbgradient_thrust_lapse = climbgradient_delta_pressure*(1-(0.43+0.014* bypass_ratio)*np.sqrt(climbgradient_mach_number))
		oswald_efficiency = 1/(math.pi* aspect_ratio* parasite_drag+(1/ span_efficiency))
		return 2* climbgradient_massfraction[self.num]/climbgradient_thrust_lapse*( climbgradient_gradient[self.num]/100+2*math.sqrt( climbgradient_zerodrag[self.num]/(math.pi* climbgradient_oswaldfactor[self.num]* bypass_ratio)))

def to_field_length_list(aspect_ratio):
	to_pressure = 101325*(1+(500*-0.0065)/(288.15))**(-9.81/(-0.0065*287))
	to_density = to_pressure/(287* to_temp)
	to_velocity = np.sqrt((9000*2)/(to_density* to_cL))
	to_mach_number = to_velocity/np.sqrt(1.4*287* to_temp)
	to_total_pressure = (1+(1.4-1)/2*np.power(to_mach_number,2))* to_pressure
	to_delta_pressure = to_total_pressure/101325
	to_thrust_lapse = to_delta_pressure*(1-(0.43+0.014* bypass_ratio)*np.sqrt(to_mach_number))
	return (1.15*to_thrust_lapse*np.sqrt(2 * wing_loading/( to_field_length * 0.85 * to_density * 9.81 * math.pi *  to_oswald_efficiency * aspect_ratio))+ 2*(4*11)/ to_field_length)

def find_design_point(target_pos, lines_arr):
	target = lines_arr[target_pos]
	intersections = []
	yvals = []
	for line in lines_arr[2:]: #first two lines are removed from the intersection; first one is itself and second one is parallel
		intrsctn = intersection(target[0],target[1],line[0],line[1])
		intersections.append(intrsctn)
		yvals.append(intrsctn[1])
	return intersections[yvals.index(max(yvals))] #return the point in the intersections list that has the maximum y-value of the intersection.

def find_cg(fuselage_length, nose_cone_length, cabin_length, m_fuel):
	#LEMAC calculation
	xc_oew = 0.25
	M_empennage = 0.017
	M_fuselage = 0.101
	M_equipment = 0.089
	M_wing = 0.122
	M_Nacelle = 0.0056
	M_Prop = 0.0225
	fuselage_group = np.array([[M_empennage, M_fuselage, M_equipment],[0.9*fuselage_length, 0.4*fuselage_length, 0.4*fuselage_length]])
	wing_group = np.array([[M_wing, M_Nacelle, M_Prop],[0.4*MAC, -3, -3]])
	fus_sum = fuselage_group.prod(axis=0).sum()
	wing_sum = wing_group.prod(axis=0).sum()
	M_fus_sum = fuselage_group[0].sum()
	M_wing_sum = wing_group[0].sum()
	fus_pos = fus_sum/M_fus_sum
	wing_pos = wing_sum/M_wing_sum
	print(fus_pos, wing_pos)
	X_LEMAC = fus_pos + MAC*(wing_pos/MAC * M_wing_sum/M_fus_sum - xc_oew*(1+M_wing_sum/M_fus_sum))
	X_TEMAC = X_LEMAC + MAC

	#CG location
	m_Payload = 18960/max_to_mass
	cg_matrix = np.array([[m_OE, m_Payload, m_fuel],[X_LEMAC + xc_oew*MAC, nose_cone_length + 0.5*cabin_length, X_LEMAC+0.4*MAC]])
	moments = cg_matrix.prod(axis=0)
	print(moments)
	cg_positions = np.array([[cg_matrix[1][0], cg_matrix[0][0]],
				 [(moments[0]+moments[1])/(cg_matrix[0][0]+cg_matrix[0][1]),(cg_matrix[0][0]+cg_matrix[0][1])],
				 [(moments[0]+moments[1]+moments[2])/(cg_matrix[0][0]+cg_matrix[0][1]+cg_matrix[0][2]),(cg_matrix[0][0]+cg_matrix[0][1]+cg_matrix[0][2])],
				 [(moments[0]+moments[2])/(cg_matrix[0][0]+cg_matrix[0][2]),(cg_matrix[0][0]+cg_matrix[0][2])]]) #OEW, WOE+WP, WOE+WP+WF, WOE+WF
	# print(cg_positions)
	#plotting the cgs
	a,b = zip(*np.vstack([cg_positions,[cg_matrix[1][0], cg_matrix[0][0]]])) #added the starting point for proper plotting
	plt.subplot(2,3,3)
	plt.title("Cg Positions")
	plt.plot(a,b, 'bo-')
	plt.plot((X_LEMAC, X_TEMAC),(1,1),'ro-')
	plt.ylim(0,1.1)
	print("\n\033[1m\033[4m Cg Locations & Mass Fractions \033[0m")
	print("OEW: {0} \nWOE+WP: {1} \nWOE+WP+WF: {2} \nWOE+WF: {3} ".format(cg_positions[0],cg_positions[1],cg_positions[2],cg_positions[3] ))
	return cg_positions

def empennage_size(l_fus, cg_aft, l_MAC, S_wing, b):
	htail_aero_centre_location = l_fus - 4
	htail_moment_arm_cg_aft = htail_aero_centre_location - cg_aft 
	htail_c_v = 0.95
	htail_area = (htail_c_v * l_MAC * S_wing) / (htail_moment_arm_cg_aft)

	vtail_aero_centre_location = htail_aero_centre_location - 2
	vtail_moment_arm_cg_aft = vtail_aero_centre_location - cg_aft
	vtail_c_v = 0.066
	vtail_area = (vtail_c_v * b * S_wing) / (vtail_moment_arm_cg_aft)
	print("\n\033[1m\033[4m Empennage \033[0m")
	print("Hor.Tail Location: %.5f [m] \nHor.Tail Area: %.5f [m^2] \nVer.Tail Location: %.5f [m] \nVer.Tail Area: %.5f [m^2]" % (htail_aero_centre_location, htail_area, vtail_aero_centre_location, vtail_area))
	return htail_aero_centre_location, htail_area, vtail_aero_centre_location, vtail_area

def planform_print(span, root_c, tip_c,sweep_quart):
	plt.subplot(2,3,4)
	plt.title("Wing Planform")
	x = [0,0,span, span,0]
	y = [root_c, 0, 0.25*root_c + np.tan(sweep_quart)*span - 0.25*tip_c, 0.25*root_c + np.tan(sweep_quart)*span + 0.75*tip_c,root_c]
	plt.plot(x,y, 'ro-')
	plt.gca().set_aspect('equal', 'box')

def matchingdiag_print(lines, labels, design_point):
	plt.subplot(2,3,1)
	plt.title("Matching Diagram")
	for i in range(len(lines)): #plotting all lines
		plt.plot(lines[i][0], lines[i][1], label = labels[i])
	plt.plot(design_point[0], design_point[1],"ro")
	plt.gca().set_aspect('auto','box')
	plt.ylim(0,1)
	plt.legend()
	plt.grid()


def optimisation(clmax_landing, max_to_mass):
    #matching diagram
	x_const = [100*i for i in range(0,91)]
	global lines, labels, design_point
	lines = [([min_speed_list(clmax_landing)]*91, x_const), ([field_length_list(clmax_landing)]*91, x_const),(wing_loading, cruise_speed_list())]
	labels = ["Minimum speed","Landing Field Length","Cruise speed","Climb Gradient 1","Climb Gradient 2","Climb Gradient 3","Climb Gradient 4","Climb Gradient 5","Takeoff Field Length"]
	
	gradient = climb_gradient(0)
	lines.append((wing_loading, 0.5*gradient.climb_gradient_list()))
	for i in range(1,5):
		gradient = climb_gradient(i)
		lines.append((wing_loading, gradient.climb_gradient_list()))
	lines.append((wing_loading, to_field_length_list(aspect_ratio)))
	
	design_point = find_design_point(0,lines) #make minimum speed line the target line to intersect with
	thrust_max = float(design_point[1][0])*max_to_mass*9.81 /1000

	##HLDs and Control surfaces 
	S =  max_to_mass*9.81 / float(design_point[0][0]) #<- column and row position avoids deprecation warning

	#wing parameters calculated from wing area
	span = np.sqrt(aspect_ratio*S)
	chord_root = 2*S / ((1+ taper_ratio)*span)
	chord_tip = chord_root *  taper_ratio
	y_1 = 0.10*span/2 #position of the beginning of the HLD; 15% of the half-span
	print("\nS: %.5f [m^2] \nSpan: %.5f [m] \nRoot Chord: %.5f [m] \nTip Chord: %.5f [m] \ny_1: %.5f [m] \nHLD margin: %.5f" % (S, span, chord_root,  chord_tip, y_1,  hld_margin))

	#sweep angle relations (probably should be a function lol)
	global sweep_LE, sweep_sixc, sweep_half
	sweep_LE = np.tan(sweep_quarter) + 0.25 * (2*chord_root)/(span)*(1-taper_ratio)
	sweep_sixc = np.tan(sweep_LE) - 0.6 * (2*chord_root)/(span)*(1-taper_ratio)
	sweep_half = np.tan(sweep_LE) - 0.5 * (2*chord_root)/(span)*(1-taper_ratio)
	#will need out of loop and this is easier
	global t_r
	t_r = chord_root*t_cratio
 
	#HLD and Control surfaces placement
	delta_CLmax =  clmax_landing - cl_leadingedge -  CLmax_wingclean

	S_ratio = delta_CLmax/(0.9 *  delta_clmax * np.cos(sweep_sixc))

	S_wf = S_ratio * S
	y_2 = np.min(np.roots([-(chord_root - chord_tip)/(span), chord_root, ((chord_root - chord_tip)/(span) * np.power(y_1,2) - chord_root*y_1 - S_wf/2)])) +  hld_margin
	print("y_2 for HLD:", y_2)
	inter = 2* (chord_root - chord_tip)/span
	C_lp = -4 * ( C_lalpha +  C_d0)/(S * np.power(span,2)) * ((chord_root/3 * np.power(span/2,3) - inter * np.power(span/2,4)/4))

	C_lda = -( P * C_lp)/( dalpha) * (span/(2* stall_speed))
	b_2 = np.roots([2/3 * (chord_root - chord_tip)/span, -chord_root/2, 0, chord_root/2 * np.power(y_2, 2) - 2/3 * (chord_root - chord_tip)/span * np.power(y_2, 3) + (C_lda*S*span)/(2* C_lalpha* tau)])
	print("b_2 for alieron (select the reasonable one):",b_2)

	cruise_density = (101325*(1+(-0.0065* cruise_altitude/288.15))**(-9.81/(-0.0065*287))) /(287* cruise_temp)
	cruise_velocity =  cruise_minmach*np.sqrt(1.4* cruise_temp*287)
	cruise_oswald_efficiency = 4.61*(1-0.045*np.power(aspect_ratio,0.68))*np.power(math.cos( sweep_quarter),0.15) - 3.1
	c_Drag =  c_d0initial + np.power(c_L_cruise,2)/(np.pi*aspect_ratio*cruise_oswald_efficiency)
	Drag = c_Drag * 0.5*cruise_density*np.power(cruise_velocity,2) * S

	SAR =  specific_fuel_energy* efficiency_tf/(Drag)
	print("SAR:",SAR)
	frame = inspect.currentframe()
	name,_,_,argvalue = inspect.getargvalues(frame)
	iteratedvalue = argvalue[name[0]]
	print(name, iteratedvalue)
	diff = span/2 - b_2[1]
	print("Diff:", diff)
	return [SAR, iteratedvalue, S, span, chord_root, chord_tip, S_wf, y_1, y_2, b_2[1], thrust_max, diff]

######################################################
#class 1 weight estimation
M_OE, M_f, M_MTO, m_f = Class_1_est(liftoverdrag, cruise_altitude, cruise_speed, jet_eff, specific_fuel_energy, R_nominal, R_diversion, t_E, f_con, m_OE, M_pl)

#iterating the design (matching diagram, wing sizing, hld and control surfaces)
results = []
diffs = []
current = 2500
iterator = 10
previous = 0
while True: #iterating the design
	var = current/1000
	run = optimisation(var, M_MTO) #<-- optimisation function call
	diffs.append(run[-1])
	results.append(run[:-1])
	current += iterator
 
	if(run[-1]<0): #this optimization ensures the hld and alierons fill the whole wing. Optimal SAR is only optimal if aspect ratio is optimal!
		optimal = previous
		matchingdiag_print(lines, labels, design_point)
		break
	else:
		previous = run
# print(SARs)
#optimal result
# optimaldiff = min(diffs) #optimal is found when the SAR is maximum in the iterated range
# optimal = results[diffs.index(optimaldiff)] #print the optimal results based on optimal aspect ratio

labels = [["Optimal SAR:",'Iterated value:', 'S:', 'Span:', 'Chord_root:', 'Chord_tip:', 'S_wf:','y_1 (HLD):', 'y_2 (HLD):', 'b_2 (Aileron):', 'Maximum Thrust:', 'Diff:'],["[m/kg]", "","[m^2]", "[m]","[m]","[m]","[m^2]","[m]","[m]","[m]",'[kN]','[m]']]
print("\n\033[1m\033[4m Optimal Results [m] \033[0m")
print("Iterated variable: {:>18}".format(optimisation.__code__.co_varnames[0]))
for i in range(len(optimal)):
	print("{:24} {:.5f} {:16}".format(labels[0][i],optimal[i],labels[1][i]))
planform_print(optimal[3]/2,optimal[4],optimal[5], sweep_quarter)

#mean aerodynamic chord
MAC = 2/3 * optimal[4] * ((1+taper_ratio+taper_ratio**2)/(1+taper_ratio))
#finding the fuselage dimensions
S_wfuselage, l_fuselage, l_cabin, l_ncone = fuselage(83.1) #output: fuselage wetted surface area, fuselage length, cabin length, nose cone length
cg_positions = find_cg(float(l_fuselage), l_ncone, l_cabin,m_f)
cg_aft = np.max(cg_positions[:,0])
x_htail, htail_area, x_vtail, vtail_area = empennage_size(l_fuselage, cg_aft, MAC,optimal[2],optimal[3])

#new drag estimation (fast estimation)
S_wwing = 1.07*2*optimal[2]
S_wHT = 1.05 * 2 * htail_area
S_wVT = 1.05 * 2 * vtail_area
S_wnacelles = 1 #todo
c_d0new = 1.15 * (1/optimal[2] * (S_wfuselage * 0.08 + S_wwing * 0.007 + S_wnacelles*0.06 + S_wHT * 0.008 + S_wVT * 0.008))
print("\nOLD c_d0: {0}, NEW c_d0: {1}".format(c_d0initial, c_d0new)) #cd0 is wrong; 10x too large

#Weight estimation
n_max = 2.5 #max loading factor estimmation
n_ult = 1.5*n_max
b_s = optimal[3]/np.cos(sweep_half)
ZFW = (M_MTO - M_f)*9.81
M_wing = (6.67e-3 * np.power(b_s,0.75)*(1+np.sqrt(1.905/b_s)*np.power(n_ult,0.55)*np.power((b_s/t_r)/(ZFW/optimal[2]),0.30)))*ZFW/9.81 #CHANGE TO A CONSISTENT METHOD FROM BOOK
M_fuselage = 1 #todo
M_powerplant = 1 #todo
M_empennage = 1 #todo
plt.tight_layout()
plt.show() #uncomment to show dashboard