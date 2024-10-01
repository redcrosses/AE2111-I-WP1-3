import numpy as np
from consts import *
import math
import matplotlib.pyplot as plt
from intersect import intersection
import inspect

wing_loading = np.arange(0.1,9100,100) #<- 0.1 avoids the division by zero warning

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

def planform_print(span, root_c, tip_c,sweep_quart):
	plt.clf()
	x = [0,0,span/2, span/2,0]
	y = [root_c, 0, 0.25*root_c + np.tan(sweep_quart)*span - 0.25*tip_c, 0.25*root_c + np.tan(sweep_quart)*span + 0.75*tip_c,root_c]
	plt.plot(x,y, 'ro-')
	plt.gca().set_aspect('equal', 'box')
	plt.show()

def mainloop(aspect_ratio): #I have it set to vary aspect ratio. you can change it to something else
	x_const = [100*i for i in range(0,91)]
	lines = [([min_speed_list(clmax_landing)]*91, x_const), ([field_length_list(clmax_landing)]*91, x_const),(wing_loading, cruise_speed_list())]
	labels = ["Minimum speed","Landing Field Length","Cruise speed","Climb Gradient 1","Climb Gradient 2","Climb Gradient 3","Climb Gradient 4","Climb Gradient 5","Takeoff Field Length"]
	
	gradient = climb_gradient(0)
	lines.append((wing_loading, 0.5*gradient.climb_gradient_list()))
	for i in range(1,5):
		gradient = climb_gradient(i)
		lines.append((wing_loading, gradient.climb_gradient_list()))
	lines.append((wing_loading, to_field_length_list(aspect_ratio)))
	
	design_point = find_design_point(0,lines) #make minimum speed line the target line to intersect with

	# for i in range(len(lines)): #plotting all lines
	# 	plt.plot(lines[i][0], lines[i][1], label = labels[i])
	# plt.plot(design_point[0], design_point[1],"ro")
	# plt.ylim(0,1)
	# plt.legend()
	# plt.grid()
	# plt.show() #uncomment to show the matching diagram

	##HLDs and Control surfaces 
	S =  max_to_mass*9.81 / float(design_point[0][0]) #<- column and row position avoids deprecation warning

	#wing parameters calculated from wing area
	span = np.sqrt(aspect_ratio*S)
	chord_root = 2*S / ((1+ taper_ratio)*span)
	chord_tip = chord_root *  taper_ratio
	y_1 = 0.15*span/2 #position of the beginning of the HLD; 15% of the half-span
	print("\nS: %.5f [m^2] \nSpan: %.5f [m] \nRoot Chord: %.5f [m] \nTip Chord: %.5f [m] \ny_1: %.5f [m] \nHLD margin: %.5f" % (S, span, chord_root,  chord_tip, y_1,  hld_margin))

	delta_CLmax =  clmax_landing -  CLmax_wingclean

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
	c_Drag =  C_D0_cruise + np.power( c_L_cruise,2)/(np.pi*aspect_ratio*cruise_oswald_efficiency)
	Drag = c_Drag * 0.5*cruise_density*np.power(cruise_velocity,2) * S

	SAR =  specific_fuel_energy* jet_efficiency/(Drag)
	frame = inspect.currentframe()
	name,_,_,argvalue = inspect.getargvalues(frame)
	iteratedvalue = argvalue[name[0]]
	print(name, iteratedvalue)
	return [SAR, iteratedvalue, S, span, chord_root, chord_tip, S_wf, y_1, y_2, b_2[1]]

results = []
SARs = []
for iterator in range(10,160,10): #iterating the design
	var = iterator/10
	run = mainloop(var)
	SARs.append(run[0])
	results.append(run[1:])
# print(SARs)

optimalSAR = max(SARs) #optimal is found when the SAR is maximum in the iterated range
optimal = results[SARs.index(optimalSAR)] #print the optimal results based on optimal aspect ratio
optimal.insert(0,optimalSAR)
labels = [["Optimal SAR:",'Iterated value:', 'S:', 'Span:', 'Chord_root:', 'Chord_tip:', 'S_wf:','y_1 (HLD):', 'y_2 (HLD):', 'b_2 (Aileron):'],["[m/kg]", "", "","[m^2]","[m]","[m]","[m^2]","","","",""]]
print("\n###{:^36}###".format("RESULTS"))

print("Iterated variable: {:>18}".format(mainloop.__code__.co_varnames[0]))
for i in range(len(optimal)):
	print("{:24} {:.5f} {:16}".format(labels[0][i],optimal[i],labels[1][i]))

planform_print(optimal[3],optimal[4],optimal[5], sweep_quarter)

#to iterate a different parameter, try to change the mainloop function argument to the desired one. Then, in the main loop, change the iterated range to the desired one together with the divisor for the var.