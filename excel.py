import numpy as np
import consts as c
import math
import matplotlib.pyplot as plt



wing_loading = np.arange(0,9100,100)

def min_speed_list():
	landing_airdensity = 101325/(287*(c.landing_temp_diff+288.15))
	return c.clmax_landing*((c.approach_speed/1.23)**2)*landing_airdensity/(2*c.landing_massfraction)


def field_length_list():
	landing_airdensity = 101325/(287*(c.landing_temp_diff+288.15))
	return 1/c.landing_massfraction*c.landing_fieldlengthreq/0.45*landing_airdensity*c.clmax_landing/2


def cruise_speed_list():
	cruise_pressure = 101325*(1+(-0.0065*c.cruise_altitude/288.15))**(-9.81/(-0.0065*287))
	cruise_totalpressure = cruise_pressure*(1+(1.4-1)/2*c.cruise_minmach**2)**3.5
	cruise_deltapressure = cruise_totalpressure/101325
	cruise_thrustlapse = cruise_deltapressure*(1-(0.43+0.014*c.bypass_ratio)*math.sqrt(c.cruise_minmach))
	zerolift_drag = c.friction_coefficient*c.wetted_ratio
	oswald_efficiency = 1/(math.pi*c.bypass_ratio*c.parasite_drag+(1/c.span_efficiency))
	cruise_density = cruise_pressure/(287*c.cruise_temp)
	cruise_speed = math.sqrt(1.4*287*c.cruise_temp)*c.cruise_minmach
	return c.cruise_massfraction/cruise_thrustlapse*((zerolift_drag*0.5*cruise_density*cruise_speed**2)/(c.cruise_massfraction*wing_loading)+(c.cruise_massfraction*wing_loading)/(math.pi*c.bypass_ratio*oswald_efficiency*0.5*cruise_speed**2*cruise_density))



class climb_gradient():
	def __init__(self, num):
		self.num=num

	def climb_gradient_list(self):
		climbgradient_maxcl = math.sqrt(c.climbgradient_zerodrag[self.num]*math.pi*c.bypass_ratio*c.climbgradient_oswaldfactor[self.num])
		climbgradient_best_speed = np.sqrt((wing_loading*2)/(1.225*climbgradient_maxcl))
		climbgradient_mach_number = climbgradient_best_speed/np.sqrt(1.4*287*288.15)
		climbgradient_total_pressure = (1+(1.4-1)/2*np.power(climbgradient_mach_number,2))*101325
		climbgradient_delta_pressure = climbgradient_total_pressure/101325
		climbgradient_thrust_lapse = climbgradient_delta_pressure*(1-(0.43+0.014*c.bypass_ratio)*np.sqrt(climbgradient_mach_number))
		oswald_efficiency = 1/(math.pi*c.bypass_ratio*c.parasite_drag+(1/c.span_efficiency))
		return 2*c.climbgradient_massfraction[self.num]/climbgradient_thrust_lapse*(c.climbgradient_gradient[self.num]/100+2*math.sqrt(c.climbgradient_zerodrag[self.num]/(math.pi*c.climbgradient_oswaldfactor[self.num]*c.bypass_ratio)))


def to_field_length_list():
	to_pressure = 101325*(1+(500*-0.0065)/(288.15))**(-9.81/(-0.0065*287))
	to_density = to_pressure/(287*c.to_temp)
	to_velocity = np.sqrt((9000*2)/(to_density*c.to_cL))
	to_mach_number = to_velocity/np.sqrt(1.4*287*c.to_temp)
	to_total_pressure = (1+(1.4-1)/2*np.power(to_mach_number,2))*c.to_pressure
	to_delta_pressure = to_total_pressure/101325
	to_thrust_lapse = to_delta_pressure*(1-(0.43+0.014*c.bypass_ratio)*np.sqrt(to_mach_number))
	return (1.15*to_thrust_lapse*np.sqrt(2 * wing_loading/(c.to_field_length * 0.85 * to_density * 9.81 * math.pi * c.to_oswald_efficiency * c.aspect_ratio))+ 2*(4*11)/c.to_field_length)


x_const = [100*i for i in range(0,91)]


plt.plot(wing_loading, cruise_speed_list(), label = "Cruise speed")
plt.plot([min_speed_list()]*91, x_const, label = "Minimum speed")
plt.plot([field_length_list()]*91, x_const, label = "Landing Field Length")
gradient = climb_gradient(0)
plt.plot(wing_loading, 0.5*gradient.climb_gradient_list(), label = "Climb Gradient 1")
for i in range(1,5):
	gradient = climb_gradient(i)
	plt.plot(wing_loading, gradient.climb_gradient_list(), label = "Climb Gradient " + str(i+1))
plt.plot(wing_loading, to_field_length_list(), label = "Takeoff Field Length")

plt.ylim(0,1)
plt.legend()
plt.grid()
plt.show()

