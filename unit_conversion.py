#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 14:00:36 2024

@author: shiyu
"""

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

print(convert_units(4,'m^2', False))  # Converts 10 feet to meters (default to_metric=True)

