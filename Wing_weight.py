#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 14:44:08 2024

@author: shiyu
"""
import numpy as np

def Weight_wing(S_w,W_fw,Aspect_ratio,sweep,q,taper,thick_chord):
    return 0.036*(S_w**0.758)*(W_fw**0.0035)*(((Aspect_ratio/(np.cos(sweep)**2)))**0.6)*(q**0.006)*(taper**0.04)*(((100*thick_chord)/np.cos(sweep))**-0.3)
    