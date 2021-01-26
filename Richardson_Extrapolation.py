# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:27:04 2020

@author: Matthijs Schrage
"""

import numpy as np

def truncation_error(Q_h, Q_2h, Q_4h):
	
	two_power_of_p = np.abs( (Q_2h - Q_4h ) / ( Q_h - Q_2h ) )
	
	print("2^p = " + str(two_power_of_p))
	
	p = np.log2(two_power_of_p)
	
	print("p = " + str(p))
	
	max_error = (1 / ( (2 ** round(p)) - 1) ) * (Q_h - Q_2h)
	
	print("Error: " + str(max_error))
	
	return round(p)

truncation_error(1.85193712, 1.85193809, 1.85195370)




