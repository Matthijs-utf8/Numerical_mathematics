# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 16:29:45 2020

@author: Matthijs Schrage
"""

import numpy as np
import sympy as sym

x = sym.Symbol("x")

#Some function to help us with syntax
def sin(x):
	return sym.sin(x)

def cos(x):
	return sym.cos(x)

def sqrt(x):
	return sym.sqrt(x)

def linear(f, coordinates, target):
	
	#Assert that we only interpolate bwteen two points
	assert len(coordinates) == 2
	
	#This is the formula (hard-coded) for linear interpolation
	L_x = ( (target - coordinates[1]) / (coordinates[0] - coordinates[1]) ) * f(coordinates[0]) + ( (target - coordinates[0]) / (coordinates[1] - coordinates[0]) ) * f(coordinates[1]) 
	
	print("Linear interpolation in " + str(target) + ": " + str(L_x) )

#This function is used in the langrangian interpolation to help make the code more compact and flexible
def helper(i):
	
	# If len(coordinates) == 3:
	i = i % 3
	
	return i

def lagrangian(f, coordinates, target):
	
	coordinates = sorted(coordinates)
	
	L_n = 0
	
	for i in range(len(coordinates)):
		
		#For every coordinate, we calculate the lagrangian interpolation
		L_i_n = ( ( (x - coordinates[helper(i+1)]) * (x - coordinates[helper(i+2)] ) ) / ( (coordinates[i] - coordinates[helper(i+1)]) * (coordinates[i] - coordinates[helper(i+2)]) ) ) * f(coordinates[i])
		L_n += L_i_n
	
	#It is not stricly necessary to make a function, but it is just more neat
	L_n = sym.lambdify(x, L_n)
	
	print("Lagrangian interpolation in " + str(target) + ": " + str(L_n(target) ) )

fx2 = sym.lambdify(x, x*(-1 + x**2))
# linear(fx2, [2,3], 2.5)
lagrangian(fx2, [-1, 0, 1], -0.5)

def lagr(coordinates, f_values, k, n):
	
	assert len(coordinates) == len(f_values)
	assert k in coordinates
	
	coordinates.remove(k)
	
	L_k_n = 1
	
	for i in coordinates:
		
		L_k_n *= (x - i) / (k - i)
	
	print(L_k_n)

# lagr([0,1,2,3], [0.5*np.sqrt(2), np.e, np.pi, 2*np.sqrt(3)], 2, 3)



