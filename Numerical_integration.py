# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:00:08 2020

@author: Matthijs Schrage
"""

import numpy as np
import sympy as sym

#Some function to help us with syntax
def sin(x):
	return sym.sin(x)

def cos(x):
	return sym.cos(x)

def sqrt(x):
	return sym.sqrt(x)

x = sym.Symbol("x")

points = [0, 1/5, 1/3, 3/2]

def rectangle(sequence, bounds, method="composite_left"):
	
	partitions = len(sequence) - 1
	area = 0
	
	if method == "composite_left":
	
		for n in range(partitions):
			
			area += sequence[n]
	
	elif method == "composite_right":
			
		for n in range(partitions):
			
			area += sequence[n+1]
		
	area *= ( (bounds[1] - bounds[0]) / partitions )
	
	print("Area " + method + " rectangle integration: " + str(area))
	
	return area

def trapezoidal(sequence, bounds):
	
	partitions = len(sequence) - 1
	area = 0
	
	for i in range(partitions):
			
		area += ( (bounds[1] - bounds[0]) / partitions ) * ( sequence[i] + sequence[i+1] ) /2
	
	print("Area trapezoidal integration: " + str(area))

def midpoint(function, stepsize, bounds):
	
	partitions = int(np.ptp(bounds) / stepsize)
	area = 0
	
	print(partitions)
	
	for i in range(partitions):
		
		print((bounds[0] + (i)*stepsize))
		print((bounds[0] + (i+1)*stepsize))
		print("_____________________")
		
		area += ( (bounds[0] + (i+1)*stepsize) - (bounds[0] + (i)*stepsize) ) * function(((bounds[0] + (i+1)*stepsize) + (bounds[0] + (i)*stepsize))/2)
		
	print("Area midpoint integration: " + str(area))

rectangle(points, [0,1], method="composite_left")
rectangle(points, [0,1], method="composite_right")
trapezoidal(points, [0,1])

# f = sym.lambdify(x, 1-cos(x), 'numpy')
# midpoint(f, np.pi, [0,np.pi])




