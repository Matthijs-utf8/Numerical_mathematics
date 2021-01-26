# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:54:53 2020

@author: Matthijs Schrage
"""

import numpy as np
import sympy as sym

x = sym.Symbol("x")

def simple_difference_formula(f, point, step, method="central"):
	
	if method == "forward":
		
		f_prime = ( f(point + step) - f(point) ) / step
	
	elif method == "backward":
		
		f_prime = ( f(point) - f(point - step) ) / step
	
	elif method == "central":
		
		f_prime = ( f(point + step) - f(point - step) ) / (2 * step)
	
	else:
		raise ValueError("This method does not exist")
	
	print("f_prime evaluated at " + str(point) + " with stepsize " + str(step) + " and method " + method + ": " + str(f_prime) )


# fx = sym.lambdify( x, sym.sin( 3*x + 2 ) )
# simple_difference_formula(fx, 1, 0.07, method="central")

def general_difference_formula(stepsizes, deriv_order):
	
	#The derivative dictates what the system of equations should look like
	if deriv_order == 0:
		prod = np.array([1,0,0])
	elif deriv_order == 1:
		prod = np.array([0,1,0])
	elif deriv_order == 2:
		prod = np.array([0,0,1])
	
	#Only works for 3x3 systems
	matrix = np.array([[1,1,1],
							 [stepsizes[0], stepsizes[1], stepsizes[2]],
							 [0.5 * stepsizes[0]**2, 0.5 * stepsizes[1]**2, 0.5 * stepsizes[2]**2 ] ] )
	
	#Return the alphas via text
	alphas = np.linalg.solve(matrix, prod)
	print("Alphas: " + str(alphas))
	print("constant k in ke/h: " + str(np.sum(np.abs(alphas) ) ) )
	
general_difference_formula([-1,0,1], 1)