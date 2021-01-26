# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 10:40:54 2020

@author: Matthijs Schrage
"""
import numpy as np
import sympy as sym

x = sym.Symbol("x")
y = sym.Symbol("y")
variables = [x, y]

#Some function to help us with syntax
def sin(x):
	return sym.sin(x)

def cos(x):
	return sym.cos(x)

def sqrt(x):
	return sym.sqrt(x)

#Define a function that takes in a certain formula and lambifies it and it's first two derivatives
def make_function(function):
	
	first_derivative = sym.diff(function, x)
	second_derivative = sym.diff(first_derivative, x)
	
	f = sym.lambdify(x, function, 'numpy')
	f_prime = sym.lambdify(x, first_derivative, 'numpy')
	f_double_prime = sym.lambdify(x, second_derivative, 'numpy')
	
	return f, f_prime, f_double_prime

def bisection(function, interval, error):
	
	# Check if there exists at least one zero solution on the interval
	assert interval[0] < interval[1]
	assert function(interval[0]) * function(interval[1]) < 0
	
	#This is the formula to calculate the number of iterations we need to take to be inside a certain specfified error
	iterations = np.log2( ( interval[1] - interval[0] ) / error)
	
	print("Nr of iterations in bisection: " + str(int(np.floor(iterations) ) ) )
	print("Maximum error in bisection: " + str(error))
	
	#Iterate over the interval until we reach the required number of iterations
	for n in range(int(np.floor(iterations))):
		midpoint = 0.5 * (interval[0] + interval[1])
		
	# 	products = np.array([function(interval[0]) * function(midpoint), function(midpoint) * function(interval[1])])
		half = np.argmin([function(interval[0]) * function(midpoint), function(midpoint) * function(interval[1])])
		
		#Return the half that has the negative sign
		if half == 0:
			interval = np.array([interval[0], midpoint])
		elif half == 1:
			interval = np.array([midpoint, interval[1]])
		else:
			raise ValueError("Is ff iets mis gegaan in bisection")
	
	print("p lies somewhere between: " + str(interval))
	
	return interval

# fx = make_function( x * ( x**2 + 5*x + 6 ) * (x + 3.5) )[0]
# fx = make_function( x - 0.5 * sin(x) - np.pi )[0]
# bisection(fx, [-10, 10], 0.01)

def error_sequence(labda, alpha, e0, nr_of_iterations):
	
	#Check value of alpha
	assert 1 <= alpha >= 2
	
	#Initiate the error list
	errors = []
	error = e0
	errors.append(error)
	
	#Calculate error for the number if iterations specified
	for _ in range(nr_of_iterations):
		
		error = (np.abs(error) ** alpha) * labda
		
		errors.append(error)
	
	print("Error sequence: " + str(errors))
	
	return errors

def fixed_point_iteration1(h, iterations, p0):
	
	approximations = []
	approximations.append(p0)
	
	p = 5 ** (1/4)
	
	print("\np = " + str(p))
	
	# This is the equation used for the fixed point iteration
	# CAN BE DIFFERENT, CHECK VERY CAREFULLY
	gx, gx_prime = make_function(x + h(x) * ( (x**4) - 5 ))[:2]
	
	print("|g'(p)| = " + str( np.abs( gx_prime(p) ) ) )
	
	for _ in range(iterations):
		
		p0 = gx(p0)
		
		approximations.append(p0)
	
	print("Approximations fixed_point: \n" + str(approximations) )
	
	return p0

h1 = make_function( -1 / (x ** 5) )[0]
h2 = make_function( -1 / (x ** 2) )[0]
h3 = make_function( -1 / (4 * x ** 3) )[0]

fixed_point_iteration1(h1, 2, 1.5)
fixed_point_iteration1(h2, 2, 1.5)
fixed_point_iteration1(h3, 2, 1.5)

def fixed_point_iteration2(gx, gx_prime, p, iterations):
	
	approximations = []
	approximations.append(p)
	
	print("|g'(p)| = " + str( np.abs( gx_prime(p) ) ) )
	
	for _ in range(iterations):
		
		p = gx(p)
		
		approximations.append(p)
	
	print("Approximations fixed_point: \n" + str(approximations) )

# gx1, gx1_prime = make_function( ( np.pi * np.sqrt(2) / 8 ) / sin(x) )[:2]
# gx2, gx2_prime = make_function( x - ( x * sin(x) - (  np.pi * np.sqrt(2) / 8 ) ) )[:2]

# gx1, gx1_prime = make_function( sqrt(x + 6) )[:2]
# gx2, gx2_prime = make_function( ( (x**2) + 6 ) / (2 * x - 1) )[:2]

# fixed_point_iteration2(gx1, gx1_prime, 3, 1)
# fixed_point_iteration2(gx2, gx2_prime, 3, 1)

def newton_raphson(fx, fx_prime, x, iterations, p):
	
	approximations = []
	approximations.append(x)
	
	for _ in range(iterations):
		
		x = x - ( fx(x) / fx_prime(x) )
		approximations.append(x)
	
	print("Approximations newton_raphson: \n" + str(approximations) )
	
	errors = np.abs( np.array(approximations) - p )
	
	print("Errors newton_raphson: \n" + str(errors))

# fx1, fx1_prime = make_function( 1 + (2*x**2)/(3*x-2) )[:2]

# newton_raphson(fx1, fx1_prime, 1/4, 1, np.sqrt(5))

def jacobian(system, p, iterations=2):
	
	J = sym.zeros(len(system),len(variables))
	
	for i, fs in enumerate(system):
		for j, s in enumerate(variables):
			J[i,j] = sym.diff(fs, s)
	
	print("Jacobian matrix in formula form: \n" + str(np.array(J)))
	
	for _ in range(iterations):
		
		J_p = np.array(J.subs( [ ( x, p[0] ), ( y, p[1] ) ] ) ).astype("float64")
		F_p = np.array(system.subs( [ ( x, p[0] ), ( y, p[1] ) ] ) ).astype("float64")
		h = np.linalg.solve(J_p,- F_p)
		p = (p + h.T)[0]
		
		print("Jacobian: \n" + str(J_p))
		
		print("-F_pn: \n" + str(-F_p))
		
		print("hn: \n" + str(h))
		
		print("||hn||: " + str(np.linalg.norm(h)))
		
		print("New p: \n" + str(p))
		
		print("__________")

# z = sym.Matrix([ (1/np.pi)*sin(np.pi * x), cos(np.pi*x) + 2*y ])

# jacobian(z, np.array([1,1]), iterations=1)
