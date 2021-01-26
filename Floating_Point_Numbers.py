# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 14:05:16 2020

@author: matth
"""
from decimal import Decimal
import numpy as np

u = 2
v = 3
w = 5
precision = 4

# round(float(number), precision)

print( ( round(float(round(float(355/113), 9) - round(float(np.pi), 9)), 9) - (355/113 - np.pi) ) / (355/113 - np.pi) )

# print( round(float(2/3), 4))


