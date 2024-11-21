import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import time

# Question 1
# solve eqs 8 and 9
# BC: m(r=0) = 0 , rho(r=0) = rho_c
# intergrate from 0 (or as close as possible) to rs such that rho(rs) = 0

# 10 inital values of rho to calculate solutions for (given)
rho_c = np.linspace(1/10.,2.5*(10**6),10)  # already made unitless (p/p0)
initial_m = 0.
initial_r = 1e-10

# function which returns equations 8 and 9 to be solved with solve_ivp
def diffEQNs(r, f):
   