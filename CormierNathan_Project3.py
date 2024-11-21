import numpy as np
from scipy import integrate
from matplotlib import pyplot

# Question 1
# solve eqs 8 and 9
# BC: m(r=0) = 0 , rho(r=0) = rho_c
# intergrate from 0 (or as close as possible) to rs such that rho(rs) = 0

# given values
u_e = 2

# 10 inital values of rho to calculate solutions for (given)
rho_c = np.linspace(1/10,2.5*(10^6),10)
a=1

def diffEQNs():
    return a