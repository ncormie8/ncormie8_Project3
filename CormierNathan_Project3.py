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

first_rho = 1/10.
last_rho = 2.5*(10**6)
rho_c = np.linspace(first_rho,last_rho,10)
print(rho_c)
a=1
def diffEQNs():
    
    return a