import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

# Question 1
# solve eqs 8 and 9
# BC: m(r=0) = 0 , rho(r=0) = rho_c
# intergrate from 0 (or as close as possible) to rs such that rho(rs) = 0

# given values
u_e = 2.
p0 = (9.47*(10**5))*(u_e)      # [cm^-3]
R0 = (7.72*(10**8))/(u_e)      # [cm]
M0 = (5.67*(10**33))/(u_e**2)  # [g]

# 10 inital values of rho to calculate solutions for (given)
rho_c = np.linspace(1/10.,2.5*(10**6),10)

# gamma
# x = (rho/p0)**2
# gam = 

# function which returns equations 8 and 9 to be solved with solve_ivp
def diffEQNs(r, f, gamma, ue):
    ''''''
    m, rho = f
    gamma = ((rho/p0)**(2/3.))/((3.*(1+(rho/p0)**(2/3.)))**(1/2.))
    drho_dr = (-m*rho)/(gamma*(r**2))
    dm_dr = (r**2)*rho 
    return drho_dr , dm_dr

