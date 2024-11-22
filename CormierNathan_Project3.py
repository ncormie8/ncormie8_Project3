import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import time

# Question 1
# solve eqs 8 and 9
# BC: m(r=0) = 0 , rho(r=0) = rho_c
# intergrate from 0 (or as close as possible) to rs such that rho(rs) = 0

# 10 inital values of rho to calculate solutions for (given)
initial_density = np.linspace(1/10.,2.5*(10**6),10)  # already made unitless (p/p0)
initial_mass = 0.
initial_radius = 1e-10

# beta function definiton
def odes(r, f, gamma, mew_e):
   # constants
   M0 = 5.67e33/(mew_e**2)  # [g]
   R0 = 7.72e8/mew_e    # [cm]
   p0 = 9.74e5*mew_e  # [g/cm^3]

   mass, density = f   # setting mass and density as functions odes() in terms of r

   # definining differential equations  
   drho_dr = -(mass/M0)*(density/p0)/gamma*((r/R0)**2)
   dm_dr = ((r/R0)**2)*(density/p0)

   return [dm_dr,drho_dr]



