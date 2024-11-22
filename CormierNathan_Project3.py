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
def odes(r, f, mew_e):
   mass = f[0]   # setting mass and density as functions odes() in terms of r
   density = f[1]

   x = (density)**(1/3)              # defining x
   gamma = (x**2)/(3*(1+x**2)**(1/2))   # defining gamma in terms of x

   # definining differential equations  
   drho_dr = -(mass)*(density)/gamma*((r)**2)
   dm_dr = ((r)**2)*(density)
   
   return [dm_dr,drho_dr]

soln = integrate.solve_ivp(odes,[initial_radius,100],[0.,initial_density[0]],args=(2.0,))
print(soln)

plt.plot(soln.t, soln.y[0,:], 'r', soln.t, soln.y[1,:], 'g')
plt.xlabel('Radius r')
plt.ylabel('Red = m(r) ||||| Green = rho(r)')
plt.show()