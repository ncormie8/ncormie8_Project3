import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import time

# Question 1
# solve eqs 8 and 9
# BC: m(r=0) = 0 , rho(r=0) = rho_c
# intergrate from 0 (or as close as possible) to rs such that rho(rs) = 0

# M0 = 5.67e33/(2.**2)  # [g]
# R0 = 7.72e8/2.    # [cm]
# p0 = 9.74e5*2.  # [g/cm^3]

# 10 inital values of rho to calculate solutions for (given)
initial_density = np.linspace(1/10.,2.5*(10**6),10)  # already made unitless (p/p0)
initial_mass = 0.
initial_radius = 1e-10

# beta function definiton
def odes(r, f):
   mass = f[0]   # setting mass and density as functions odes() in terms of r
   density = f[1]

   # Logic to stop bad values of density being input
   if density >= 0:
      x = (density)**(1/3)              # defining x
      gamma = (x**2)/(3*(1+x**2)**(1/2))   # defining gamma in terms of x

      # definining differential equations  
      drho_dr = -(mass)*(density)/(gamma*(r**2))
      dm_dr = (r**2)*(density)
   
      return [dm_dr,drho_dr]

# setting event condition to stop integration when density is zero


# initializing an empty arrays which can store the solution values mass and radius
solved_masses = np.zeros(np.size(initial_density))
solved_radii = np.zeros(np.size(initial_density))

# setting evaluation time for integration
t_eval = np.linspace(initial_radius,10,1000000000)

for i in range(np.size(initial_density)):
   def zerodensity(r, f):
      return f[0]
   
   zerodensity.terminal = True
   zerodensity.direction = -1  

   # setting the initial condition
   f0 = [0., initial_density[i]]
   soln = integrate.solve_ivp(odes,[initial_radius,10],f0,events=zerodensity,rtol=1e-8,atol=1e-10)
   print(soln.t_events[0])
   print(soln.y_events[0])
   time.sleep(1)