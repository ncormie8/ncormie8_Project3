import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import time

# Question 1
# solve eqs 8 and 9
# BC: m(r=0) = 0 , rho(r=0) = rho_c
# intergrate from 0 (or as close as possible) to rs such that rho(rs) = 0

ue = 2.

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
def zerodensity(r, f):
   return  f[1] 

zerodensity.terminal = True
zerodensity.direction = -1  


# initializing an empty arrays which can store the solution values mass and radius
solved_masses = np.zeros(np.size(initial_density))
solved_radii = np.zeros(np.size(initial_density))

# setting evaluation time for integration
t_eval = np.linspace(initial_radius,100,10000000)

for i in range(np.size(initial_density)):
   
   # setting the initial condition
   f0 = [0., initial_density[i]]
   
   # solving eqns 8 and 9 for the given initial condition
   soln = integrate.solve_ivp(odes,[initial_radius,100],f0,t_eval=t_eval,events=zerodensity,rtol=1e-8,atol=1e-10)
   
   # saving the dimensionless values of mass and radius to arrays
   solved_radii[i] = soln.t[-1]
   solved_masses[i] = soln.y[0,-1]


# Question 2

# Transform your results from the ODE solution into physical (not dimensionless)
# units and plot R as a function of M (R on y, M on x). Can you estimate the Chandraskhar limit
# from your plot? How does it compoare to Mch = 5.836/(ue^2)?

# relative mass, radius, and density for converting to real units
M0 = 5.67e33/(ue**2)  # [g]
R0 = 7.72e8/(ue)    # [cm]
p0 = 9.74e5*(ue)  # [g/cm^3]

# similar loop as in question 2
for j in range(np.size(initial_density)):
   
   # setting the initial condition for the current loop iteration
   f0 = [0., initial_density[j]]
   soln = integrate.solve_ivp(odes,[initial_radius,100],f0,t_eval=t_eval,events=zerodensity,rtol=1e-8,atol=1e-10)
   
   # extracting the values of mass and radius, and giving them units
   plotting_mass = np.multiply(M0,soln.y[0,:])
   plotting_radius = np.multiply(R0,soln.t[:])

   # plotting mass on x, radius on y
   plt.plot(plotting_mass,plotting_radius,'-')
   plt.xlabel('Mass [g]')
   plt.ylabel('Stellar radius [cm]')
   plt.title('White dwarf radius in terms of mass for initial density')
   plt.show()