import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import pandas as pd
import time                #used for time.sleep() to debug outputs from poor logic


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
   soln2 = integrate.solve_ivp(odes,[initial_radius,100],f0,t_eval=t_eval,events=zerodensity,rtol=1e-8,atol=1e-10)
   
   # extracting the values of mass and radius, and giving them units
   plotting_mass = np.multiply(M0,soln2.y[0,:])
   plotting_radius = np.multiply(R0,soln2.t[:])

   # plotting mass on x, radius on y
   # plt.plot(plotting_mass,plotting_radius,'-')
   # plt.xlabel('Mass [g]')
   # plt.ylabel('Stellar radius [cm]')
   # plt.title('White dwarf radius in terms of mass for initial density ~'+str(np.round(initial_density[j],2)))
   # plt.show()

# My result for the Chandrasekhar limit was estimated from the plots to be
# 2.86e33. This is roughly 2 times as large compared to the result obtained by
# Kippenhahn & Weigert (1990) of 5.836/ue^2 ~= 1.459. This is likely due to
# the use of poor starting values for the central density of white dwarves
# used in the integration solution.

# Question 3

# Pick 3 values of initial central density and run solve_ivp again, choosing a
# different integration method from the one used in part one (RK45). How different
# are you results?

# initializing an empty arrays which can store the 
# solution values mass and radius for quesiton 3
solved_masses3 = np.zeros(np.size(initial_density))
solved_radii3 = np.zeros(np.size(initial_density))

# using the first 3 values of initial central density
initial_densities3 = initial_density[0:3]

for k in range(np.size(initial_densities3)):
   
   # setting the initial condition
   f0 = [0., initial_densities3[k]]
   
   # solving eqns 8 and 9 for the given initial condition, now using BDF as the integration method
   soln3 = integrate.solve_ivp(odes,[initial_radius,100],f0,t_eval=t_eval,events=zerodensity,rtol=1e-8,atol=1e-10, method='RK23')
   
   # saving the dimensionless values of mass and radius to arrays
   solved_radii3[k] = soln3.t[-1]
   solved_masses3[k] = soln3.y[0,-1]

   # comparison to be printed to the terminal
   print('Initial Central Density = '+str(initial_densities3[k]))
   print('RK45 output radius = '+str(solved_radii[k])+' | RK23 output radius = '+str(solved_radii3[k]))
   print('RK45 output mass = '+str(solved_masses[k])+' | RK23 output mass = '+str(solved_masses3[k])+'\nALL QUANTITES ARE UNITLESS\n')

# WRITTEN ANSWER TO QUESTION 3:
# The results from the method used in question 1 (RK45) when compared to those obtained using
# the method in question 3 (RK23) are very similar to one another. This is not much of a surprise
# as they are both explicit Runge-Kutta methods of iterative numerical analysis. While the solutions
# are practically identical, the RK23 method differs slightly from RK45 for the masses moreso than
# for the radii.

# Question 4

# Plot the observed data obtained by Tremblay et al. (2017) from the given file wd_mass_radius.csv.
# Your plot must include Plot these observed data with their error bars on your computed mass-radius relation, paying attention to the units. How well do the observations agree with your calculations?

# Note - Measurements are in units of the Sun's mass and radius

csv_out = pd.read_csv(r'C:\Users\natha\Desktop\UWO\2024-2025\1st Semester\Physics 3926 - Computer simulations\Python\wd_mass_radius.csv',skiprows=1)
