import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import pandas as pd        # used for reading the given csv file
import time                # used for time.sleep() to debug outputs from poor/problematic code


# Question 1
# solve eqs 8 and 9
# BC: m(r=0) = 0 , rho(r=0) = rho_c
# intergrate from 0 (or as close as possible) to rs such that rho(rs) = 0

# given constant
ue = 2.

# 10 inital values of rho to calculate solutions for (given)
initial_density = np.linspace(1/10.,2.5*(10**6),10)  # already made unitless (p/p0)
initial_mass = 0.
initial_radius = 1e-10

# beta function definiton
def odes(r, f):
   '''This function is used to describe the differential mass and density equations 
   of a white dwarf star, and returns the equations to be solved by solve_ivp.'''
   
   mass = f[0]   # setting mass and density as functions odes() in terms of input f
   density = f[1]

   # Logic to stop bad values of density from being used
   if density >= 0:
      x = (density)**(1/3)              # defining x
      gamma = (x**2)/(3*(1+x**2)**(1/2))   # defining gamma in terms of x

      # definining differential equations  
      dm_dr = (r**2)*(density)
      drho_dr = -(mass)*(density)/(gamma*(r**2))
      
   
      return [dm_dr,drho_dr]

# setting event condition to stop integration when density is zero
def zerodensity(r, f):
   '''Function used with solve_ivp to stop integration when the density crosses
   zero, going from a positive to a negative value.'''
   return  f[1] 

zerodensity.terminal = True
zerodensity.direction = -1  

# initializing an empty arrays which can store the solution values mass and radius
solved_masses = np.zeros(np.size(initial_density))
solved_radii = np.zeros(np.size(initial_density))

# setting evaluation time for integration
t_eval = np.linspace(initial_radius,10,10000000)

for i in range(np.size(initial_density)):
   
   # setting the initial condition
   f0 = [0., initial_density[i]]
   
   # solving eqns 8 and 9 for the given initial condition
   soln = integrate.solve_ivp(odes,[initial_radius,10],f0,t_eval=t_eval,events=zerodensity,rtol=1e-8,atol=1e-10)
   
   # saving the dimensionless values of mass and radius to arrays
   solved_radii[i] = soln.t[-1]
   solved_masses[i] = soln.y[0,-1]

   # printing the calculated mass and radius for each initial density
   print('For initial central density '+str(initial_density[i]))
   print('Calculated radius = '+str(solved_radii[i]))
   print('Calculated mass = '+str(solved_masses[i])+'\n')


# Question 2

# Transform your results from the ODE solution into physical (not dimensionless)
# units and plot R as a function of M (R on y, M on x). Can you estimate the Chandraskhar limit
# from your plot? How does it compoare to Mch = 5.836/(ue^2)?

# relative mass, radius, and density for converting to real units
M0 = 5.67e33/(ue**2)  # [g]
R0 = 7.72e8/(ue)    # [cm]
p0 = 9.74e5*(ue)  # [g/cm^3]

# similar loop as in question 2 but now for plotting
for j in range(np.size(initial_density)):
   
   # setting the initial condition for the current loop iteration
   f0 = [0., initial_density[j]]
   soln2 = integrate.solve_ivp(odes,[initial_radius,10],f0,t_eval=t_eval,events=zerodensity,rtol=1e-8,atol=1e-10)
   
   # extracting the values of mass and radius, and giving them units
   plotting_mass = np.multiply(M0,soln2.y[0,:])
   plotting_radius = np.multiply(R0,soln2.t[:])

   # plotting mass on x, radius on y
   plt.plot(plotting_mass,plotting_radius,'-',color='k')
   plt.xlabel('Mass [g]')
   plt.ylabel('Stellar radius [cm]')
   plt.title('White dwarf radius in terms of mass for initial density ~'+str(np.round(initial_density[j],2)))
   plt.show()

# WRITTEN ANSWER QUESTION 2:
# My result for the Chandrasekhar limit was estimated from the plots to be
# ~2.017 (2.86e33 [g] / M0 [g]). This is roughly 38% larger than the result obtained by
# Kippenhahn & Weigert (1990) of 5.836/ue^2 ~= 1.459. This is likely due to
# the use of inappropriate starting values for the central density of white dwarves
# used in the integration solution.

# Question 3

# Pick 3 values of initial central density and run solve_ivp again, choosing a
# different integration method from the one used in part one (RK45). How different
# are your results?

# initializing an empty arrays which can store the 
# solution values mass and radius for quesiton 3
solved_masses3 = np.zeros(np.size(initial_density))
solved_radii3 = np.zeros(np.size(initial_density))

# using the first 3 values of initial central density
initial_densities3 = initial_density[0:3]

for k in range(np.size(initial_densities3)):
   
   # setting the initial condition
   f0 = [0., initial_densities3[k]]
   
   # solving eqns 8 and 9 for the given initial condition, now using RK23 as the integration method
   soln3 = integrate.solve_ivp(odes,[initial_radius,10],f0,t_eval=t_eval,events=zerodensity,rtol=1e-8,atol=1e-10, method='RK23')
   
   # saving the dimensionless values of mass and radius to arrays
   solved_radii3[k] = soln3.t[-1]
   solved_masses3[k] = soln3.y[0,-1]

   # comparison to be printed to the terminal
   print('Initial Central Density = '+str(f0[1]))
   print('RK45 output radius = '+str(solved_radii[k])+' | RK23 output radius = '+str(solved_radii3[k]))
   print('RK45 output mass = '+str(solved_masses[k])+' | RK23 output mass = '+str(solved_masses3[k])+'\nALL QUANTITES ARE UNITLESS\n')

# WRITTEN ANSWER TO QUESTION 3:
# The results from the method used in question 1 (RK45) when compared to those obtained using
# the method in question 3 (RK23) are very similar. This is not much of a surprise as they are 
# both explicit Runge-Kutta methods of iterative numerical analysis. While the solutions
# are practically identical, the RK23 method differs slightly from RK45 for the masses more so than
# for the radii.

# Question 4

# Plot the observed data obtained by Tremblay et al. (2017) from the given file wd_mass_radius.csv.
# Your plot must include Plot these observed data with their error bars on your computed mass-radius relation, paying attention to the units. How well do the observations agree with your calculations?

# Note - Measurements are in units of the Sun's mass and radius

# sun constants
sun_mass = 1.989e33      #[g]
sun_radius = 69.634e9    # [cm]

# reading the csv file from its raw string path and assigning it to csv_out
csv_out = pd.read_csv(r'C:\Users\natha\Desktop\UWO\2024-2025\1st Semester\Physics 3926 - Computer simulations\Python\wd_mass_radius.csv')

# isolating only the needed data into easily manipulable np array
csv_data = np.array(csv_out)

# array data slicing for plotting with unit conversion included
Msun_wd = csv_data[:,0]        # white dwarf sun masses 
Msun_unc_wd = csv_data[:,1]    # white dwarf sun masses uncertainty
Rsun_wd = csv_data[:,2]        # white dwarf sun radii
Rsun_unc_wd = csv_data[:,3]    # white dwarf sun radii uncertainty

# testing our mass-radius relation for the following initial central densities (guessed until I found some that looked good)
test_densities=[0.5,1,1.5,2,2.5,3,3.5,4,4.5,5]

for l in range(np.size(test_densities)):
   
   # setting the initial condition for the current loop iteration
   f0 = [0., test_densities[l]]
   soln4 = integrate.solve_ivp(odes,[initial_radius,10],f0,t_eval=t_eval,events=zerodensity,rtol=1e-8,atol=1e-10)
   
   # extracting the values of mass and radius, and giving them units
   plotting_mass4 = np.multiply(M0/sun_mass,soln4.y[0,:])
   plotting_radius4 = np.multiply(R0/sun_radius,soln4.t[:])
   
   # loop for plotting all wd values with bidirectional error bars
   for m in range(np.size(Msun_wd)):
      # plotting each radius-mass pair 
      x_mass = Msun_wd[m]
      x_unc = Msun_unc_wd[m]
      y_radius = Rsun_wd[m]
      y_unc = Rsun_unc_wd[m]

      # plotting each radius mass pair with mass and radius error bars
      plt.plot(x_mass,y_radius,'--')
      plt.errorbar(x_mass,y_radius,xerr=x_unc,yerr=y_unc,ecolor='b',elinewidth=0.5)
   
   # plot formatting
   # plotting mass on x, radius on y
   plt.plot(plotting_mass4,plotting_radius4,'-',color='k')
   plt.xlabel('Solar masses')
   plt.ylabel('Solar radii')
   plt.title('White dwarf solar radii in terms of solar mass for initial density ~'+str(np.round(f0[1],2)))
   plt.show()

# WRITTEN ANSWER TO QUESTION 4:
# After experimenting with varying values for initial central density, I found that the a range
# from 0.5 to 5 with evenly spaced intervals of 0.5 yielded results for final mass and radius
# similar to those observed by Tremblay et al. 2017. However, when using the values for initial
# central density from question 1, only the first initial central density of 0.1 came close to
# the Tremblay et al. 2017 observation values with all others having much higher masses and
# significantly lower radii.