import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

# Question 1
# solve eqs 8 and 9
# BC: m(r=0) = 0 , rho(r=0) = rho_c
# intergrate from 0 (or as close as possible) to rs such that rho(rs) = 0

# given values

# p0 = (9.47*(10**5))*(u_e)      # [g cm^-3]
# R0 = (7.72*(10**8))/(u_e)      # [cm]
# M0 = (5.67*(10**33))/(u_e**2)  # [g]

# 10 inital values of rho to calculate solutions for (given)
rho_c = np.linspace(1/10.,2.5*(10**6),10)  # already made unitless (p/p0)
initial_m = 0.
initial_r = 1e-10

# function that returns dimensionless quantitues for r, m, and rho
def dimensionless(r,m,rho):
   u_e = 2.
   p0 = (9.47*(10**5))*(u_e)      # [g cm^-3]
   R0 = (7.72*(10**8))/(u_e)      # [cm]
   M0 = (5.67*(10**33))/(u_e**2)  # [g]
   
   R = np.abs(r/R0)
   M = np.abs(m/M0)
   Rho = np.abs(rho/p0)
   return [R,M,Rho]

# function that calculates gamma for an input value of rho
def gammaX(rho):
   x = dimensionless(0,0,rho)[2]**(1/3)
   gamma = (x**2)/(3*(1+x**2)**(0.5))
   return gamma

# function which returns equations 8 and 9 to be solved with solve_ivp
def diffEQNs(r, f):
   '''r is independent variable,
   rho and m are dependent variables coorespoinding to indices [0] and [1] of f
   and must be unitless inputs (rho = rho/p0 , m = m/m0),
   '''
   
   
   # setting each ode to an element of f
   m , rho = f
   
   # defining gamma
   gam = gammaX(rho)

   # defining each ode
   drho_dr = (-m*rho)/(gam*(r**2))
   dm_dr = (r**2)*rho

   # returning eqns as appended array for easier operation
   out = np.append(dm_dr,drho_dr)
   return out

f0 = [initial_m,rho_c[0]]            # m(r=0) = 0 | rho(r=0) = rho_c 

soln = integrate.solve_ivp(diffEQNs, (initial_r,100), f0)

plt.plot(soln.t, soln.y[0,:], 'r', soln.t, soln.y[1,:], 'g')
plt.xlabel('Radius r')
plt.ylabel('Red: dmdr soln , Green: drho/dr soln')
plt.show()