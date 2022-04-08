# Main code for generating the main sequence stars

import os
import numpy as np
import pandas as pd
from standards import Standards as s
import tools as tools
import scipy.integrate as int


def solveODEs(r_init, r_fin, y0, steps):
    '''Function to solve the ODEs (equation 2 in the project description)'''

    r_values = np.linspace(r_init, r_fin, steps)
    t_span = (r_init, r_fin)

    solutions = int.solve_ivp(tools.dydr,t_span,y0,'RK45', t_eval=r_values) # outputs an array 'y' which contain the solutions to the ODEs which are described dydr from tools.py

    rho = solutions.y[0]
    T = solutions.y[1]
    M = solutions.y[2]
    L = solutions.y[3]
    #tau = solutions.y[4]
    r = solutions.t # t here mean the 'time points' which in this case are the radius values or rvals

    #print(solutions.t) for testing

    return (r,rho,T,M,L)


# Setting the Initial conditions
r = 1 # starting radius for integration
steps = h = 1000 # inital step size

rho_c = 77750 # (in kg/m^3) we will vary rho_c to satisfy the surface boundary condition on the luminosity and surface temperature
Tc = 10486 # (core temperature in K) we will choose Tc to be the MS paramter
M = (4/3)*np.pi*(r**3)*rho_c # calculating the intial mass using the initial radius and densities specified
L = (4/3)*np.pi*(r**3)*rho_c*tools.epsilon(rho_c, Tc) # calculating the intial luminosity using the initial radius and densities specified
tau = tools.kappa(rho_c, Tc)*rho_c*r # calculating the intial tau using the initial radius and densities specified
del_tau = tools.del_tau(rho_c,Tc,r,M,L) # intial value for the opacity proxy

del_tau_threshold = 1e-5 # threshold value for the opacity proxy
MassLimit = s.Msun*(1e3) # mass limit for integration

y = np.zeros(5)
y = rho_c, Tc, M, L, tau

r_values = [r] # will store all radius values
y_values = [y] # will store values for rho, T, M, l and tau respectively

# ------------------------
# Section 2.2.1 (Step 2)
# ------------------------
# Now we account for the two main complications which arises from the behaviour near r = 0 and the practical determination of R_star

# We want to integrate the ODEs until del_tau (defined in eqution 16) is << 1 and we introdue a mass limit. We stop when M > 10^3 solar masses
while (np.any(del_tau) >  del_tau_threshold) and (M < MassLimit):

    y5,hnew = tools.RKF45Method_adaptive(tools.dydr,r,y,h)

    r_values.append(r + h)
    y_values.append(y5)

    # Channging variales around for next iteration of the 'while' loop
    # I.e., hnew becomes h for the next iteration and y5 becomes y for the next iteration
    h = hnew
    y = y5
    r += h

    rho = y[0]
    T = y[1]
    M = y[2]
    L = y[3]
    tau = y[4]

    del_tau = tools.del_tau(rho,T, r, M, L )
