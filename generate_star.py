# Main code for generating the main sequence stars

import os
import numpy as np
import pandas as pd
from standards import Standards as s
import tools as tools
import scipy.integrate as int

# Initial conditions
M0 = 0
L0 = 0
rho0 = rho_c = 50 # we will vary rho_c to satisfy the surface boundary condition on the luminosity and surface temperature
T0 = T_c = 10000 # we will choose Tc to be the MS paramter
y0 = np.zeros(4)
y0 = rho0,T0,M0,L0

r_init = 0.0000001 # intial value for the range of radii
r_fin = 10000000 # final values for the range of radii
steps = 10000 # number of steps we want the integration to go through
del_tau = 1 # intial value for del_tau
M_vals = [0]

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


y0 = rho0,T0,M0,L0
test = solveODEs(0.01, 1000, y0, 100)



# # ------------------------
# # Section 2.2.1 (Step 2)
# # ------------------------
#
# # Now we account for the two main complications which arises from the behaviour near r = 0 and the practical determination of R_star
#
# # We want to integrate the ODEs until del_tau (defined in eqution 16) is << 1 and we introdue a mass limit. We stop when M > 10^3 solar masses
# if (abs(del_tau) > 1e-5) and (max(M_vals) < 1e3*s.Msun):
#
#     solns = solveODEs(r_init, r_fin, y0, steps) # calling the solveODEs function to solve the ODEs using RK45
#
#     # Storing the results of the ODE integration
#     r_vals = solns[0]
#     rho_vals = solns[1]
#     T_vals = solns[2]
#     M_vals = solns[3]
#     L_vals = solns[4]
#     tau_vals = solns[5]
#
#     # Recalculating the value for del_tau (as defined by equation 16)
