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
r_initial = 1
r_final = 1e5
steps = 1000

rho_c = 77750 # (in kg/m^3) we will vary rho_c to satisfy the surface boundary condition on the luminosity and surface temperature
Tc = 2e7 # (core temperature in K) we will choose Tc to be the MS paramter
M0 = (4/3)*np.pi*(r_initial**3)*rho_c
L0 = (4/3)*np.pi*(r_initial**3)*rho_c*tools.epsilon(rho_c, Tc)
tau0 = tools.kappa(rho_c, Tc)*rho_c*r_initial

y0 = np.zeros(5)
y0 = rho_c,Tc,M0,L0,tau0

del_tau = tools.del_tau(rho_c,Tc,r_initial,M0,L0) # intial value for del_tau
M_vals = [M0]

test = tools.RK4Method(r_initial,r_final,y0,steps)



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
