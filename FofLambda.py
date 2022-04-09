import os
import numpy as np
import pandas as pd
from standards import Standards as s
import tools as tools

lam_vals = 5

# Edit dydr to account for gravitational scaling
def dydr_gscaled(r,y) :
    """Defining the 5 ODEs to solve using RK4. Takes radius(r) and y as inputs where y is a 5 column
       matrix representing [rho,Temp,Mass,Luminosity,OpticalDepth]

       Edited to account for gravitational scaling based on lambda
    """
    rho = y[0]
    T = y[1]
    M = y[2]
    L = y[3]
    tau = y[4]

    # The five ODEs we are trying to solve (equation 2)
    dydr = np.zeros(5)
    dydr[0] = tools.drho_dr_scaled(rho, T, r, M,L, lam_vals)
    dydr[1] = tools.dT_dr_scaled(rho, T, r, L, M, lam_vals)
    dydr[2] = tools.dM_dr(rho, r)
    dydr[3] = tools.dL_dr(rho , T, r)
    dydr[4] = tools.dtau_dr(rho, T)

    return dydr

while lam_vals >= 3:
    exec(open("generate_star_gravmod.py").read())
    lam_vals -= 1
