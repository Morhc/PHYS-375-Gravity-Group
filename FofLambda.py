import os
import numpy as np
import pandas as pd
from standards import Standards as s
import tools as tools

lam_val = 0

# Edit dydr to account for gravitational scaling
def dydr_gscaled(r,y) :
    """Defining the 5 ODEs to solve using RK4. Takes radius(r) and y as inputs where y is a 5 column
        matrix representing [rho, Temp, Mass, Luminosity, OpticalDepth].
        Edited to account for gravitational scaling based on lambda.
    INPUTS:
        r - The radius per step of the star.
        y - The values being integrated (rho, T, M, L, tau).
    OUTPUTS:
        dydr - The integrated values for rho, T, M, L, tau.
    """

    rho, T, M, L, tau = y


    # The five ODEs we are trying to solve (equation 2)
    dydr = np.zeros(5)
    dydr[0] = tools.drho_dr_scaled(rho, T, r, L, M, lam_val)
    dydr[1] = tools.dT_dr_scaled(rho, T, r, L, M, lam_val)
    dydr[2] = tools.dM_dr(rho, r)
    dydr[3] = tools.dL_dr(rho, T, r)
    dydr[4] = tools.dtau_dr(rho, T)

    return dydr
