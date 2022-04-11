 #script: tools.py
#purpose: To have a hub for all useful functions.

import os
import numpy as np
import pandas as pd
from standards import Standards as s

def load_data():
    """Loads in Broderick's data. Warning, the data is in CGMS units, NOT SI units.
    OUTPUTS:
        highmass - A DataFrame of the high mass star's information.
        lowmass - A DataFrame of the low mass star's information.
        summary - A summary of the star information.
    """

    #all the data for a single high mass star
    highmass = pd.read_csv(os.path.join(s.data_folder, 'highmass_star.txt'), sep='\s+', header=0, skiprows=[1,2])
    #all the data for a single low mass star
    lowmass = pd.read_csv(os.path.join(s.data_folder, 'lowmass_star.txt'), sep='\s+', header=0, skiprows=[1,2])
    #a summary of the two high and low mass stars
    summary = pd.read_csv(os.path.join(s.data_folder, 'star_summary.txt'), sep='\s+', header=0, skiprows=[1,2])

    #note: there appears to be two dL/dr columns with the same data

    return highmass, lowmass, summary

def dydr(r,y) :
    """Defining the 5 ODEs to solve using RK4. Takes radius(r) and y as inputs where y is a 5 column
       matrix representing [rho,Temp,Mass,Luminosity,OpticalDepth, lambda].

       DEFUNCT.
    """
    rho, T, M, L, tau = y

    # The five ODEs we are trying to solve (equation 2)
    dydr = np.array([drho_dr(rho, T, r, L, M), dT_dr(rho, T, r, L, M), dM_dr(rho, r), dL_dr(rho, T, r), dtau_dr(rho, T)])

    return dydr

def epsilon(rho, T):
    """Calculates epsilon using Equations 8 and 9 from the project description.
    INPUTS:
        rho - The density of the star at some radius.
        T - The temperature of the star at some radius.
    OUPUTS:
        eps - The total energy generation of the star, in W/kg.
    """

    rho5, T6 = rho/1e5, T/1e6

    epp = (1.07e-7)*rho5*(s.X**2)*np.power(T6, 4) #PP-chain, in W/kg
    ecno = (8.24e-26)*rho5*s.X*s.X_CNO*np.power(T6, 19.9) #CNO cycle, in W/kg

    eps = epp + ecno

    return eps

def kappa_es():
    """Calculates the Rosseland mean opacity for electron scattering via Equation 10 from the project
    description.
    OUTPUTS:
        kappaes - The Rosseland mean opactiy for electron scattering, in m^2/kg.
    """

    kappaes = 0.02*(1+s.X)

    return kappaes

def kappa_ff(rho, T):
    """Calculates the Rosseland mean opacity for a Kramers-like aproximation to free-free interaction
    via Equation 11 in the project description.
    INPUTS:
        rho - The density of the star at some radius.
        T - The temperature of the star at some radius.
    OUTPUTS:
        kappaff - The Rosseland mean opactiy for free-free, in m^2/kg.
    """

    rho3 = rho/1e3

    kappaff = (1e24)*(s.Z+0.0001)*np.power(rho3, 0.7)*np.power(T, -3.5)

    return kappaff

def kappa_Hminus(rho, T):
    """Calculates the Rosseland mean opacity due to Hminus at sufficiently low temperatures near surface
    via Equation 12 in the project description.
    INPUTS:
        rho - The density of the star at some radius.
        T - The temperature of the star at some radius.
    OUTPUTS:
        kappaHm - The Rosseland mean opactiy for Hminus, in m^2/kg.
    """

    rho3 = rho/1e3
    kappaHm = (2.5e-32)*(s.Z/0.02)*np.power(rho3, 0.5)*np.power(T, 9)

    return kappaHm

def kappa(rho, T):
    """Calculates the radiative opacity via Equation 14 in the project description.
    INPUTS:
        rho - The density of the star at some radius.
        T - The temperature of the star at some radius.
    OUTPUTS:
        kappa - The radiative opacity, in m^2/kg.

    """

    kappaes = kappa_es()
    kappaff = kappa_ff(rho, T)
    kappaHm = kappa_Hminus(rho, T)

    kappa = np.power(1/kappaHm + 1/max(kappaes, kappaff), -1)

    return kappa

def P_degen(rho):
    """Calculates the non-relativistic degenerate pressure, from Equation 5 in the project description.
    INPUTS:
        rho - The density of the star at some radius.
    OUTPUTS:
        Pdeg - The degenerate pressure, in Pa.
    """

    Pdeg = np.power(3*(np.pi**2), 2/3)*np.power(s.hbar,2)*np.power(rho/s.mp, 5/3)/(5*s.me)

    return Pdeg

def P_ideal(rho,T):
    """Calculates the pressure from an ideal gas, from Equation 5 in the project description.
    INPUTS:
        rho - The density of the star at some radius.
        T - The temperature of the star at some radius.
    OUTPUTS:
        Pgas - The ideal pressure, in Pa.
    """

    Pgas = rho*s.k*T/(s.mu*s.mp)

    return Pgas

def P_photongas(T):
    """Calculates the pressure from the photon gas, from Equation 5 in the project description.
    INPUTS:
        T - The temperature of the star at some radius.
    OUTPUTS:
        Ppho - The photon pressure, in Pa.
    """

    Ppho = s.a*np.power(T, 4)/3

    return Ppho

def pressure(rho, T):
    """Calculates the pressure via Equation 5 in the projection description.
    INPUTS:
        rho - The density of the star at some radius.
        T - The temperature of the star at some radius.
    OUTPUTS:
        P - The pressure, in Pa.
        Pdeg - The degenerate pressure, in Pa.
        Pgas - The ideal pressure, in Pa.
        Ppgo - The photon pressure, in Pa.
    """

    Pdeg = P_degen(rho)
    Pgas = P_ideal(rho,T)
    Ppho = P_photongas(T)

    P = Pdeg + Pgas + Ppho

    return P

def dP_dT(rho, T):
    """Calculates the partial derivative of the pressure wrt temperature via Equation 7 in the project
    description.
    INPUTS:
        rho - The density of the star at some radius.
        T - The temperature of the star at some radius.
    OUTPUTS:
        dP_dT - The partial derivative of the pressure wrt temperature, in Pa/K.
    """

    dP_dT = rho*s.k/(s.mu*s.mp) + s.a*np.power(T, 3)*4/3

    return dP_dT

def dP_drho(rho, T):
    """Calculates the partial derivative of the pressure wrt rho via Equation 7 in the project
    description.
    INPUTS:
        rho - The density of the star at some radius.
        T - The temperature of the star at some radius.
    OUTPUTS:
        dP_drho - The partial derivative of the pressure wrt rho, in Pa m^3/kg [J/kg].
    """

    A = np.power(3*(np.pi**2), 2/3)*np.power(s.hbar,2)*np.power(rho/s.mp, 2/3)/(3*s.me*s.mp)
    B = s.k*T/(s.mu*s.mp)

    dP_drho = A + B

    return dP_drho

def dlnP_dlnT(P, T):
    """Calculates the derivative of lnP wrt lnT. NOTE: WRONG
    INPUTs:
        P - The pressure at some radius.
        T - The temperature at some radius.
    OUTPUTS:
        dlnP_dlnT - The derivative.
    """

    dlnP_dlnT = np.log(P)/np.log(T)

    return dlnP_dlnT

def dtau_dr(rho, T):
    """Calculates the opatical depth via Equation 2.
    INPUTS:
        rho - The density at some radius.
        T - The temperature at some radius.
    OUTPUTS:
        tau - The optical depth.
    """

    tau = kappa(rho, T)*rho

    return tau

def dL_dr_PP(rho, T, r):
    """Calculates the luminosity via Equation 2, but replaces epsilon with the equation for eps_pp.
    INPUTS:
        rho - The density at some radius.
        T - The temperature at some radius.
        r - The radius.
    OUTPUTS:
        dL_dr_pp - The luminosity from the PP chain at some radius, in W/m.
    """

    rho5, T6 = rho/1e5, T/1e6

    epp = (1.07e-7)*rho5*(s.X**2)*np.power(T6, 4) #PP-chain, in W/kg

    dL_dr_pp = 4*np.pi*np.power(r, 2)*rho*epp

    return dL_dr_pp

def dL_dr_CNO(rho, T, r):
    """Calculates the luminosity via Equation 2, but replaces epsilon with the equation for eps_cno.
    INPUTS:
        rho - The density at some radius.
        T - The temperature at some radius.
        r - The radius.
    OUTPUTS:
        dL_dr_cno - The luminosity from the CNO cycle at some radius, in W/m.
    """

    rho5, T6 = rho/1e5, T/1e6

    ecno = (8.24e-26)*rho5*s.X*s.X_CNO*np.power(T6, 19.9) #CNO cycle, in W/kg

    dL_dr_cno = 4*np.pi*np.power(r, 2)*rho*ecno


    return dL_dr_cno

def dL_dr(rho, T, r):
    """Calculates the luminosity via Equation 2.
    INPUTS:
        rho - The density at some radius.
        T - The temperature at some radius.
        r - The radius.
    OUTPUTS:
        dL_dr - The luminosity at some radius, in W/m.
    """

    dL_dr = 4*np.pi*(r**2)*rho*epsilon(rho, T)

    return dL_dr

def dM_dr(rho, r):
    """Calculates the mass via Equation 2.
    INPUTS:
        rho - The density at some radius.
        r - The radius.
    OUTPUTS:
        dM_dr - The mass at some radius, in kg/m.
    """

    dM_dr = 4*np.pi*np.power(r, 2)*rho

    return dM_dr

def dT_dr(rho, T, r, L, M):
    """Calculates the temperature via Equation 2.
    INPUTS:
        rho - The density at some radius.
        T - The temperature at some radius.
        r - The radius.
        L - The luminosity at some radius.
        M - The mass at some radius.
    OUTPUTS:
        dT_dr - The temperature at some radius, in K/m.

    DEFUNCT.
    """

    A = 3*kappa(rho, T)*rho*L/(16*np.pi*s.a*s.c*np.power(T, 3)*np.power(r, 2))
    B = (1 - 1/s.gamma)*T*s.G*M*rho/(pressure(rho, T)*np.power(r, 2))

    dT_dr = -np.minimum(A, B)

    return dT_dr

def drho_dr(rho, T, r, L, M):
    """Calculates the density via Equation 2.
    INPUTS:
        rho - The density at some radius.
        T - The temperature at some radius.
        r - The radius.
        L - The luminosity at some radius.
        M - The mass at some radius.
    OUTPUTS:
        drho_dr - The density at some radius, in kg/m^4.

    DEFUNCT.
    """

    A = s.G*M*rho/np.power(r, 2)
    B = dP_dT(rho, T)*dT_dr(rho, T, r, L, M)
    C = dP_drho(rho, T)

    drho_dr = -(A+B)/C

    return drho_dr

def luminosity_check(R, L, T):
    """This function is used as a check to see how close the ODE solved luminosity is to the
    theoretical value, taken from Equation 17 in the project description.
    INPUTS:
        R - The radius of the star.
        L - The luminosity of the star.
        T - The surface temperature of the star.

    OUTPUTS:
        f - If this value is 0 then we have a perfect match.
    """

    theoretical =  4*np.pi*s.sb*(R**2)*(T**4)
    f = (L-theoretical) / np.sqrt(L*theoretical)

    return f


# ---------------------------------
# Modification of Gravity Section:
# ---------------------------------

def dT_dr_scaled(rho, T, r, L, M, lam):
    """Calculates the temperature via Equation 2, with a gravitational modification according to Equation 19.
    INPUTS:
        rho - The density at some radius.
        T - The temperature at some radius.
        r - The radius.
        L - The luminosity at some radius.
        M - The mass at some radius.
    OUTPUTS:
        dT_dr - The temperature at some radius, in K/m.
    """

    A = 3*kappa(rho, T)*rho*L/(16*np.pi*s.a*s.c*np.power(T, 3)*np.power(r, 2))
    B = (1 - 1/s.gamma)*T*s.G*M*rho*(1+lam/r)/(pressure(rho, T)*np.power(r, 2))

    dT_dr = -np.minimum(A, B)

    return dT_dr

def drho_dr_scaled(rho, T, r, L, M, lam):
    """Calculates the density via Equation 2, with a gravitational modification according to Equation 19.
    INPUTS:
        rho - The density at some radius.
        T - The temperature at some radius.
        r - The radius.
        L - The luminosity at some radius.
        M - The mass at some radius.
    OUTPUTS:
        drho_dr - The density at some radius, in kg/m^4.
    """

    A = s.G*M*rho*(1+lam/r)/np.power(r, 2)
    B = dP_dT(rho, T)*dT_dr(rho, T, r, L, M)
    C = dP_drho(rho, T)

    drho_dr = -(A+B)/C

    return drho_dr
