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
    '''Determining the value of epislon using equations 8 and 9'''

    rho_5 = rho/1e5
    T_6 = T/1e6

    e_pp = (1.07e-7)*rho_5*np.power(s.X,2)*np.power(T_6,4) # output vale is units of W/kg
    e_cno = (8.24e-26)*rho_5*s.X*s.X_CNO*( (T_6)**(19.9) ) # output vale is units of W/kg

    e = e_pp + e_cno

    return (e)

def kappa_es():
    '''Calcualte and return the value for kappa_es in m**2/kg'''
    return 0.02*(1 + s.X)

def kappa_ff(rho, T):
    '''Calcualte and return the value for kappa_ff in m**2/kg'''
    rho_3 = rho/1e3

    return ( 1e24*(s.Z + 0.0001)*( (rho_3)**(0.7) )*( (T)**(-3.5) ) )

def kappa_Hminus(rho, T):
    '''Calcualte and return the value for kappa_Hminus in m**2/kg'''
    rho_3 = rho/1e3

    return ( (2.5e-32)*(s.Z / 0.02)*( (rho_3)**(0.5) )*( (T)**(9) ) )

def kappa(rho, T):
    '''Determing the value for kappa using equation 10,11,12,13 and 14'''

    kappa_es_value = kappa_es()
    kappa_ff_value = kappa_ff(rho, T) # Z would need to be a globally defined variable # units of m^2/kg
    kappa_Hminus_value = kappa_Hminus(rho,T)

    # splitting the term into two
    A = 1/kappa_Hminus_value
    B = 1/max(kappa_es_value,kappa_ff_value)

    return ( np.power(A+B,-1) )

def P_degen(rho):
    '''Calculate the degenerate pressure'''
    return ( ((3*np.pi**2)**(2/3)/5)*(s.hbar**2/s.me)*(rho/s.mp)**(5/3))

def P_ideal(rho,T):
    '''Calculate the ideal gas pressure'''
    return (rho*s.k*T/(s.mu*s.mp))

def P_photongas(T):
    '''Calculate the photon gas pressure'''
    return ((s.a*T**4)/3)

def pressure(rho, T):
    '''Calcualting the pressure as described by equation 5 and 6 '''

    P_deg = P_degen(rho)
    P_id = P_ideal(rho,T)
    P_pg = P_photongas(T)

    return P_deg + P_id + P_pg

def dP_dT(rho, T):
    '''This function will calcualte the partial derivate of pressure with respect to temperature (see equation 5)'''
    # Taking the partial derivative wrt temperature of eqution 5
    A = rho*s.k / (s.mu*s.mp)
    B = (4/3)*s.a*T**3

    return A + B

def dP_drho(rho, T):
    '''This function will calcualte the partial derivate of pressure with respect to density (see equation 5)'''
    # Taking the partial derivative wrt density of eqution 5
    A = ( (3*np.pi**2)**(2/3) / 3 )*( s.hbar**2 / (s.me*s.mp) )
    B = ( rho / s.mp)**(2/3)
    C = s.k*T/ (s.mu*s.mp)

    return ( (A*B) + C )

def dlnP_dlnT(P, T):
    '''This function will calculate the partial derivative of d(ln(P))/d(ln(T))'''
    return np.log(P)/np.log(T)


def dtau_dr(rho, T):
    '''This function will calculate the optical depth DE'''
    '''Calls the kappa function for the opacity'''

    return kappa(rho, T)*rho

def dL_dr_PP(rho,T,r):
    '''This function will calcualte dL/dr using only e_pp'''
    rho_5 = rho/1e5
    T_6 = T/1e6

    e_pp = (1.07e-7)*rho_5*np.power(s.X,2)*np.power(T_6,4) # output vale is units of W/kg

    return ( 4*np.pi*(r**2)*rho*e_pp )

def dL_dr_CNO(rho,T,r):
    '''This function will calcualte dL/dr using only e_cno'''
    rho_5 = rho/1e5
    T_6 = T/1e6

    e_cno = (8.24e-26)*rho_5*s.X*s.X_CNO*( (T_6)**(19.9) ) # output vale is units of W/kg

    return ( 4*np.pi*(r**2)*rho*e_cno )

def dL_dr(rho, T, r):
    '''This function will calcualte the luminosity DE'''

    return 4*np.pi*(r**2)*rho*epsilon(rho, T)

def dM_dr(rho, r):
    '''This function will calculate the mass DE'''
    return 4*np.pi*(r**2)*rho

def dT_dr(rho, T, r, L, M):
    '''This function will calcualte the termperature DE'''

    A = 3*kappa(rho, T)*rho*L/(16*np.pi*(s.c*s.a)*(T**3)*(r**2))
    B = (1 - 1/s.gamma)*(T/pressure(rho,T))*(s.G*M*rho/(r**2))

    min_value = -np.minimum(A,B)

    return min_value

def drho_dr(rho, T, r, L, M):
    '''This fucntion will calcualte the density DE'''

    A = ( (s.G*M*rho) / (r**2) )
    B = dP_dT(rho, T)*dT_dr(rho, T, r, L, M)
    C = dP_drho(rho, T)

    return -(A+B)/C

def luminosity_check(R, L, T):
    '''This function will check L(R_star) with the theoretical luminosity
    See equation 17 in the project description
    '''
    theoretical = 4*np.pi*s.sb*(R**2)*(T**4)
    f = (L-theoretical) / np.sqrt(L*theoretical)

    return f


# ---------------------------------
# Modification of Gravity Section:
# ---------------------------------

def dT_dr_scaled(rho, T, r, L, M, lam): # Built upon dT_dr function
    '''This function will calculate the termperature DE with scaled lambda factor for gravity
        Extrapolation of dT_dr function
    '''

    A = 3*kappa(rho, T)*rho*L/(16*np.pi*(s.sb*4)*(T**3)*(r**2))
    B = (1-1/s.gamma)*(T/pressure(rho,T))*(( ((s.G * M) / (r ** 2)) + ( (s.G * M * lam) / (r ** 3) ) )*rho) # y is replaced with Mass
    min_value = -np.minimum(A,B)

    return (min_value)

def drho_dr_scaled(rho, T, r, L, M, lam):
    '''This fucntion will calculate the density DE with scaled lambda factor for gravity
        Extrapolation of drho_dr function
    '''

    A = (( ((s.G * M) / (r ** 2)) + ( (s.G * M * lam) / (r ** 3) ) )*rho)
    B = (dP_dT(rho, T) )*(dT_dr_scaled(rho, T, r, L, M, lam))
    C = dP_drho(rho, T)

    return ( (-(A+B))/C )
