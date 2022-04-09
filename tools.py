#script: tools.py
#purpose: To have a hub for all useful functions.

import os
import numpy as np
import pandas as pd
from standards import Standards as s

def load_data():
    """Loads in the data. Warning, the data is in CGMS units, NOT SI units.
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
       matrix representing [rho,Temp,Mass,Luminosity,OpticalDepth, lambda]
    """
    rho = y[0]
    T = y[1]
    M = y[2]
    L = y[3]
    tau = y[4]

    # The five ODEs we are trying to solve (equation 2)
    dydr = np.zeros(5)
    dydr[0] = drho_dr(rho, T, r, M,L)
    dydr[1] = dT_dr(rho, T, r, L, M)
    dydr[2] = dM_dr(rho, r)
    dydr[3] = dL_dr(rho , T, r)
    dydr[4] = dtau_dr(rho, T)

    return dydr

def RKF45Method_adaptive(f,r,y,h):
    '''Second RK4Method which invovled adpative step sizes
       I followed the algorithm outlined here -> https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method#Implementing_an_RK4(5)_Algorithm
       The coefficietns for 'A' were taken from the table title "COEFFICIENTS FOR RK4(5), FORMULA 2 Table III in Fehlberg"

       Input(s):
       - 'r' - array of radius
       - 'y' - 5 column matrix representing [rho,Temp,Mass,Luminosity,OpticalDepth]
       - 'h' - the stepsize

       Output(s):
       - 'y5' - 5th order solutions to ODE
       - 'hnew' - the new adpted stepsize
    '''
    K1 = h*f(r,y)
    K2 = h*f( r + (1/4)*h, y + (1/4)*K1 )
    K3 = h*f( r + (3/8)*h, y + (3/32)*K1 + (9/32)*K2 )
    K4 = h*f( r + (12/13)*h, y + (1932/2197)*K1 + (-7200/2197)*K2 + (7296/2197)*K3 )
    K5 = h*f( r + (1)*h, y + (439/216)*K1 + (-8)*K2 + (3680/513)*K3 + (-845/4104)*K4 )
    K6 = h*f( r + (1/2)*h, y + (-8/27)*K1 + (2)*K2 + (-3544/2565)*K3 + (1859/4104)*K4 + (-11/40)*K5 )

    # The truncation error (TE)
    TE = abs( (1/360)*K1 + (0)*K2 + (-128/4275)*K3 + (-2197/75240)*K4 + (1/50)*K5 + (2/55)*K6 )

    # the error estimate of the RKF45 method
    y5 = y + (25/216)*K1 + (1408/2565)*K3 + (2197/4101)*K4 - (1/5)*K5
    y4 = y + (16/135)*K1 + (6656/12825)*K3 + (28561/56430)*K4 - (9/50)*K5 + (2/55)*K6
    error = abs( y5 - y4)

    # Calcuating the new step-sizes
    hnew = 0.9*h*( (min(error)/TE)**(1/5) )
    return(y5,hnew)


def RK4Method(ri,rf,y0,h):
    """Performing the RK4 algorithm.
       y0 are the initial conditions for rho, temp, mass, luminosity and  tau (in that order)
    """

    #Setting up r and y
    steps = int((rf-ri)/h)
    r = np.linspace(ri,rf,steps)
    y = np.zeros( (5, len(r)) )

    #initial conditions

    y[0,0] = y0[0]   # intial rho
    y[0,1] = y0[1]  # inital temp
    y[0,2] = y0[2] # inital mass - centre boundary condition
    y[0,3] = y0[3]  # intial luminosity?- centre boundary condition
    y[0,4] = y0[4]  # inital optical depth (tau)



    #for loop to calculate y matrix using RK4 method. Should output rho,Temp,Mass,Luminosity,OpticalDepth at all radii
    for i in range(r.size - 1):

        k0 = h*dydr(r[i],y[i,:])[0]
        k1 = h*dydr(r[i],y[i,:])[1]
        k2 = h*dydr(r[i],y[i,:])[2]
        k3 = h*dydr(r[i],y[i,:])[3]
        k4 = h*dydr(r[i],y[i,:])[4]

        s0 = h*dydr(r[i]+h/2,[y[i,0]+k0/2,y[i,1]+k1/2,y[i,2]+k2/2,y[i,3]+k3/2,y[i,4]+k4/2])[0]
        s1 = h*dydr(r[i]+h/2,[y[i,0]+k0/2,y[i,1]+k1/2,y[i,2]+k2/2,y[i,3]+k3/2,y[i,4]+k4/2])[1]
        s2 = h*dydr(r[i]+h/2,[y[i,0]+k0/2,y[i,1]+k1/2,y[i,2]+k2/2,y[i,3]+k3/2,y[i,4]+k4/2])[2]
        s3 = h*dydr(r[i]+h/2,[y[i,0]+k0/2,y[i,1]+k1/2,y[i,2]+k2/2,y[i,3]+k3/2,y[i,4]+k4/2])[3]
        s4 = h*dydr(r[i]+h/2,[y[i,0]+k0/2,y[i,1]+k1/2,y[i,2]+k2/2,y[i,3]+k3/2,y[i,4]+k4/2])[4]

        l0 = h*dydr(r[i]+h/2,[y[i,0]+s0/2,y[i,1]+s1/2,y[i,2]+s2/2,y[i,3]+s3/2,y[i,4]+s4/2])[0]
        l1 = h*dydr(r[i]+h/2,[y[i,0]+s0/2,y[i,1]+s1/2,y[i,2]+s2/2,y[i,3]+s3/2,y[i,4]+s4/2])[1]
        l2 = h*dydr(r[i]+h/2,[y[i,0]+s0/2,y[i,1]+s1/2,y[i,2]+s2/2,y[i,3]+s3/2,y[i,4]+s4/2])[2]
        l3 = h*dydr(r[i]+h/2,[y[i,0]+s0/2,y[i,1]+s1/2,y[i,2]+s2/2,y[i,3]+s3/2,y[i,4]+s4/2])[3]
        l4 = h*dydr(r[i]+h/2,[y[i,0]+s0/2,y[i,1]+s1/2,y[i,2]+s2/2,y[i,3]+s3/2,y[i,4]+s4/2])[4]

        p0 = h*dydr(r[i]+h,[y[i,0]+l0,y[i,1]+l1,y[i,2]+l2,y[i,3]+l3,y[i,4]+l3])[0]
        p1 = h*dydr(r[i]+h,[y[i,0]+l0,y[i,1]+l1,y[i,2]+l2,y[i,3]+l3,y[i,4]+l3])[1]
        p2 = h*dydr(r[i]+h,[y[i,0]+l0,y[i,1]+l1,y[i,2]+l2,y[i,3]+l3,y[i,4]+l3])[2]
        p3 = h*dydr(r[i]+h,[y[i,0]+l0,y[i,1]+l1,y[i,2]+l2,y[i,3]+l3,y[i,4]+l3])[3]
        p4 = h*dydr(r[i]+h,[y[i,0]+l0,y[i,1]+l1,y[i,2]+l2,y[i,3]+l3,y[i,4]+l3])[4]

        y[i+1,0] = y[i+1,0] + 1/6*(k0+2*s0+2*l0+p0)
        y[i+1,2] = y[i+1,1] + 1/6*(k1+2*s1+2*l1+p1)
        y[i+1,3] = y[i+1,2] + 1/6*(k2+2*s2+2*l2+p2)
        y[i+1,4] = y[i+1,3] + 1/6*(k3+2*s3+2*l3+p3)
        y[i+1,5] = y[i+1,4] + 1/6*(k4+2*s4+2*l4+p4)

        return y

def epsilon(rho, T):
    '''Determining the value of epislon using equations 8 and 9'''

    X_CNO = 0.03*s.X # X would need to be a globally defined variable
    rho_5 = rho/1e5
    T_6 = T/1e6

    e_pp = (1.07e-7)*rho_5*np.power(s.X,2)*np.power(T_6,4) # output vale is units of W/kg
    e_cno = (8.24e-26)*rho_5*s.X*X_CNO*( (T_6)**(19.9) ) # output vale is units of W/kg

    e = e_pp + e_cno

    return (e)

def kappa_es():
    '''Calcualte and return the value for kappa_es in m**2/kg'''
    return (0.028*(1 + s.X))

def kappa_ff(rho, T):
    '''Calcualte and return the value for kappa_ff in m**2/kg'''
    rho_3 = rho/1e3

    return ( (1.0e24)*(s.Z + 0.0001)*( (rho_3)**(0.7) )*( (T)**(-3.5) ) )

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
    return ( ( (3*np.pi**2)**(2/3) / 5 )*( s.hbar**2 / s.me )*( rho / s.mp)**(5/3) )

def P_ideal(rho,T):
    '''Calculate the ideal gas pressure'''
    return (rho*s.k*T/s.mu*s.mp)

def P_photongas(T):
    '''Calculate the photon gas pressure'''
    return ((s.a*T**4)/3)

def pressure(rho, T):
    '''Calcualting the pressure as described by equation 5 and 6 '''

    P_deg = P_degen(rho)
    P_id = P_ideal(rho,T)
    P_pg = P_photongas(T)

    return ( P_deg + P_id + P_pg )

def dP_dT(rho, T):
    '''This function will calcualte the partial derivate of pressure with respect to temperature (see equation 5)'''
    # Taking the partial derivative wrt temperature of eqution 5
    A = rho*s.k / (s.mu*s.mp)
    B = (4/3)*s.a*T**3

    return ( A+B )

def dP_drho(rho, T):
    '''This function will calcualte the partial derivate of pressure with respect to density (see equation 5)'''
    # Taking the partial derivative wrt density of eqution 5
    A = ( (3*np.pi**2)**(2/3) / 3 )*( s.hbar**2 / (s.me*s.mp) )
    B = ( rho / s.mp)**(2/3)
    C = s.k*T/ (s.mu*s.mp)

    return ( (A*B) + C )

def dlnP_dlnT(rho, T):
    '''This function will calculat the partial derivative of d(ln(P))/d(ln(T))'''
    C1 = rho*s.k / (s.mu*s.mp)
    C2 = (1/3)*s.a*T

    numerator = C1*T + 4*C2*T**4
    denominator = pressure(rho, T)
    return ( numerator/denominator )

def dtau_dr(rho, T):
    '''This function will calculate the optical depth DE'''
    '''Calls the kappa function for the opacity'''

    return ( kappa(rho, T)*rho )

def dL_dr_PP(rho,T,r):
    '''This function will calcualte dL/dr using only e_pp'''
    rho_5 = rho/1e5
    T_6 = T/1e6

    e_pp = (1.07e-7)*rho_5*np.power(s.X,2)*np.power(T_6,4) # output vale is units of W/kg

    return ( 4*np.pi*(r**2)*rho*e_pp )

def dL_dr_CNO(rho,T,r):
    '''This function will calcualte dL/dr using only e_cno'''
    X_CNO = 0.03*s.X # X would need to be a globally defined variable
    rho_5 = rho/1e5
    T_6 = T/1e6

    e_cno = (8.24e-26)*rho_5*s.X*X_CNO*( (T_6)**(19.9) ) # output vale is units of W/kg

    return ( 4*np.pi*(r**2)*rho*e_cno )

def dL_dr(rho , T, r):
    '''This function will calcualte the luminosity DE'''

    return ( 4*np.pi*(r**2)*rho*epsilon(rho, T) )

def dM_dr(rho, r):
    '''This function will calculate the mass DE'''
    return ( 4*np.pi*(r**2)*rho )

def dT_dr(rho, T, r, L, M):
    '''This function will calcualte the termperature DE'''

    A = 3*kappa(rho, T)*rho*L/(16*np.pi*(s.sb*4)*(T**3)*(r**2))
    B = (1-1/s.gamma)*(T/pressure(rho,T))*(s.G*M*rho/(r**2))

    min_value = -np.minimum(A,B)

    return (min_value)

def drho_dr(rho, T, r, M,L):
    '''This fucntion will calcualte the density DE'''

    A = ( (s.G*M*rho) / (r**2) )
    B = (dP_dT(rho, T) )*(dT_dr(rho, T, r, L, M))
    C = dP_drho(rho, T)

    return ( (-(A+B))/C )

def del_tau(rho,T,r,M,L):
    '''This fucntion will calcualte the opacity proxy'''
    return kappa(rho,T)*(rho**2) / (abs(drho_dr(rho, T, r, M,L)))

def luminosity_check(r, L, T):
    '''This function will check L(R_star) with the theoretical luminosity
    See equation 17 in the project description
    '''
    theoretical = 4*np.pi*s.sb*(r**2)*(T**4)
    f = (L-theoretical) / np.sqrt(L*theoretical)

    return(f)


# ---------------------------------
# Modification of Gravity Section:
# ---------------------------------

# Lambda Guess': global variables
lam_guess = 1000 # Guess for constraints of lambda, from negative int to positive int
lam_lim_guess = 2 # Guess for scaling factor, limit is the length of lambda_values / lambda_limit_guess

# Function to scale g
def g_scaled(r,y,lam_guess,lam_lim_guess): # Based on equation 19
    '''
        Inputs (given as a function call of RK4):
            r - Array (Radius),
            y - Array (for y[2]: Mass)
            lam_guess - Int (Global Variable)
            lam_lim_guess - Int (Global Variable)

        Returns:
            g scaling as r^-3 below limit
            g non-scaled above chosen limit
    '''

    # Lambda constraints
    lambda_vals = np.linspace(-lam_guess, lam_guess, (len(r))) # values of lambda ranging from guesses
    limit = (len(lambda_vals)) / lam_lim_guess # limit arbitrarily chosen, can change according to limit consideration...?

    # Container for final values of g (scaled or unscaled)
    final_g = []
    n = 0 # Counter for iteration loop

    # Logic for lambda limit:
    while n <= (len(r) - 1):
        for i in (lambda_vals[0:limit]):
            final_g.append( ((s.G * y[2][i]) / (r[i] ** 2)) + ((s.G * y[2][i]) * lambda_vals[i] / (r[i] ** 3)) )
        for j in (lambda_vals[limit:]):
            final_g.append(((s.G * y[2][j]) / (r[j] ** 2)))

        n += 1

    return final_g

# Edit drho_dr and dTdr from above to account for gravity scaling

def dT_dr_scaled(rho, T, r, L, M): # Built upon dT_dr function
    '''This function will calculate the termperature DE with scaled lambda factor for gravity
        Extrapolation of dT_dr function
    '''

    A = 3*kappa(rho, T)*rho*L/(16*np.pi*(s.sb*4)*(T**3)*(r**2))
    B = (1-1/s.gamma)*(T/pressure(rho,T))*((g_scaled(r,M,lam_guess,lam_lim_guess)*rho)) # y is replaced with Mass
    min_value = -np.minimum(A,B)

    return (min_value)

def drho_dr_scaled(rho, T, r, M,L):
    '''This fucntion will calculate the density DE with scaled lambda factor for gravity
        Extrapolation of drho_dr function
    '''

    A = ((g_scaled(r,M,lam_guess,lam_lim_guess)*rho))
    B = (dP_dT(rho, T) )*(dT_dr_scaled(rho, T, r, L, M))
    C = dP_drho(rho, T)

    return ( (-(A+B))/C )

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
    lamb = y[5]

    # The five ODEs we are trying to solve (equation 2)
    dydr = np.zeros(5)
    dydr[0] = drho_dr_scaled(rho, T, r, M,L)
    dydr[1] = dT_dr_scaled(rho, T, r, L, M)
    dydr[2] = dM_dr(rho, r)
    dydr[3] = dL_dr(rho , T, r)
    dydr[4] = dtau_dr(rho, T)

    return dydr

# Finally, run RK4 method for new scaled values
