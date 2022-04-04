#script: tools.py
#purpose: To have a hub for all useful functions.

import os
import numpy as np
import pandas as pd
from standards import Standards as s

def load_data():
    """Loads in the data. Warning, the data is in CGMS units, NOT SI units.
    """

    highmass = pd.read_csv(os.path.join(s.data_folder, 'highmass_star.txt'), sep='\s+', header=0, skiprows=[1,2])
    lowmass = pd.read_csv(os.path.join(s.data_folder, 'lowmass_star.txt'), sep='\s+', header=0, skiprows=[1,2])
    summary = pd.read_csv(os.path.join(s.data_folder, 'star_summary.txt'), sep='\s+', header=0, skiprows=[1,2])

    return highmass, lowmass, summary

def dydr(r,y) :
    """Defining the 5 ODEs to solve using RK4. Takes radius(r) and y as inputs where y is a 5 column
       matrix representing [rho,Temp,Mass,Luminosity,OpticalDepth]
    """

    dydr = np.zeros(5)
    dydr[0] = -(G*y[2]*y[0]/r**2 + partialP/partialT*dydr[1])/partialP/partialrho #Determine PDEs used (Could be as simple as subbing eq 7 from project instructions)
    dydr[1] = -np.min( (3*kappa(y[0], y[1])*y[0]*y[3]/(16*np.pi*s.a*s.c*y[1]**3*r**2)), (1-1/s.gamma)*(y[1]*s.G*y[2]*y[0])/(P*r**2) ) # gamma is defined in standards.py. Someone please check whether I called the standards class correctly - DP
    dydr[2] = 4*np.pi*r**2*y[0] # mass differential equation
    dydr[3] = 4*np.pi*r**2*y[0]*epsilon(y[0], y[1]) # Added an epsilon fucntion which takes rho and temperature. From what I understood y[0] and y[1] should contain those values -DP
    dydr[4] = kappa(y[0], y[1])*y[0] # Added a kappa fucntion which takes rho and temperature. From what I understood y[0] and y[1] should contain those values -DP

    return dydr

def RK4Method(ri,rf,h):
    """Performing the RK4 algorithm.
    """

    #Setting up r and y
    steps = int(rf-ri)/h
    r = np.linspace(ri,rf,steps)
    y = np.zeros(len(r),5)

    #initial conditions
    y[0,0] =   # intial rho?
    y[0,1] =   # inital temp?
    y[0,2] =   # inital mass?
    y[0,3] =   # intial luminosity?
    y[0,4] =   # inital optical depth (tau)?



    #for loop to calculate y matrix using RK4 method. Should output rho,Temp,Mass,Luminosity,OpticalDepth at all radii
    for i in range(0,len(r)-1)

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

    X_CNO = 0.03*X # X would need to be a globally defined variable
    rho_5 = rho/1e5
    T_6 = T/1e6

    e_pp = (1.07e-7)*rho_5*np.power(X,2)*np.power(T_6,4) # output vale is units of W/kg
    e_cno = (8.24e-26)*rho_5*X*X_CNO*np.power(T_6, 19.9) # output vale is units of W/kg

    e = e_pp + e_cno

    return (e)

def kappa(rho, T):
    '''Determing the value for kappa using equation 10,11,12,13 and 14'''
    rho_3 = rho/1e3

    kappa_es = 0.028*(1 + X)  # units of m^2/kg
    kappa_ff = (1.0e24)*(Z + 0.0001)*np.power(rho_3,0.7)*np.power(T,-3.5) # Z would need to be a globally defined variable # units of m^2/kg
    kappa_Hminus = (2.5e-32)*(Z / 0.02)*np.power(rho_3,0.5)*np.power(T,9) # units of m^2/kg

    # splitting the term into two
    A = 1/kappa_Hminus
    B = 1/np.max(kappa_es,kappa_ff)

    return ( np.power(A+B,-1) )

def pressure(rho, T):
    '''Calcualting the pressure as described by equation 5 and 6 '''
    '''NOT DONE'''

    P_degenerate = (s.hbar**2)*(3*np.pi**2)**(2/3)
