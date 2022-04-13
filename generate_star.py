#script: generate_star.py
#purpose: creates a star!

import os
import numpy as np
import pandas as pd
from standards import Standards as s
import tools as tools
import scipy.integrate as int
import matplotlib.pyplot as plt

import FofLambda as grav

np.seterr(all='ignore')

def solveODEs(r_init, r_fin, y0):
    '''Function to solve the ODEs (equation 2 in the project description)'''

    r_values = np.linspace(r_init, r_fin, 10000)
    t_span = (r_init, r_fin)

    solutions = int.solve_ivp(grav.dydr_gscaled, t_span, y0, 'RK45', t_eval=r_values) # outputs an array 'y' which contain the solutions to the ODEs which are described dydr from tools.py

    rho = solutions.y[0]
    T = solutions.y[1]
    M = solutions.y[2]
    L = solutions.y[3]
    tau = solutions.y[4]
    r = solutions.t # t here mean the 'time points' which in this case are the radius values or rvals

    return r, rho, T, M, L, tau

def bisection(rhoc_min, rhoc_max, r_initial, r_final, y, steps):
    '''This will use the bisection method to vary and determine the value for
    for rho_c such that the lumnosity boundary condition is  met
    '''

    rhoc_mid = (rhoc_min + rhoc_max)/2 # determine mid point for rhoc given the min and max inputs

    y_min = np.copy(y) # new intial value array for rhoc_min
    y_min[0] = rhoc_min

    y_mid = np.copy(y) # new intial value array for rhoc_mid
    y_mid[0] = rhoc_mid

    y_max = np.copy(y) # new intial value array for rhoc_max
    y_max[0] = rhoc_max

    # Now we want to use rhoc min,mid, and max to solve the ODEs individually
    solutions_min = solveODEs(r_initial, r_final, y_min, steps)
    solutions_mid = solveODEs(r_initial, r_final, y_mid, steps)
    solutions_max = solveODEs(r_initial, r_final, y_max, steps)

    # Need to find R_star, L_star, and T_star for the three solutions
    tau_min = solutions_min[5]
    tau_mid = solutions_mid[5]
    tau_max = solutions_max[5]

    i_min = np.argmin(tau_min[-1] - tau_min - 2/3)
    i_mid = np.argmin(tau_mid[-1] - tau_mid - 2/3)
    i_max = np.argmin(tau_max[-1] - tau_max - 2/3)

    R_Star_min = solutions_min[0][i_min]
    R_Star_mid = solutions_mid[0][i_mid]
    R_Star_max = solutions_max[0][i_max]

    L_Star_min = solutions_min[4][i_min]
    L_Star_mid = solutions_mid[4][i_mid]
    L_Star_max = solutions_max[4][i_max]

    T_Star_min = solutions_min[2][i_min]
    T_Star_mid = solutions_mid[2][i_mid]
    T_Star_max = solutions_max[2][i_max]

    # Using the luminosity_check (equation 17) for the valoues found using rhoc min,max, mid
    fmin = tools.luminosity_check(R_Star_min, L_Star_min, T_Star_min)
    fmid = tools.luminosity_check(R_Star_mid, L_Star_mid, T_Star_mid)
    fmax = tools.luminosity_check(R_Star_max, L_Star_max, T_Star_max)

    # Note: if 'f' is > 0 then our L(R_star) is greater than the theoretical

    if fmin*fmax > 0:
        print ("Change initial conditions, no rho_c possible.")
        quit()
    else:
        #initial condition for rhoc
        current = rhoc_mid
        #arbitrary initializing
        previous, diff = 0, 2

        while abs(diff) > 1:
            if fmin*fmid < 0: # Either fmin or fmid overshot while the other undershot
                rhoc_max = rhoc_mid

                rhoc_mid = (rhoc_min+rhoc_max)/2

                # Now we re-solve ODEs using the new rhoc_mid
                y_mid = np.copy(y) # new intial value array for rhoc_mid
                y_mid[0] = rhoc_mid

                soln = solveODEs(r_initial, r_final, y_mid, steps)
                tau = soln[5]
                i = np.argmin(tau[-1] - tau - 2/3)


                R_Star_mid = soln[0][i]
                T_Star_mid = soln[2][i]
                L_Star_mid = soln[4][i]

                fmid = tools.luminosity_check(R_Star_mid, L_Star_mid, T_Star_mid)

            elif fmid*fmax < 0: # Either fmid or fmax overshot while the other undershot
                rhoc_min = rhoc_mid
                rhoc_mid = (rhoc_min+rhoc_max)/2

                # Now we re-solve ODEs using the new rhoc_mid
                y_mid =  np.copy(y) # new intial value array for rhoc_mid
                y_mid[0] = rhoc_mid

                soln = solveODEs(r_initial, r_final, y_mid, steps)
                tau = soln[5]
                i = np.argmin(tau[-1] - tau - 2/3)


                R_Star_mid = soln[0][i]
                T_Star_mid = soln[2][i]
                L_Star_mid = soln[4][i]

                fmid = tools.luminosity_check(R_Star_mid, L_Star_mid, T_Star_mid)

            elif fmid == 0:
                return rhoc_mid

            previous = current
            current = rhoc_mid
            diff = current - previous

    return rhoc_mid

def get_f(rhoc, Tc):

    #the range of radii, in kg/m^3
    r0, r_max = 1e-6, 50*s.Rsun
    rs = np.linspace(r0, r_max, 10000)

    #Step 1: Select rhoc and Tc
    #Step 2: Integrate equations (2) to obtain Rstar, L(Rstar), T(Rstar)
    M = 4*np.pi*(r0**3)*rhoc/3
    L0 = 4*np.pi*(r0**3)*rhoc*tools.epsilon(rhoc, Tc)
    tau0 = tools.kappa(rhoc, Tc)*rhoc*r0

    ic = np.array([rhoc, Tc, M, L0, tau0])

    i = 0
    deltau = 1
    tolerance = 1e-3
    #r, rho, T, M, L, tau = solveODEs(r0, r_max, initial_conditions)
    while (deltau > tolerance) and (M < 1e3 * s.Msun):

        ic = tools.RKF45Method_adaptive(grav.dydr_gscaled, rs[i], ic, 1000)
        rho, T, M, L, tau = ic
        deltau = tools.del_tau(rho,T,rs[i],M,L)

        if i%1000 == 0:
            print(deltau)
            print(ic)

        i+=1
    #Equation 4 in the project description says we should find when this is zero to define the surface
    #If there was no np.abs then it would gravitate towards the final index
    i = np.argmin(np.abs(tau[-1] - tau - 2/3))

    #Get values at the surface
    Rstar, Tstar, Mstar, Lstar = r[i], T[i], M[i], L[i]

    #Following Equation 17 in the project description
    theoretical = 4*np.pi*s.sb*(Rstar**2)*(Tstar**4)
    f_rhoc = (Lstar-theoretical) / ((Lstar*theoretical)**(1/2))

    return f_rhoc, r, rho, T, M, L, tau, i

def new_bisection(Tc):
    """A bisection method to find the central density to use to solve Equations 2.
    Source: https://personal.math.ubc.ca/~pwalls/math-python/roots-optimization/bisection/
    INPUTS:
        Tc - The central temperature.
    OUTPUTS:
        mid_rhoc - The central density solved for.
    """

    #starting range of densities, in kg/m^3
    min_rhoc, max_rhoc = 300, 500000

    fmid, mid_rhoc = 10, 10 #arbitrary starting point / irrelevant
    tolerance = 1e-5 #tolerance for the bisection method

    count = 0 #for the sake of processing time
    while (abs(fmid) > tolerance) and (count < 200):

        mid_rhoc = (max_rhoc + min_rhoc)/2

        fmin = get_f(min_rhoc, Tc)[0]
        fmid = get_f(mid_rhoc, Tc)[0]
        fmax = get_f(max_rhoc, Tc)[0]

        #print(fmin, fmid, fmax, mid_rhoc)

        #nothing will change if this happens
        if fmin == fmid or fmid == fmax:
            break

        if fmin*fmax > 0:
            print('No roots are possible, function does not cross zero.')
            mid_rhoc = None
            break

        if abs(fmid) < tolerance: break #close enough to zero, we got it

        #if the root is in-between the minimum and the mid-point
        if fmin*fmid < 0:
            max_rhoc = mid_rhoc

        #if the root is in-between the mid-point and the maximum
        elif fmid*fmax < 0:
            min_rhoc = mid_rhoc

        count +=1


    return mid_rhoc

def create_star(Tc):
    """Creates a star given some input central temperature.
    """

    rhoc = new_bisection(Tc)
    _, r, rho, T, M, L, tau, i = get_f(rhoc, Tc)
    R_star, rho_c, M_star, L_star, T_star = r[i], rho[0], M[i], L[i], T[i]

    #get the values up to the radius
    r, rho_r, T_r, M_r, L_r = r[:i], rho[:i], T[:i], M[:i], L[:i]

    #Calculate the pressures
    P_r = tools.pressure(rho_r, T_r)
    Pdeg = tools.P_degen(rho_r)
    Pgas = tools.P_ideal(rho_r, T_r)
    Ppho = tools.P_photongas(T_r)

    #Calculating the mean opacities
    kappa_values = []
    for i in range(0, len(rho_r)):
        kappa_values.append(tools.kappa(rho_r[i], T_r[i]))

    kappa_ff_values = tools.kappa_ff(rho_r, T_r)
    kappa_es_values = np.full(len(kappa_ff_values), tools.kappa_es())
    kappa_Hminus_values = tools.kappa_Hminus(rho_r, T_r)

    # Saving all values related to dL/dr
    dL_dr_values = tools.dL_dr(rho_r , T_r, r)
    dL_dr_pp_values = tools.dL_dr_PP(rho_r , T_r, r)
    dL_dr_cno_values = tools.dL_dr_CNO(rho_r , T_r, r)

    # Saving all values for dln(P)/dln(T) partial derivative
    dlnP_dlnT_values = tools.dlnP_dlnT(P_r, T_r)

    df = pd.DataFrame(data=zip(r/R_star, rho_r, L_r, M_r, T_r, Pgas,
                               Pdeg, Ppho, P_r, dlnP_dlnT_values,
                               kappa_values, kappa_ff_values, kappa_Hminus_values, kappa_es_values,
                               dL_dr_values, dL_dr_pp_values, dL_dr_cno_values),
                      columns=['r/R', 'rho(r)', 'L(r)', 'M(r)', 'T(r)', 'Pgas', 'Pdeg', 'Ppho', 'P', 'dlogP/dlogT',
                               'kappa', 'kappaff', 'kappaHm', 'kappaes', 'dL/dr', 'dLpp/dr', 'dLcno/dr'])

    #to g/cm^3
    df['rho(r)'] = df['rho(r)']/1000
    #to g
    df['M(r)'] = df['M(r)']*1000
    #to erg/s
    df['L(r)'] = df['L(r)'] / 1e-7

    #to erg/cm^3
    df['P'] = df['P'] * 0.1
    df['Pdeg'] = df['Pdeg'] * 0.1
    df['Pgas'] = df['Pgas'] * 0.1
    df['Ppho'] = df['Ppho'] * 0.1

    #to erg/s/cm
    df['dL/dr'] = df['dL/dr'] / 1e-7 / 100
    df['dLpp/dr'] = df['dLpp/dr'] / 1e-7 / 100
    df['dLcno/dr'] = df['dLcno/dr'] / 1e-7 / 100

    #to cm^2/g
    df['kappa'] = df['kappa'] * 10
    df['kappaff'] = df['kappaff'] * 10
    df['kappaes'] = df['kappaes'] * 10
    df['kappaHm'] = df['kappaHm'] * 10

    return df, (rho_c/1000, Tc, R_star*100, M_star*1000, L_star/1e-7, T_star)


    """
    #set the initial conditions
    r_initial = 1
    r_final = 1e11
    steps = 100000 #number of steps

    # intializing an array that contains the intial values for density, temperature, mass, luminosity and tau
    y = np.array([0, Tc, 0, 0, 0])

    #our solution to f(rho_c)
    rho_c = bisection(rhoc_min, rhoc_max, r_initial, r_final, y, steps)

    # re-defining the previos 'y' array to now contain the newly determine rho_c value
    y_new =  np.array([rho_c, Tc, 0, 0, 0])

    # Solving the ODEs
    final = solveODEs(r_initial, r_final, y_new, steps)

    # Storing the solutions to the ODEs appropriately
    r, rho_r, T_r, M_r, L_r, tau_r = final

    i = np.argmin(tau_r[-1] - tau_r - 2/3)

    """
